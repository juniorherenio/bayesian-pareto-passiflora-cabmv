# =============================================================================
# Bayesian Bivariate Mixed Models — Genetic Correlations Between Traits
#
# OBJECTIVE:
#   Estimate rg (genetic correlation coefficient) between target trait pairs
#   utilizing a multi-trait bivariate MCMCglmm framework with an unstructured
#   covariance matrix specification: us(trait):animal
#
# BIOLOGICALLY RELEVANT TARGET PAIRS:
#   1. comp_mm       × diam_mm          (Fruit morphology dimensions)
#   2. comp_mm       × massa_fruto_g    (Fruit size × weight interaction)
#   3. esp_casca_mm  × rend_polpa_pct   (Pathosystem trade-off dynamics)
#   4. esp_casca_mm  × massa_polpa_g    (Rind thickness × pulp mass allocation)
#   5. massa_fruto_g × rend_polpa_pct   (Primary yield component interaction)
#   6. sst_brix      × rend_polpa_pct   (Fruit quality parameters link)
#   7. sst_brix      × esp_casca_mm     (Soluble solids × rind thickness relation)
#
# DISEASE COVARIATE: P_severo_media_sc (POM) — consistent with 04e_v2 script
# =============================================================================

library(here)
library(MCMCglmm)
library(Matrix)
library(coda)
library(dplyr)

cat("\n", strrep("=", 60), "\n")
cat("04f_bivariate.R — Genetic Correlations\n")
cat(format(Sys.time(), "[%H:%M:%S]"), "\n")
cat(strrep("=", 60), "\n\n")

# =============================================================================
# 1. SETUP AND INVERSE RELATIONSHIP MATRIX INITIALIZATION
# =============================================================================

dataset <- readRDS(here("01_dados/processados/dataset_pom.rds"))
Amat    <- readRDS(here("01_dados/processados/Amat.rds"))

dataset <- as.data.frame(dataset)
dataset$bloco     <- as.factor(dataset$bloco)
dataset$familia   <- as.factor(dataset$familia)
dataset$parcela   <- as.factor(dataset$parcela)
dataset$individuo <- factor(as.factor(dataset$individuo),
                            levels = rownames(Amat))

dataset$estacao_ano <- as.character(dataset$estacao_ano)
dataset$estacao_ano[dataset$estacao_ano == "verao_2018"] <- "verao_2019"
dataset$estacao_ano <- factor(dataset$estacao_ano,
                              levels = c("inverno_2018","primavera_2018",
                                         "verao_2019","outono_2019",
                                         "inverno_2019"))

dat <- dataset[!is.na(dataset$estacao_ano), ]
dat <- droplevels(dat)

dat$estacao <- dat$estacao_ano
levels(dat$estacao) <- c("i18","p18","v19","o19","i19")

dat$animal <- dat$individuo
dat$pe     <- dat$individuo
dat$P_severo_media_sc[is.na(dat$P_severo_media_sc)] <- 0

ind_obs     <- levels(droplevels(dat$individuo))
Amat_sub    <- Amat[ind_obs, ind_obs]
Ainv_sparse <- as(solve(Amat_sub), "dgCMatrix")
dat$individuo <- factor(dat$individuo, levels=rownames(Amat_sub))
dat$animal    <- dat$individuo
dat$pe        <- dat$individuo

cat(sprintf("Analytical Sample Profile: %d obs | %d genotypes available\n\n", nrow(dat), length(ind_obs)))

# Core multi-trait evaluation list matrix
pares <- list(
  list(t1="comp_mm",        t2="diam_mm",        nome="comp_diam"),
  list(t1="comp_mm",        t2="massa_fruto_g",  nome="comp_massa"),
  list(t1="esp_casca_mm",  t2="rend_polpa_pct",  nome="casca_rend"),
  list(t1="esp_casca_mm",  t2="massa_polpa_g",   nome="casca_polpa"),
  list(t1="massa_fruto_g", t2="rend_polpa_pct",  nome="massa_rend"),
  list(t1="sst_brix",      t2="rend_polpa_pct",  nome="brix_rend"),
  list(t1="sst_brix",      t2="esp_casca_mm",     nome="brix_casca")
)

# MCMC parameters definition
NITT   <- 150000; BURNIN <- 50000; THIN <- 100

dir.create(here("04_resultados/mcmc/chains"),
           showWarnings=FALSE, recursive=TRUE)
dir.create(here("04_resultados/parametros_geneticos/bayesiano"),
           showWarnings=FALSE, recursive=TRUE)

# =============================================================================
# 2. POSTERIOR GENETIC CORRELATION (rg) EXTRACTION ENGINE
# =============================================================================

extrair_rg <- function(m, t1, t2) {
  vc     <- as.data.frame(m$VCV)
  nm_v1  <- paste0("trait", t1, ":trait", t1, ".animal")
  nm_v2  <- paste0("trait", t2, ":trait", t2, ".animal")
  nm_cov <- paste0("trait", t1, ":trait", t2, ".animal")
  
  # Structural string order safeguard check
  if (!nm_cov %in% names(vc)) {
    nm_cov <- paste0("trait", t2, ":trait", t1, ".animal")
  }
  
  if (all(c(nm_v1, nm_v2, nm_cov) %in% names(vc))) {
    rg_post <- vc[[nm_cov]] / sqrt(vc[[nm_v1]] * vc[[nm_v2]])
    return(rg_post)
  }
  return(NULL)
}

# =============================================================================
# 3. MULTI-TRAIT BIVARIATE ITERATION ROUTINE
# =============================================================================

t0_total   <- proc.time()
resultados <- list()

for (par in pares) {
  t1   <- par$t1
  t2   <- par$t2
  nome <- par$nome
  
  cat(strrep("=", 50), "\n")
  cat(sprintf("EVALUATING PAIR: %s × %s\n\n", t1, t2))
  
  # Extraction baseline variances
  Vp1 <- var(dat[[t1]], na.rm=TRUE)
  Vp2 <- var(dat[[t2]], na.rm=TRUE)
  Vp  <- c(Vp1, Vp2)
  
  # Parameterize Inverse Wishart prior matrices for 'us(trait):animal' covariance structure
  prior_biv <- list(
    G = list(
      G1 = list(V = diag(2) * mean(Vp) * 0.2, nu = 3), # Additive genetic matrix prior
      G2 = list(V = diag(2) * mean(Vp) * 0.2, nu = 3), # Permanent environmental prior
      G3 = list(V = diag(2) * mean(Vp) * 0.05, nu = 3) # Plot spatial prior
    ),
    R = list(V = diag(2) * mean(Vp) * 0.4, nu = 3)     # Residual covariance prior matrix
  )
  
  # Multivariate fixed formula specification
  fixed_f <- as.formula(
    paste0("cbind(", t1, ",", t2, ") ~ trait - 1 +",
           " trait:bloco + trait:estacao + trait:P_severo_media_sc")
  )
  
  cat(sprintf("[%s] Running bivariate model estimation...\n",
              format(Sys.time(), "%H:%M:%S")))
  t0 <- proc.time()
  set.seed(2026)
  
  m <- tryCatch(
    MCMCglmm(
      fixed    = fixed_f,
      random   = ~ us(trait):animal + us(trait):pe + us(trait):parcela,
      rcov     = ~ us(trait):units,
      ginverse = list(animal = Ainv_sparse),
      prior    = prior_biv,
      data     = dat,
      nitt     = NITT,
      burnin   = BURNIN,
      thin     = THIN,
      family   = c("gaussian","gaussian"),
      verbose  = FALSE
    ),
    error = function(e) {
      cat("   ✗ Fitting Error:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (is.null(m)) next
  
  elapsed <- (proc.time() - t0)["elapsed"] / 60
  cat(sprintf("[%s] Model convergence reached — runtime: %.1f min | DIC=%.1f\n",
              format(Sys.time(), "%H:%M:%S"), elapsed, m$DIC))
  
  # Export MCMC chain checkpoint
  saveRDS(m, here(sprintf(
    "04_resultados/mcmc/chains/biv_%s.rds", nome)))
  
  # Process posterior correlation vector
  rg_post <- extrair_rg(m, t1, t2)
  
  if (!is.null(rg_post)) {
    rg_mean <- mean(rg_post, na.rm=TRUE)
    rg_025  <- quantile(rg_post, 0.025, na.rm=TRUE)
    rg_975  <- quantile(rg_post, 0.975, na.rm=TRUE)
    P_pos   <- mean(rg_post > 0, na.rm=TRUE)
    ess_rg  <- effectiveSize(rg_post)
    
    cat(sprintf("   rg coefficient = %.3f [95%% HPD: %.3f, %.3f] | Posterior Probability P(rg>0)=%.3f | ESS=%.0f\n",
                rg_mean, rg_025, rg_975, P_pos, ess_rg))
    
    # Biometric interpretation mapping
    interpretacao <- case_when(
      abs(rg_mean) < 0.2                    ~ "negligible",
      abs(rg_mean) < 0.4 & rg_mean > 0      ~ "weak positive",
      abs(rg_mean) < 0.4 & rg_mean < 0      ~ "weak negative",
      abs(rg_mean) < 0.7 & rg_mean > 0      ~ "moderate positive",
      abs(rg_mean) < 0.7 & rg_mean < 0      ~ "moderate negative",
      rg_mean > 0                           ~ "strong positive",
      TRUE                                  ~ "strong negative"
    )
    cat(sprintf("   Biometric Interpretation: %s\n", interpretacao))
    
    resultados[[nome]] <- list(
      t1=t1, t2=t2,
      rg=rg_mean, rg_025=rg_025, rg_975=rg_975,
      P_pos=P_pos, ESS=ess_rg,
      DIC=m$DIC,
      interpretacao=interpretacao
    )
  }
  cat("\n")
}

# =============================================================================
# 4. FINAL INTEGRATED SUMMARY MATRIX
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("GENETIC CORRELATIONS MATRIX — INTEGRATED SUMMARY\n")
cat(strrep("=", 60), "\n\n")

rg_tab <- do.call(rbind, lapply(names(resultados), function(nm) {
  r <- resultados[[nm]]
  data.frame(
    par           = nm,
    trait1        = r$t1,
    trait2        = r$t2,
    rg            = round(r$rg, 3),
    rg_025        = round(r$rg_025, 3),
    rg_975        = round(r$rg_975, 3),
    P_pos         = round(r$P_pos, 3),
    ESS           = round(r$ESS, 0),
    interpretacao = r$interpretacao,
    row.names     = NULL
  )
}))

print(rg_tab[, c("trait1","trait2","rg","rg_025","rg_975",
                 "P_pos","interpretacao")],
      row.names=FALSE)

# Export output data structures
write.csv(rg_tab,
          here("04_resultados/parametros_geneticos/bayesiano/rg_bivariado.csv"),
          row.names=FALSE)

t_total <- (proc.time() - t0_total)["elapsed"] / 60
cat(sprintf("\n✓ Bivariate routine module completed successfully — total runtime: %.1f minutes\n", t_total))
cat("  Exported Matrix:       rg_bivariado.csv\n")
cat("  Next pipeline module:  05_indice_resiliencia.R\n")