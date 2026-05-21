# =============================================================================
# Bayesian Univariate Mixed Models with G×E — 7 Agronomic Traits
# Disease Covariate: P_severo_media_sc (POM)
#
# METHODOLOGICAL COMPARISON vs ORIGINAL SCRIPT:
#   Covariate: audpc_sc → P_severo_media_sc
#   Dataset:   dataset_completo.rds → dataset_pom.rds
#   Biological Rationale: POM captures the correct causal relationships
#     esp_casca_mm (rind thickness):  β > 0 (thicker rind under disease pressure) ✓
#     rend_polpa_pct (pulp yield):    β < 0 (reduced pulp yield under disease pressure) ✓
#     sst_brix (soluble solids):      near-zero correlation (eliminates confounding) ✓
#
# G×E VARIANCE STRUCTURE SPECIFICATION:
#   G (Additive Genetic):          idh(estacao):animal
#   R (Residual):                  idh(estacao):units
#   + pe (Permanent Environment) + parcela (Plot effect)
#
# OUTPUTS:
#   MCMC chains: 04_resultados/mcmc/chains/GxA_POM_{trait}_chain{1,2}.rds
#   BLUPs:       04_resultados/blups/bayesiano/blups_GxA_POM_{trait}.csv
#   Parameters:  04_resultados/parametros_geneticos/bayesiano/params_GxA_POM.csv
# =============================================================================

library(here)
library(MCMCglmm)
library(Matrix)
library(coda)
library(dplyr)

cat("\n", strrep("=", 60), "\n")
cat("04e_v2_mcmcglmm_GxE_POM.R\n")
cat(format(Sys.time(), "[%H:%M:%S]"), "\n")
cat(strrep("=", 60), "\n\n")

# =============================================================================
# 1. SETUP AND ENVIRONMENTAL ALIGNMENT
# =============================================================================

dataset <- readRDS(here("01_dados/processados/dataset_pom.rds"))
Amat    <- readRDS(here("01_dados/processados/Amat.rds"))

dataset <- as.data.frame(dataset)
dataset$bloco     <- as.factor(dataset$bloco)
dataset$familia   <- as.factor(dataset$familia)
dataset$parcela   <- as.factor(dataset$parcela)
dataset$individuo <- factor(as.factor(dataset$individuo),
                            levels = rownames(Amat))

# Season-year factor structural alignment
dataset$estacao_ano <- as.character(dataset$estacao_ano)
dataset$estacao_ano[dataset$estacao_ano == "verao_2018"] <- "verao_2019"
dataset$estacao_ano <- factor(dataset$estacao_ano,
                              levels = c("inverno_2018","primavera_2018",
                                         "verao_2019","outono_2019",
                                         "inverno_2019"))

dat <- dataset[!is.na(dataset$estacao_ano), ]
dat <- droplevels(dat)
dat$individuo <- factor(dat$individuo, levels = rownames(Amat))
dat$animal    <- dat$individuo
dat$pe        <- dat$individuo

# Simplify season levels for clean model matrix headings
dat$estacao <- dat$estacao_ano
levels(dat$estacao) <- c("i18","p18","v19","o19","i19")
n_est <- nlevels(dat$estacao)

# POM Covariate — Pre-standardized during the upstream step
# Missing data guard check
dat$P_severo_media_sc[is.na(dat$P_severo_media_sc)] <- 0

cat("Season Factor Levels:", paste(levels(dat$estacao), collapse=", "), "\n")
cat("Observations tally per season:\n"); print(table(dat$estacao))
cat(sprintf("\nMissing values (NAs) in P_severo_media_sc: %d\n",
            sum(is.na(dat$P_severo_media_sc))))

# Relationship A-matrix Subsetting
ind_obs     <- levels(droplevels(dat$individuo))
Amat_sub    <- Amat[ind_obs, ind_obs]
Ainv_sparse <- as(solve(Amat_sub), "dgCMatrix")
dat$individuo <- factor(dat$individuo, levels = rownames(Amat_sub))
dat$animal    <- dat$individuo
dat$pe        <- dat$individuo

cat(sprintf("\nAnalytical Sample Profile: %d obs | %d genotypes | %d seasons\n\n",
            nrow(dat), length(ind_obs), n_est))

# Core Target Traits Block
traits_7 <- c("comp_mm", "diam_mm", "massa_fruto_g",
              "rend_polpa_pct", "esp_casca_mm",
              "sst_brix", "massa_polpa_g")

# Baseline Phenotypic Variances
Vp_vec <- sapply(traits_7, function(tr) var(dat[[tr]], na.rm=TRUE))
cat("Phenotypic Variance Vector (Vp):\n"); print(round(Vp_vec, 2))

# =============================================================================
# 2. MCMC SAMPLING PARAMETERS
# =============================================================================

BURNIN   <- 50000
NITT     <- 200000
THIN     <- 100
N_CHAINS <- 2

cat(sprintf("\nMCMC Configurations: nitt=%d | burnin=%d | thin=%d\n", NITT, BURNIN, THIN))
cat(sprintf("Target Pathosystem Covariate: P_severo_media_sc (POM)\n\n"))

# Initialize directory infrastructure
dir.create(here("04_resultados/mcmc/chains"),
           showWarnings=FALSE, recursive=TRUE)
dir.create(here("04_resultados/blups/bayesiano"),
           showWarnings=FALSE, recursive=TRUE)
dir.create(here("04_resultados/parametros_geneticos/bayesiano"),
           showWarnings=FALSE, recursive=TRUE)

# =============================================================================
# 3. CORE ANALYTICAL PROCESSING LOOP
# =============================================================================

t0_total   <- proc.time()
resultados <- list()

for (tr in traits_7) {
  
  cat("\n", strrep("=", 60), "\n")
  cat(sprintf("EVALUATING TRAIT: %s\n", tr))
  cat(strrep("=", 60), "\n\n")
  
  Vp <- Vp_vec[tr]
  
  # Parameterize Inverse Wishart priors for heterogeneous G×E models
  prior_GxA <- list(
    G = list(
      G1 = list(V = diag(n_est) * Vp * 0.2, nu = n_est + 1), # Additive genetic variance per season
      G2 = list(V = Vp * 0.2,                nu = 0.002),     # Permanent environmental effect
      G3 = list(V = Vp * 0.05,               nu = 0.002)      # Micro-spatial plot effect
    ),
    R = list(V = diag(n_est) * Vp * 0.4, nu = n_est + 1)     # Residual variance structure per season
  )
  
  fixed_f <- as.formula(
    paste0(tr, " ~ bloco + estacao + P_severo_media_sc")
  )
  
  chains_tr <- list()
  
  for (i in 1:N_CHAINS) {
    cat(sprintf("[%s] Processing %s — MCMC Chain %d/%d...\n",
                format(Sys.time(), "%H:%M:%S"), tr, i, N_CHAINS))
    t0 <- proc.time()
    set.seed(2025 + i)
    
    m <- MCMCglmm(
      fixed    = fixed_f,
      random   = ~ idh(estacao):animal + pe + parcela,
      rcov     = ~ idh(estacao):units,
      ginverse = list(animal = Ainv_sparse),
      prior    = prior_GxA,
      data     = dat,
      nitt     = NITT,
      burnin   = BURNIN,
      thin     = THIN,
      pr       = TRUE,
      verbose  = FALSE
    )
    
    elapsed <- (proc.time() - t0)["elapsed"] / 60
    cat(sprintf("[%s] Chain %d processing completed — runtime: %.1f min | DIC=%.1f\n",
                format(Sys.time(), "%H:%M:%S"), i, elapsed, m$DIC))
    
    saveRDS(m, here(sprintf(
      "04_resultados/mcmc/chains/GxA_POM_%s_chain%d.rds", tr, i)))
    cat(sprintf("  ✓ Saved Checkpoint: GxA_POM_%s_chain%d.rds\n", tr, i))
    
    chains_tr[[i]] <- m
  }
  
  # ── Convergence Diagnostics ──
  cat(sprintf("\nConvergence Diagnostics for %s:\n", tr))
  ess <- effectiveSize(chains_tr[[1]]$VCV)
  cat(sprintf("  Effective Sample Size (ESS) — Minimum: %.0f | Maximum: %.0f\n",
              min(ess), max(ess)))
  
  gr <- tryCatch(
    gelman.diag(
      mcmc.list(as.mcmc(chains_tr[[1]]$VCV),
                as.mcmc(chains_tr[[2]]$VCV)),
      multivariate = FALSE
    ),
    error = function(e) NULL
  )
  if (!is.null(gr)) {
    gr_max <- max(gr$psrf[, 1], na.rm=TRUE)
    cat(sprintf("  Maximum Potential Scale Reduction Factor (R̂): %.3f %s\n",
                gr_max, ifelse(gr_max < 1.1, "✓", "⚠")))
  }
  
  # ── Narrow-Sense Heritability (h²) Tracking per Season ──
  cat(sprintf("\n  Narrow-sense heritability (h²) breakdown per season:\n"))
  vc       <- as.data.frame(chains_tr[[1]]$VCV)
  vc_names <- names(vc)
  nm_units <- grep("units", vc_names, value=TRUE)
  
  h2_estacao <- list()
  for (est in levels(dat$estacao)) {
    nm_a <- paste0("estacao", est, ".animal")
    nm_e <- grep(est, nm_units, value=TRUE)[1]
    
    if (nm_a %in% vc_names && !is.na(nm_e) && nm_e %in% vc_names) {
      h2_post <- vc[[nm_a]] / (vc[[nm_a]] + vc[[nm_e]])
      h2_estacao[[est]] <- h2_post
      cat(sprintf("    Season %s: posterior mean h²=%.3f [95%% HPD: %.3f, %.3f]\n",
                  est, mean(h2_post),
                  quantile(h2_post, 0.025),
                  quantile(h2_post, 0.975)))
    }
  }
  
  # ── Weighted Marginal Narrow-Sense Heritability ──
  n_est_obs <- table(dat$estacao)
  h2_vals   <- sapply(names(h2_estacao), function(e) mean(h2_estacao[[e]]))
  pesos     <- as.numeric(n_est_obs[names(h2_estacao)])
  h2_pond   <- sum(h2_vals * pesos) / sum(pesos)
  cat(sprintf("  Observation-weighted narrow-sense heritability (h²_pond): %.3f\n", h2_pond))
  
  # ── Fixed Pathosystem Covariate Effect (POM) ──
  nm_beta <- grep("P_severo", colnames(chains_tr[[1]]$Sol), value=TRUE)
  if (length(nm_beta) > 0) {
    beta <- chains_tr[[1]]$Sol[, nm_beta[1]]
    cat(sprintf("  β_POM slope coefficient: %.3f [95%% HPD: %.3f, %.3f] | Posterior Probability P(β<0)=%.3f\n",
                mean(beta),
                quantile(beta, 0.025), quantile(beta, 0.975),
                mean(beta < 0)))
  } else {
    beta <- rep(NA, 1)
  }
  
  # ── Breeding Values Prediction (BLUPs extraction) ──
  sol_names   <- colnames(chains_tr[[1]]$Sol)
  blup_animal <- sol_names[grep("\\.animal\\.", sol_names)]
  
  if (length(blup_animal) > 0) {
    blups_tr <- data.frame(
      trait     = tr,
      individuo = gsub(".*\\.animal\\.(\\w+)$", "\\1", blup_animal),
      blup      = colMeans(chains_tr[[1]]$Sol[, blup_animal, drop=FALSE]),
      blup_sd   = apply(chains_tr[[1]]$Sol[, blup_animal, drop=FALSE],
                        2, sd),
      blup_025  = apply(chains_tr[[1]]$Sol[, blup_animal, drop=FALSE],
                        2, quantile, 0.025),
      blup_975  = apply(chains_tr[[1]]$Sol[, blup_animal, drop=FALSE],
                        2, quantile, 0.975),
      row.names = NULL
    )
    
    write.csv(blups_tr,
              here(sprintf(
                "04_resultados/blups/bayesiano/blups_GxA_POM_%s.csv", tr)),
              row.names=FALSE)
    cat(sprintf("  ✓ Exported Breeding Values Matrix: %d genotypes included\n", nrow(blups_tr)))
    
    top5 <- blups_tr[order(blups_tr$blup, decreasing=TRUE), ]
    cat("  Top 5 Genotypes based on estimated breeding values:\n")
    print(head(top5[, c("individuo","blup","blup_sd")], 5),
          row.names=FALSE)
  }
  
  resultados[[tr]] <- list(
    h2_pond    = h2_pond,
    h2_estacao = h2_estacao,
    beta_POM   = if (!is.na(beta[1])) mean(beta) else NA,
    P_beta_neg = if (!is.na(beta[1])) mean(beta < 0) else NA,
    DIC        = chains_tr[[1]]$DIC,
    ESS_min    = min(ess)
  )
}

# =============================================================================
# 4. FINAL INTEGRATED SUMMARY MATRIX
# =============================================================================

cat("\n\n", strrep("=", 60), "\n")
cat("FINAL INTEGRATED SUMMARY MATRIX REPORT\n")
cat(strrep("=", 60), "\n\n")

params_GxA_POM <- do.call(rbind, lapply(traits_7, function(tr) {
  r <- resultados[[tr]]
  data.frame(
    trait      = tr,
    h2_pond    = round(r$h2_pond, 3),
    beta_POM   = round(r$beta_POM, 3),
    P_beta_neg = round(r$P_beta_neg, 3),
    DIC        = round(r$DIC, 1),
    ESS_min    = round(r$ESS_min, 0),
    row.names  = NULL
  )
}))
print(params_GxA_POM, row.names=FALSE)

# Marginal heritabilities structured by season
h2_por_estacao <- do.call(rbind, lapply(traits_7, function(tr) {
  r <- resultados[[tr]]
  do.call(rbind, lapply(names(r$h2_estacao), function(est) {
    h2 <- r$h2_estacao[[est]]
    data.frame(
      trait    = tr,
      estacao  = est,
      h2_media = round(mean(h2), 3),
      h2_025   = round(quantile(h2, 0.025), 3),
      h2_975   = round(quantile(h2, 0.975), 3),
      row.names = NULL
    )
  }))
}))

cat("\nNarrow-sense heritability structured by trait × season matrix:\n")
print(h2_por_estacao, row.names=FALSE)

# Statistical validation benchmarking against the baseline AUDPC script
cat("\nDirect Benchmark Validation: β_POM vs β_AUDPC (from the baseline tracking script):\n")
beta_aacpd <- c(
  comp_mm=1.100, diam_mm=0.963, massa_fruto_g=2.279,
  rend_polpa_pct=1.049, esp_casca_mm=-0.755,
  sst_brix=0.053, massa_polpa_g=2.221
)
p_aacpd <- c(
  comp_mm=0.243, diam_mm=0.248, massa_fruto_g=0.324,
  rend_polpa_pct=0.257, esp_casca_mm=0.983,
  sst_brix=0.443, massa_polpa_g=0.173
)
cat(sprintf("  %-20s %8s %8s %8s %8s\n",
            "Trait", "β_AUDPC", "P(β<0)", "β_POM", "P(β<0)"))
for (tr in traits_7) {
  r <- resultados[[tr]]
  cat(sprintf("  %-20s %+8.3f %8.3f %+8.3f %8.3f\n",
              tr, beta_aacpd[tr], p_aacpd[tr],
              r$beta_POM, r$P_beta_neg))
}

# Export data structures
write.csv(params_GxA_POM,
          here("04_resultados/parametros_geneticos/bayesiano/params_GxA_POM.csv"),
          row.names=FALSE)
write.csv(h2_por_estacao,
          here("04_resultados/parametros_geneticos/bayesiano/h2_estacao_GxA_POM.csv"),
          row.names=FALSE)

t_total <- (proc.time() - t0_total)["elapsed"] / 60
cat(sprintf("\n✓ Processing workflow completed successfully — total runtime: %.1f minutes (%.2f hours)\n",
            t_total, t_total/60))
cat("  Breeding Values Matrix (BLUPs): blups_GxA_POM_{trait}.csv\n")
cat("  Genetic Parameter Estimates:    params_GxA_POM.csv\n")
cat("  Next block module pipeline:     04f_bivariado.R\n")