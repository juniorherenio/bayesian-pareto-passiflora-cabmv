# =============================================================================
# DIC Screening — M_AUDPC vs M_POM vs M_POM+
# Agronomic Traits: esp_casca_mm (rind thickness) and rend_polpa_pct (pulp yield)
# =============================================================================

library(here)
library(MCMCglmm)
library(Matrix)

dataset_pom <- readRDS(here("01_dados/processados/dataset_pom.rds"))
Amat        <- readRDS(here("01_dados/processados/Amat.rds"))

# -----------------------------------------------------------------------------
# Data Setup and Factor Level Alignment
# -----------------------------------------------------------------------------
dat <- as.data.frame(dataset_pom)
dat$bloco     <- as.factor(dat$bloco)
dat$parcela   <- as.factor(dat$parcela)
dat$individuo <- factor(as.factor(dat$individuo), levels=rownames(Amat))

dat$estacao_ano <- as.character(dat$estacao_ano)
dat$estacao_ano[dat$estacao_ano == "verao_2018"] <- "verao_2019"
dat$estacao_ano <- factor(dat$estacao_ano,
                          levels=c("inverno_2018","primavera_2018",
                                   "verao_2019","outono_2019","inverno_2019"))
dat <- dat[!is.na(dat$estacao_ano), ]
dat <- droplevels(dat)

dat$estacao <- dat$estacao_ano
levels(dat$estacao) <- c("i18","p18","v19","o19","i19")

dat$animal <- dat$individuo
dat$pe     <- dat$individuo

# -----------------------------------------------------------------------------
# Disease Covariates Transformation
# -----------------------------------------------------------------------------
dat$aacpdm_sc <- scale(dat$aacpdm)[,1]
dat$aacpdm_sc[is.na(dat$aacpdm_sc)] <- 0

# -----------------------------------------------------------------------------
# Kinship Matrix (A-matrix) Subsetting and Inversion
# -----------------------------------------------------------------------------
ind_obs     <- levels(droplevels(dat$individuo))
Amat_sub    <- Amat[ind_obs, ind_obs]
Ainv_sparse <- as(solve(Amat_sub), "dgCMatrix") # Sparse inverse relationship matrix
dat$individuo <- factor(dat$individuo, levels=rownames(Amat_sub))
dat$animal    <- dat$individuo
dat$pe        <- dat$individuo

n_est    <- nlevels(dat$estacao)
# MCMC chain settings for screening
NITT_tri <- 31000; BURNIN_tri <- 1000; THIN_tri <- 30

# -----------------------------------------------------------------------------
# Model Screening Execution Function
# -----------------------------------------------------------------------------
rodar <- function(tr, cov_nome, cov_f) {
  Vp    <- var(dat[[tr]], na.rm=TRUE)
  # Inverse Wishart prior structures specification
  prior <- list(
    G = list(
      G1 = list(V=diag(n_est)*Vp*0.2, nu=n_est+1), # Additive genetic variance per season (idh)
      G2 = list(V=Vp*0.2, nu=0.002),               # Permanent environmental effect (pe)
      G3 = list(V=Vp*0.05, nu=0.002)               # Plot effect
    ),
    R = list(V=diag(n_est)*Vp*0.4, nu=n_est+1)     # Residual variance structure per season (idh)
  )
  fixed_f <- as.formula(paste0(tr, " ~ bloco + estacao + ", cov_f))
  set.seed(2026)
  m <- tryCatch(
    MCMCglmm(
      fixed=fixed_f,
      random= ~ idh(estacao):animal + pe + parcela,
      rcov  = ~ idh(estacao):units,
      ginverse=list(animal=Ainv_sparse),
      prior=prior, data=dat,
      nitt=NITT_tri, burnin=BURNIN_tri, thin=THIN_tri,
      verbose=FALSE
    ),
    error=function(e) { cat("Execution Error:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(m)) {
    # Extract posterior beta coefficient for the specific disease predictor
    cov_col <- grep(strsplit(cov_f, " \\+")[[1]][1],
                    colnames(m$Sol), value=TRUE)[1]
    beta <- if (!is.na(cov_col)) mean(m$Sol[, cov_col]) else NA
    p_neg <- if (!is.na(cov_col)) mean(m$Sol[, cov_col] < 0) else NA
    cat(sprintf("  %-25s DIC=%7.1f | β=%+.3f | P(β<0)=%.3f\n",
                cov_nome, m$DIC, beta, p_neg))
    return(list(DIC=m$DIC, beta=beta, p_neg=p_neg))
  }
  return(list(DIC=NA, beta=NA, p_neg=NA))
}

# Target epidemiological candidate models for selection
modelos <- list(
  list(nome="AUDPC", formula="aacpdm_sc"),
  list(nome="POM",   formula="P_severo_media_sc"),
  list(nome="POM+",  formula="P_severo_media_sc + n_aval_sc + morreu_fusarium")
)

resultados <- data.frame()

# Iterate screening loop across agronomic target traits
for (tr in c("esp_casca_mm", "rend_polpa_pct")) {
  cat(sprintf("\nTarget Trait: %s\n", tr))
  cat(strrep("-", 65), "\n")
  for (m in modelos) {
    r <- rodar(tr, m$nome, m$formula)
    resultados <- rbind(resultados, data.frame(
      trait=tr, modelo=m$nome,
      DIC=round(r$DIC,1),
      beta=round(r$beta,3),
      P_beta_neg=round(r$p_neg,3)
    ))
  }
}

cat("\n\nFinal DIC Selection Matrix Table:\n")
print(resultados, row.names=FALSE)