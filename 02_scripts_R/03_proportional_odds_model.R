# =============================================================================
# Proportional Odds Model (POM) for ordinal CABMV severity data
#
# OBJECTIVE:
#   Estimate P(score ≥ threshold | individual k, time t, family f)
#   controlling for time and family effects simultaneously
#   → Eliminates temporal confounding inherent to AUDPC
#   → Provides individual disease probability relative to the population
#   → Generalizable to any ordinal scoring scale
#
# OUTPUTS:
#   prob_severo_k — P(score ≥ 3) per individual × assessment date
#   pom_ind_sc    — Standardized individual-level disease pressure index
#   dataset_pom   — Formatted dataset optimized for downstream MCMCglmm analysis
# =============================================================================

library(here)
library(dplyr)
library(ordinal)    # clmm — cumulative link mixed model
library(ggplot2)
library(tidyr)

cat("\n", strrep("=", 60), "\n")
cat("POM — Proportional Odds Model\n")
cat(format(Sys.time(), "[%H:%M:%S]"), "\n")
cat(strrep("=", 60), "\n\n")

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================

severity_long <- readRDS(here("01_dados/processados/severity_long.rds")) # Load baseline long dataset (without day variable)
dataset       <- readRDS(here("01_dados/processados/dataset_completo.rds"))
dataset       <- as.data.frame(dataset)

# Compute absolute days PRIOR to dataset merging and mutation step
severity_long <- severity_long |>
  mutate(
    dia = as.numeric(data_aval - as.Date("2018-03-02"))
  )

# Incorporate family markers and associated cross metadata
familia_info <- dataset |>
  select(individuo, familia, parcela) |>
  distinct() |>
  mutate(individuo = as.integer(as.character(individuo)))

severity_long <- severity_long |>
  mutate(individuo = as.integer(individuo)) |>
  left_join(familia_info, by="individuo") |>
  filter(!is.na(planta)) |>
  mutate(
    nota_ord    = factor(planta, levels=1:4, ordered=TRUE),
    familia     = factor(familia),
    individuo_f = factor(individuo),
    dia_sc      = scale(dia)[,1]
  )

cat("Data preparation completed:\n")
cat("  Observations:", nrow(severity_long), "\n")
cat("  Individuals:", n_distinct(severity_long$individuo), "\n")
cat("  Score frequency distribution:\n")
print(table(severity_long$planta))

# =============================================================================
# 2. POM MODEL FITTING — CUMULATIVE LINK MIXED MODEL
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("Fitting CLMM (Cumulative Link Mixed Model)\n")
cat(strrep("=", 60), "\n\n")

# Full interaction model matrix specification:
# nota_ord ~ dia_sc * familia + (1|individuo)
# Controls for: time trend + full-sib family + random individual genetic effect
# The residual variance quantifies individual specific probability of severe disease

cat("Model 1 — time + family + random individual effect...\n")
t0 <- proc.time()

m_pom <- tryCatch(
  clmm(
    nota_ord ~ dia_sc * familia + (1|individuo_f),
    data    = severity_long,
    link    = "logit",
    Hess    = TRUE
  ),
  error = function(e) {
    cat("  Error detected:", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(m_pom)) {
  cat(sprintf("  Execution time: %.1f seconds\n", (proc.time()-t0)["elapsed"]))
  cat(sprintf("  Log-likelihood: %.2f\n", logLik(m_pom)))
  cat(sprintf("  AIC: %.2f\n\n", AIC(m_pom)))
  
  cat("Fixed effects coefficients:\n")
  print(round(coef(m_pom), 4))
  
  cat("\nRandom individual effect variance component:\n")
  print(m_pom$ST)
}

# Reduced main-effects model (omitting interaction term for nested evaluation)
cat("\nModel 2 — time + family (main effects model without interaction)...\n")
m_pom2 <- tryCatch(
  clmm(
    nota_ord ~ dia_sc + familia + (1|individuo_f),
    data    = severity_long,
    link    = "logit",
    Hess    = TRUE
  ),
  error = function(e) NULL
)

if (!is.null(m_pom2)) {
  cat(sprintf("  AIC: %.2f\n", AIC(m_pom2)))
  cat(sprintf("  ΔAIC vs Model 1: %.2f\n\n",
              AIC(m_pom2) - AIC(m_pom)))
}

# Select the optimal model based on Information Criteria parsimony
m_final <- if (!is.null(m_pom) &&
               (is.null(m_pom2) || AIC(m_pom) <= AIC(m_pom2)))
  m_pom else m_pom2

# =============================================================================
# 3. EXTRACT INDIVIDUAL P(score ≥ 3) — FIXED ESTIMATION
# Maps back to 'sev_pred' (family subset) to bypass missing records
# Manually derives probability matrices from linear predictors
# =============================================================================

cat(strrep("=", 60), "\n")
cat("Extracting marginal individual probabilities\n")
cat(strrep("=", 60), "\n\n")

if (!is.null(m_final)) {
  
  # ── Extract model terms ──
  coefs      <- coef(m_final)
  thresholds <- coefs[grep("\\|", names(coefs))]
  betas      <- coefs[!grepl("\\|", names(coefs))]
  
  # ── Filter records with valid family structures ──
  sev_pred <- severity_long |>
    filter(!is.na(familia)) |>
    mutate(familia = droplevels(familia))
  
  cat(sprintf("Observations with family data: %d\n", nrow(sev_pred)))
  cat(sprintf("Observations missing family data: %d (excluded from summary)\n",
              nrow(severity_long) - nrow(sev_pred)))
  
  # ── Linear predictor calculation (Fixed terms) ──
  X   <- model.matrix(~ dia_sc * familia, data=sev_pred)[, -1]
  eta <- as.numeric(X %*% betas[names(betas) %in% colnames(X)])
  
  # ── Cumulative link probabilities mapping ──
  # P(score ≤ k) = plogis(threshold_k - eta)
  P_cum1 <- plogis(thresholds[1] - eta)
  P_cum2 <- plogis(thresholds[2] - eta)
  P_cum3 <- plogis(thresholds[3] - eta)
  
  # ── Conditional categorical probabilities calculation ──
  sev_pred$P_nota1 <- P_cum1
  sev_pred$P_nota2 <- P_cum2 - P_cum1
  sev_pred$P_nota3 <- P_cum3 - P_cum2
  sev_pred$P_nota4 <- 1 - P_cum3
  
  # ── Marginal probability of severe disease: P(score ≥ 3) = P(score=3) + P(score=4) ──
  sev_pred$P_severo <- sev_pred$P_nota3 + sev_pred$P_nota4
  
  cat(sprintf("\nProbability distribution profile for P(score ≥ 3):\n"))
  cat(sprintf("  Mean=%.3f | SD=%.3f | Range: [%.3f, %.3f]\n",
              mean(sev_pred$P_severo),
              sd(sev_pred$P_severo),
              min(sev_pred$P_severo),
              max(sev_pred$P_severo)))
  
  # ── Collapse to individual level summary statistics ──
  pom_ind <- sev_pred |>
    group_by(individuo) |>
    summarise(
      P_severo_media = mean(P_severo, na.rm=TRUE),
      P_severo_pico  = max(P_severo, na.rm=TRUE),
      P_severo_final = mean(tail(P_severo, 5), na.rm=TRUE),
      P_severo_var   = var(P_severo, na.rm=TRUE),
      n_aval         = n(),
      familia        = first(familia),
      .groups        = "drop"
    )
  
  cat(sprintf("\nIndividual P(score ≥ 3) matrix summary (%d genotypes evaluated):\n", nrow(pom_ind)))
  cat(sprintf("  Mean P_severo: population mean=%.3f | SD=%.3f\n",
              mean(pom_ind$P_severo_media),
              sd(pom_ind$P_severo_media)))
  cat(sprintf("  Peak P_severo:  population mean=%.3f | SD=%.3f\n",
              mean(pom_ind$P_severo_pico),
              sd(pom_ind$P_severo_pico)))
  
  # ── Append right-censored Fusarium wilt mortality records and standardize metrics ──
  mortos_fusarium <- c(10, 16, 20, 33, 42, 48, 61, 78, 87, 91,
                       123, 129, 142, 146, 174, 188, 196, 203,
                       205, 208, 219)
  
  pom_ind <- pom_ind |>
    mutate(
      morreu_fusarium   = as.integer(individuo %in% mortos_fusarium),
      n_aval_sc         = scale(n_aval)[,1],
      P_severo_media_sc = scale(P_severo_media)[,1],
      P_severo_pico_sc  = scale(P_severo_pico)[,1],
      P_severo_var_sc   = scale(P_severo_var)[,1],
      POM_idx           = scale(P_severo_media)[,1] +
        scale(P_severo_pico)[,1] +
        scale(P_severo_var)[,1]
    )
  
  cat(sprintf("\nFusarium wilt mortality rate in the target subset: %d out of %d genotypes\n",
              sum(pom_ind$morreu_fusarium), nrow(pom_ind)))
  
  # ── Identify contrasting phenotypic tails (Resilient vs Susceptible) ──
  cat("\nTop 10 most RESILIENT genotypes (Lowest posterior P_severo_media):\n")
  print(head(pom_ind[order(pom_ind$P_severo_media),
                     c("individuo","familia","P_severo_media",
                       "P_severo_pico","n_aval","morreu_fusarium")], 10),
        row.names=FALSE)
  
  cat("\nTop 10 most SUSCEPTIBLE genotypes (Highest posterior P_severo_media):\n")
  print(head(pom_ind[order(pom_ind$P_severo_media, decreasing=TRUE),
                     c("individuo","familia","P_severo_media",
                       "P_severo_pico","n_aval","morreu_fusarium")], 10),
        row.names=FALSE)
  
} else {
  cat("ERROR: m_final object is NULL — fit the CLMM routine prior to executing this block.\n")
}

# =============================================================================
# SECTION 4 — CROSS-METRIC CORRELATIONS: POM VS AUDPC VS BDSIP AND AGRONOMIC TRAITS
# =============================================================================

dataset_novo <- readRDS(here("01_dados/processados/dataset_aacpd_janela.rds"))
bdsip_std    <- readRDS(here("01_dados/processados/bdsip_completo.rds"))

aacpd_ind <- dataset_novo |>
  group_by(individuo) |>
  summarise(aacpd_media = mean(aacpdm, na.rm=TRUE), .groups="drop") |>
  mutate(individuo = as.integer(as.character(individuo)))

comp_df <- pom_ind |>
  left_join(aacpd_ind, by="individuo") |>
  left_join(bdsip_std |> select(individuo, BDSIP), by="individuo") |>
  filter(!is.na(aacpd_media))

cat("Cross-metric Pearson correlation coefficients:\n")
cat(sprintf("  Mean POM × Mean AUDPC: r=%.3f\n",
            cor(comp_df$P_severo_media, comp_df$aacpd_media,
                use="pairwise.complete.obs")))
cat(sprintf("  Mean POM × BDSIP score: r=%.3f\n",
            cor(comp_df$P_severo_media, comp_df$BDSIP,
                use="pairwise.complete.obs")))
cat(sprintf("  BDSIP score × Mean AUDPC: r=%.3f\n",
            cor(comp_df$BDSIP, comp_df$aacpd_media,
                use="pairwise.complete.obs")))

# Map disease metrics against primary agronomic fruit traits
dataset_pom_corr <- dataset_novo |>
  mutate(individuo = as.integer(as.character(individuo))) |>
  left_join(
    pom_ind |> select(individuo, P_severo_media,
                      P_severo_pico, n_aval_sc, morreu_fusarium),
    by="individuo"
  ) |>
  left_join(
    bdsip_std |> select(individuo, BDSIP) |>
      mutate(individuo = as.integer(individuo),
             BDSIP = as.numeric(BDSIP)),
    by="individuo"
  )

# Basic metadata assertions
cat("BDSIP vector class mapping:", class(dataset_pom_corr$BDSIP), "\n")
cat("Missing data tally (NAs) in BDSIP field:", sum(is.na(dataset_pom_corr$BDSIP)), "\n")

# Generate comparative correlation matrix
cat("\nCorrelation Profile: Disease Metrics vs Fruit Agronomic Traits:\n")
cat(sprintf("  %-20s %8s %8s %8s\n", "Trait", "POM", "AUDPC", "BDSIP"))

for (tr in c("comp_mm","diam_mm","massa_fruto_g",
             "rend_polpa_pct","esp_casca_mm","sst_brix")) {
  r_pom   <- cor(dataset_pom_corr$P_severo_media, dataset_pom_corr[[tr]],
                 use="pairwise.complete.obs")
  r_aacpd <- cor(dataset_pom_corr$aacpdm, dataset_pom_corr[[tr]],
                 use="pairwise.complete.obs")
  r_bdsip <- cor(dataset_pom_corr$BDSIP, dataset_pom_corr[[tr]],
                 use="pairwise.complete.obs")
  cat(sprintf("  %-20s %+8.3f %+8.3f %+8.3f\n",
              tr, r_pom, r_aacpd, r_bdsip))
}

# =============================================================================
# SECTION 5 — FORMAT FOR MCMCGLMM RUNS AND EXPORT DATA
# =============================================================================

dataset_pom <- dataset_novo |>
  mutate(individuo = as.integer(as.character(individuo))) |>
  left_join(
    pom_ind |> select(individuo, P_severo_media, P_severo_media_sc,
                      P_severo_pico_sc, P_severo_var_sc,
                      n_aval_sc, morreu_fusarium),
    by="individuo"
  )

# Enforce zero-imputation for scaled vectors (assigning population background mean)
for (v in c("P_severo_media_sc","P_severo_pico_sc",
            "P_severo_var_sc","n_aval_sc","morreu_fusarium")) {
  dataset_pom[[v]][is.na(dataset_pom[[v]])] <- 0
}

saveRDS(m_final,     here("04_resultados/modelos/m_pom_severidade.rds"))
saveRDS(pom_ind,     here("01_dados/processados/pom_indices.rds"))
saveRDS(dataset_pom, here("01_dados/processados/dataset_pom.rds"))
write.csv(pom_ind,   here("01_dados/processados/pom_indices.csv"),
          row.names=FALSE)

cat("✓ Data objects saved successfully.\n")
cat(sprintf("  dataset_pom table dimensions: %d rows\n", nrow(dataset_pom)))
cat(sprintf("  pom_ind population tracking count: %d genotypes\n", nrow(pom_ind)))

# =============================================================================
# 6. FILE EXPORT REDUNDANCY GUARD
# =============================================================================

saveRDS(m_final,   here("04_resultados/modelos/m_pom_severidade.rds"))
saveRDS(pom_ind,   here("01_dados/processados/pom_indices.rds"))
saveRDS(dataset_pom, here("01_dados/processados/dataset_pom.rds"))

write.csv(pom_ind,
          here("01_dados/processados/pom_indices.csv"),
          row.names=FALSE)

# Export long-format data frame mapped with individual posterior vectors
saveRDS(severity_long,
        here("01_dados/processados/severity_long_pom.rds"))

cat("✓ POM script pipeline completed.\n")
cat("  m_pom_severidade.rds — Fitted CLMM model checkpoint\n")
cat("  pom_indices.rds      — Posterior tracking of individual P(score≥3)\n")
cat("  dataset_pom.rds      — Compiled analytical table for downstream MCMCglmm execution\n")
cat("\nNext step: Execute model screening using DIC selection (M_AUDPC vs M_POM vs M_POM+).\n")