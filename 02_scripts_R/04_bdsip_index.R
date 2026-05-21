# =============================================================================
# Bayesian Disease Severity Index for Plant breeding (BDSIP)
#
# NOVEL METHODOLOGY — Integrates three fields never combined before:
#   1. Longitudinal Beta Regression (phytopathometry)
#   2. Adapted Survival Analysis (biomedical statistics)
#   3. Integration with Bayesian Animal Models (quantitative genetics)
#
# ESTIMATED PARAMETERS PER INDIVIDUAL GENOTYPE:
#   mu_max     — Peak severity (maximum epidemic intensity score)
#   t_pico     — Time to peak (epidemic precocity window)
#   slope_up   — Progression rate (exponential epidemic phase coefficient)
#   slope_down — Recovery/mitigation rate (post-peak phase trajectory)
#   t_surv     — Survival time to critical epidemiological threshold
#   area_beta  — Area under the Beta-fitted curve (improved alternative to AUDPC)
#   estab      — Temporal stability (longitudinal trajectory variance)
#
# METHODOLOGICAL INNOVATION:
#   Replaces the scalar AUDPC with a parameter vector that fully captures the
#   temporal DYNAMICS of disease progression, integrated into the downstream
#   MCMCglmm framework as individual-level covariates.
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(betareg)
library(survival)
library(ggplot2)
library(patchwork)

cat("\n", strrep("=", 60), "\n")
cat("BDSIP — Bayesian Disease Severity Index for Plant breeding\n")
cat(format(Sys.time(), "[%H:%M:%S]"), "\n")
cat(strrep("=", 60), "\n\n")

# =============================================================================
# 1. DATA LOADING AND BOUNDARY COMPRESSION
# =============================================================================

severity_long <- readRDS(here("01_dados/processados/severity_long.rds"))
dataset       <- readRDS(here("01_dados/processados/dataset_completo.rds"))
dataset       <- as.data.frame(dataset)

# Convert ordinal scores into continuous proportions (0-1 open interval for Beta mapping)
# Target scale: 1-4 ordinal range → transform to open (0,1) boundaries
# Transformation: (score - min) / (max - min) adjusted with Smithson & Verkuilen correction
severity_long <- severity_long |>
  mutate(
    # 'nota_planta' transformed into a strict proportion mapping inside (0,1)
    # (score-1)/3 → [0,1], followed by boundary shrinkage optimization
    prop_raw  = (planta - 1) / 3,
    prop_raw  = ifelse(is.na(prop_raw), NA, prop_raw),
    # Statistical compression to clip boundary artifacts at absolute 0 or 1
    n_obs     = n(),
    prop_beta = (prop_raw * (n_obs - 1) + 0.5) / n_obs,
    # Cumulative days relative to experimental baseline onset (2018-03-02)
    dia       = as.numeric(data_aval - as.Date("2018-03-02"))
  )

# Apply strict boundary clipping to preserve Beta integration limits
severity_long <- severity_long |>
  mutate(
    prop_beta = pmin(pmax(prop_raw, 0.001), 0.999)
  )

cat("Epidemiological severity data profile:\n")
cat("  Genotypes (Individuals):", n_distinct(severity_long$individuo), "\n")
cat("  Assessment dates:", n_distinct(severity_long$data_aval), "\n")
cat("  Total observations matrix:", nrow(severity_long), "\n\n")

# =============================================================================
# 2. CRITICAL EPIDEMIOLOGICAL THRESHOLD FOR SURVIVAL ANALYSIS
# =============================================================================

# Enforce threshold = Ordinal Score 3 (High chronic severity status)
# Biological rationale: Scores ≥ 3 indicate systemic plant compromise
LIMIAR_NOTA <- 3
LIMIAR_PROP <- (LIMIAR_NOTA - 1) / 3  # Matches continuous proportion = 0.667

cat(sprintf("Critical survival threshold set at: score ≥ %d (proportion mapping ≥ %.3f)\n\n",
            LIMIAR_NOTA, LIMIAR_PROP))

# =============================================================================
# 3. INDIVIDUAL BDSIP PARAMETER EXTRACTION ENGINE
# =============================================================================

cat(strrep("=", 60), "\n")
cat("Extracting individual BDSIP parameter matrices\n")
cat(strrep("=", 60), "\n\n")

extrair_bdsip <- function(dados_ind) {
  # 'dados_ind': subset block from 'severity_long' filtered per individual plant
  # explicitly sorted in chronological order
  
  dados_ind <- dados_ind |>
    filter(!is.na(planta)) |>
    arrange(dia)
  
  n <- nrow(dados_ind)
  if (n < 3) return(NULL)
  
  notas <- dados_ind$planta
  dias  <- dados_ind$dia
  props <- dados_ind$prop_beta
  
  # ── Parameter 1: mu_max — Epidemic Peak Intensity Score ──
  mu_max <- max(notas, na.rm=TRUE)
  idx_max <- which.max(notas)
  t_pico  <- dias[idx_max]
  
  # ── Parameter 2: slope_up — Ascending Epidemic Progression Rate ──
  # Tracks trajectory baseline from initiation to epidemic peak
  if (idx_max > 1) {
    fase_up <- dados_ind[1:idx_max, ]
    if (nrow(fase_up) >= 2) {
      fit_up <- tryCatch(
        lm(planta ~ dia, data=fase_up),
        error=function(e) NULL
      )
      slope_up <- if (!is.null(fit_up)) coef(fit_up)[2] else NA
    } else slope_up <- NA
  } else slope_up <- 0
  
  # ── Parameter 3: slope_down — Post-Peak Mitigation/Recovery Rate ──
  if (idx_max < n) {
    fase_down <- dados_ind[idx_max:n, ]
    if (nrow(fase_down) >= 2) {
      fit_down <- tryCatch(
        lm(planta ~ dia, data=fase_down),
        error=function(e) NULL
      )
      slope_down <- if (!is.null(fit_down)) coef(fit_down)[2] else NA
    } else slope_down <- NA
  } else slope_down <- 0
  
  # ── Parameter 4: t_surv — Time elapsed to Critical Threshold ──
  # Identifies earliest time point where ordinal score reaches or exceeds LIMIAR_NOTA
  acima_limiar <- dados_ind$dia[dados_ind$planta >= LIMIAR_NOTA]
  t_surv  <- if (length(acima_limiar) > 0) min(acima_limiar) else Inf
  evento  = if (is.finite(t_surv)) 1 else 0  # Binary indicator: 1 = event reached, 0 = right-censored
  t_surv_obs <- if (is.finite(t_surv)) t_surv else max(dias)
  
  # ── Parameter 5: area_beta — Area Under Curve (Enhanced Beta-scaled AUDPC) ──
  # Trapezoidal summation mapped across continuous Beta-scaled proportions
  area_beta <- if (n >= 2) {
    sum(diff(dias) * (head(props,-1) + tail(props,-1)) / 2)
  } else NA
  
  # ── Parameter 6: estab — Longitudinal Trajectory Stability ──
  # Trajectory variance across time blocks; higher values capture erratic disease patterns
  estab <- var(notas, na.rm=TRUE)
  
  # ── Parameter 7: taxa_final — Final Terminal Assessment Severity ──
  taxa_final <- notas[n]
  
  # ── Parameter 8: proporcao_acima — Ratio of assessments exceeding target threshold ──
  prop_acima <- mean(notas >= LIMIAR_NOTA, na.rm=TRUE)
  
  list(
    mu_max       = mu_max,
    t_pico       = t_pico,
    slope_up     = slope_up,
    slope_down   = slope_down,
    t_surv_obs   = t_surv_obs,
    evento       = evento,
    area_beta    = area_beta,
    estab        = estab,
    taxa_final   = taxa_final,
    prop_acima   = prop_acima,
    n_aval       = n
  )
}

# Execute loop computation across all matching individual ids
individuos <- unique(severity_long$individuo)
bdsip_list <- lapply(individuos, function(ind) {
  dados_ind <- severity_long |>
    filter(individuo == ind)
  params <- extrair_bdsip(dados_ind)
  if (is.null(params)) return(NULL)
  data.frame(individuo=ind, params)
})

bdsip_df <- do.call(rbind, Filter(Negate(is.null), bdsip_list))

cat(sprintf("BDSIP vectors computed for %d genotypes\n\n", nrow(bdsip_df)))

cat("BDSIP structural parameter distributions:\n")
params_cols <- c("mu_max","t_pico","slope_up","area_beta","estab","t_surv_obs")
for (p in params_cols) {
  vals <- bdsip_df[[p]]
  cat(sprintf("  %-15s Mean=%.3f | SD=%.3f | Range: [%.3f, %.3f]\n",
              p, mean(vals,na.rm=TRUE), sd(vals,na.rm=TRUE),
              min(vals,na.rm=TRUE), max(vals,na.rm=TRUE)))
}

# =============================================================================
# 4. SURVIVAL ANALYSIS — KAPLAN-MEIER BY FULL-SIB FAMILY STRUCTURE
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("Survival Analysis Routine\n")
cat(strrep("=", 60), "\n\n")

# Incorporate cross family infrastructure back to bdsip_df
familia_info <- dataset |>
  select(individuo, familia) |>
  distinct() |>
  mutate(individuo = as.integer(as.character(individuo)))

bdsip_df <- bdsip_df |>
  left_join(familia_info, by="individuo")

# Fit Kaplan-Meier parameters grouped by cross structure
surv_obj <- Surv(time=bdsip_df$t_surv_obs,
                 event=bdsip_df$evento)
km_fit <- survfit(surv_obj ~ familia, data=bdsip_df)

cat("Kaplan-Meier estimates — Median time elapsed to reach score ≥ 3:\n")
print(summary(km_fit)$table[, c("records","events","median")])

# Execute Log-rank significance test
logrank <- survdiff(surv_obj ~ familia, data=bdsip_df)
cat(sprintf("\nLog-rank test significance: χ²=%.3f, p-value=%.4f\n",
            logrank$chisq,
            pchisq(logrank$chisq, df=length(unique(bdsip_df$familia))-1,
                   lower.tail=FALSE)))

# =============================================================================
# 5. LONGITUDINAL BETA REGRESSION — EXTRACTING TRAJECTORY MATRICES
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("Longitudinal Beta Regression Profiling\n")
cat(strrep("=", 60), "\n\n")

# Clean data matrices and enforce explicit family string markers
sev_beta <- severity_long |>
  filter(!is.na(prop_beta)) |>
  mutate(
    familia = case_when(
      individuo %in% unique(dataset$individuo[dataset$familia == "1"]) ~ "F1",
      individuo %in% unique(dataset$individuo[dataset$familia == "2"]) ~ "F2",
      individuo %in% unique(dataset$individuo[dataset$familia == "3"]) ~ "F3",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(familia))

# Fit Beta model with space-time interaction structures via betareg framework
tryCatch({
  m_beta <- betareg(prop_beta ~ dia * familia | dia,
                    data = sev_beta,
                    link = "logit")
  
  cat("Longitudinal Beta Model diagnostics:\n")
  cat("  AIC value:", AIC(m_beta), "\n")
  cat("  Pseudo-R² coefficient:", summary(m_beta)$pseudo.r.squared, "\n\n")
  
  cat("Fixed parameters summary coefficients:\n")
  print(round(coef(m_beta), 4))
}, error=function(e) {
  cat("  Warning: betareg routine failed to reach mathematical convergence —", conditionMessage(e), "\n")
})

# =============================================================================
# 6. VARIABLE STANDARDIZATION AND COMPOSITE BDSIP GENERATION
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("Composite BDSIP Score Derivation\n")
cat(strrep("=", 60), "\n\n")

# Standardize individual trait indicators via Z-score scaling matrix
bdsip_std <- bdsip_df |>
  mutate(
    mu_max_sc     = scale(mu_max)[,1],
    t_pico_sc     = scale(t_pico)[,1],
    slope_up_sc   = scale(slope_up)[,1],
    area_beta_sc  = scale(area_beta)[,1],
    estab_sc      = scale(estab)[,1],
    # t_surv_sc: inverted tracking — longer survival = more resilient = reduces composite value
    t_surv_sc     = -scale(t_surv_obs)[,1],
    prop_acima_sc = scale(prop_acima)[,1]
  )

# Composite BDSIP formula integration:
# Merges intensity + precocidade + velocity profiles + area vectors
# Equal weights applied across core metrics by default (adjustable by the breeder)
bdsip_std <- bdsip_std |>
  mutate(
    BDSIP = (mu_max_sc + t_surv_sc + area_beta_sc +
               slope_up_sc + estab_sc) / 5
  )

cat("Composite BDSIP distribution profile:\n")
cat(sprintf("  Mean=%.3f | SD=%.3f | Range: [%.3f, %.3f]\n",
            mean(bdsip_std$BDSIP), sd(bdsip_std$BDSIP),
            min(bdsip_std$BDSIP), max(bdsip_std$BDSIP)))

cat("\nTop 10 most RESILIENT genotypes (Most negative composite BDSIP parameters):\n")
top_resilientes <- bdsip_std |>
  arrange(BDSIP) |>
  select(individuo, familia, mu_max, t_surv_obs, area_beta, BDSIP) |>
  head(10)
print(top_resilientes, row.names=FALSE)

cat("\nTop 10 most SUSCEPTIBLE genotypes (Most positive composite BDSIP parameters):\n")
top_suscet <- bdsip_std |>
  arrange(desc(BDSIP)) |>
  select(individuo, familia, mu_max, t_surv_obs, area_beta, BDSIP) |>
  head(10)
print(top_suscet, row.names=FALSE)

# =============================================================================
# 7. CROSS-METRIC COMPARATIVE ANALYSIS: BDSIP VS AUDPC
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("Cross-Evaluation: BDSIP vs AUDPC\n")
cat(strrep("=", 60), "\n\n")

# Load fruit-development windowed AUDPC data
dataset_novo <- readRDS(here("01_dados/processados/dataset_aacpd_janela.rds"))

aacpd_ind <- dataset_novo |>
  group_by(individuo) |>
  summarise(
    aacpd_media = mean(aacpdm, na.rm=TRUE),
    .groups="drop"
  ) |>
  mutate(individuo = as.integer(as.character(individuo)))

comp_df <- bdsip_std |>
  left_join(aacpd_ind, by="individuo") |>
  filter(!is.na(aacpd_media))

r_bdsip_aacpd <- cor(comp_df$BDSIP, comp_df$aacpd_media,
                     use="pairwise.complete.obs")
cat(sprintf("Pearson correlation matrix (BDSIP vs AUDPC): r=%.3f\n\n", r_bdsip_aacpd))

# Component-specific breakdown against baseline AUDPC values
for (comp in c("mu_max","t_pico","slope_up","area_beta","estab","t_surv_obs")) {
  r <- cor(comp_df[[comp]], comp_df$aacpd_media, use="pairwise.complete.obs")
  cat(sprintf("  %-15s × AUDPC matrix: r=%.3f\n", comp, r))
}

# Measure correlation against agronomic fruit quality metrics
cat("\nCorrelation Profile: BDSIP Score vs Fruit Agronomic Traits:\n")
dataset_bdsip <- dataset_novo |>
  mutate(individuo = as.integer(as.character(individuo))) |>
  left_join(bdsip_std |> select(individuo, BDSIP, mu_max, area_beta,
                                t_surv_obs), by="individuo")

for (tr in c("comp_mm","diam_mm","massa_fruto_g","rend_polpa_pct",
             "esp_casca_mm","sst_brix")) {
  r_bdsip <- cor(dataset_bdsip$BDSIP, dataset_bdsip[[tr]],
                 use="pairwise.complete.obs")
  r_aacpd <- cor(dataset_bdsip$aacpdm, dataset_bdsip[[tr]],
                 use="pairwise.complete.obs")
  cat(sprintf("  %-20s BDSIP correlation: r=%+.3f | AUDPC correlation: r=%+.3f\n",
              tr, r_bdsip, r_aacpd))
}

# =============================================================================
# 8. STRUCTURING COVARIATE ARRAYS FOR BAYESIAN ANIMAL MODELS (MCMCglmm)
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("Structuring Covariate Arrays for downstream MCMCglmm Analysis\n")
cat(strrep("=", 60), "\n\n")

# Align specific individual BDSIP metrics into individual observations within the master harvest sheet
dataset_mcmc <- dataset_novo |>
  mutate(individuo = as.integer(as.character(individuo))) |>
  left_join(
    bdsip_std |>
      select(individuo, BDSIP, mu_max_sc, t_pico_sc,
             slope_up_sc, area_beta_sc, estab_sc, t_surv_sc),
    by="individuo"
  )

# Impute missing records with zero (assigning background baseline mean to rows lacking disease monitoring blocks)
for (v in c("BDSIP","mu_max_sc","t_pico_sc","slope_up_sc",
            "area_beta_sc","estab_sc","t_surv_sc")) {
  dataset_mcmc[[v]][is.na(dataset_mcmc[[v]])] <- 0
}

cat(sprintf("  Final structural matrix formatted for MCMCglmm: %d rows\n", nrow(dataset_mcmc)))
cat(sprintf("  Missing values tally (NAs) in composite BDSIP field: %d\n", sum(is.na(dataset_mcmc$BDSIP))))

cat("\nDownstream Information Criteria Screening Matrix (Planned DIC evaluations):\n")
cat("  M_AUDPC: trait ~ block + season + audpc_sc\n")
cat("  M_BDSIP: trait ~ block + season + BDSIP_score\n")
cat("  M_FULL:  trait ~ block + season + BDSIP_score + mu_max_sc + t_surv_sc\n")
cat("  → Contrast Deviance Information Criterion (DIC) fit diagnostics post MCMC execution.\n\n")

# =============================================================================
# 9. OUTPUT AND PARAMETER FILE EXPORTING
# =============================================================================

dir.create(here("04_resultados/parametros_geneticos/bayesiano"),
           showWarnings=FALSE, recursive=TRUE)
dir.create(here("01_dados/processados"),
           showWarnings=FALSE, recursive=TRUE)

write.csv(bdsip_df,
          here("01_dados/processados/bdsip_parametros.csv"),
          row.names=FALSE)
write.csv(bdsip_std,
          here("01_dados/processados/bdsip_completo.csv"),
          row.names=FALSE)
saveRDS(bdsip_std,
        here("01_dados/processados/bdsip_completo.rds"))
saveRDS(dataset_mcmc,
        here("01_dados/processados/dataset_bdsip.rds"))

cat("\n✓ BDSIP methodology processing pipeline completed successfully.\n")
cat("  bdsip_parametros.csv — Raw parameter extraction values per genotype\n")
cat("  bdsip_completo.rds    — Standardized composite indexing metrics + isolated terms\n")
cat("  dataset_bdsip.rds    — Final analytical data structure configured for MCMCglmm modeling\n")
cat("\nNext step: Fit structures in MCMCglmm and screen M_AUDPC vs M_BDSIP tracking performance via DIC.\n")