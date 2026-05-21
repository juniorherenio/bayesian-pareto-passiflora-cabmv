# =============================================================================
# Bayesian Resilience Index for Parent Selection
# Breeding Context: BC1 (Passiflora edulis × Passiflora setacea) 
#   Objective: Recover the commercial fruit traits of P. edulis 
#   while maintaining the CABMV resistance introgressed from P. setacea.
#
# INDEX STRUCTURE SPECIFICATION:
#   IR_final = 0.60 × IR_comercial   (Target: P. edulis genome recovery)
#            + 0.40 × IR_resistencia (Target: Maintain P. setacea resistance)
#            + 0.05 × IR_estabilidade (Tie-breaking stability parameter)
#            - 0.10 × penalty_precision
#
# IR_comercial Component:
#   Agronomic traits that differentiate P. edulis from P. setacea.
#   Higher relative weights assigned to rind thickness (the most persistent 
#   wild-type maternal inheritance from P. setacea).
#   'esp_casca_mm' is inverted (THIN rind = approaches target P. edulis ideotype).
#
# IR_resistencia Component:
#   POM-derived 'P_severo_media' (Inverted — low severe probability = optimal)
#   + BDSIP-derived 't_surv_obs' (Survival time to score ≥ 3 — longer = superior)
#
# IR_estabilidade Component:
#   Inter-season variance (Inverted — low longitudinal variance = stable)
#   IE parameter = Stress tolerance efficiency ratio (Summer vs Winter variance)
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(stringr)

cat("\n", strrep("=", 60), "\n")
cat("05_resilience_index.R — v3\n")
cat("Breeding Context: BC1 (P. edulis × P. setacea)\n")
cat(format(Sys.time(), "[%H:%M:%S]"), "\n")
cat(strrep("=", 60), "\n\n")

# =============================================================================
# 1. DATA LOADING AND METADATA ALIGNMENT
# =============================================================================

blups_est <- readRDS(here(
  "04_resultados/blups/bayesiano/blups_GxA_POM_por_estacao.rds"))
pom_ind   <- readRDS(here("01_dados/processados/pom_indices.rds"))
bdsip_std <- readRDS(here("01_dados/processados/bdsip_completo.rds"))
dataset   <- as.data.frame(readRDS(here(
  "01_dados/processados/dataset_pom.rds")))

cat(sprintf("Estimated Breeding Values Loaded: %d records (%d traits × %d genotypes × %d seasons)\n\n",
            nrow(blups_est),
            n_distinct(blups_est$trait),
            n_distinct(blups_est$individuo),
            n_distinct(blups_est$estacao)))

blups_est <- blups_est |>
  mutate(var_pred=blup_sd^2, precisao=1/var_pred)

# Cross structure and full-sib family mapping
familia_info <- dataset |>
  select(individuo, familia, cruzamento) |>
  distinct() |>
  mutate(individuo=as.character(individuo))

# =============================================================================
# 2. SELECTION PARAMETERS — BACKCROSSING (BC1) DEFINITIONS
# =============================================================================

traits_7 <- c("comp_mm","diam_mm","massa_fruto_g",
              "rend_polpa_pct","esp_casca_mm","sst_brix","massa_polpa_g")

# Weighted narrow-sense heritabilities (h²) derived from the univariated 04e_v2 script
h2_pond <- c(
  comp_mm        = 0.372,
  diam_mm        = 0.408,
  massa_fruto_g  = 0.372,
  rend_polpa_pct = 0.378,
  esp_casca_mm   = 0.404,
  sst_brix       = 0.391,
  massa_polpa_g  = 0.364
)

# Commercial weights matrix — Backcross context (BC1 → P. edulis recurrent parent)
# 'esp_casca_mm' gets a high relative weight: thick rind is a strong wild trait linkage
# 'sst_brix' and 'massa_polpa_g' have lower relative weight adjustments
w_comercial <- c(
  comp_mm        = 0.15,  # Small fruit size found in wild P. setacea
  diam_mm        = 0.15,  # Idem
  massa_fruto_g  = 0.15,  # Reduced fruit weight found in wild P. setacea
  rend_polpa_pct = 0.20,  # Extremely scarce pulp ratio in wild P. setacea
  esp_casca_mm   = 0.15,  # Thick rind phenotype — prominent wild donor trait linkage
  sst_brix       = 0.10,  # Lower relative weight adjustment
  massa_polpa_g  = 0.10   # Lower relative weight adjustment
)
w_comercial <- w_comercial / sum(w_comercial)

# Desirable selection direction matrix (Positive = increase, negative = decrease)
direcao <- c(
  comp_mm        =  1,
  diam_mm        =  1,
  massa_fruto_g  =  1,
  rend_polpa_pct =  1,
  esp_casca_mm   = -1,  # THIN rind = approaches target P. edulis recurrent phenotype
  sst_brix       =  1,
  massa_polpa_g  =  1
)

# Multi-trait index integration coefficients
alpha_comercial   <- 0.60  # P. edulis recurrent genome recovery weight
beta_resistencia  <- 0.40  # Introgression donor resistance retention weight
gamma_estab       <- 0.05  # Tie-breaking environmental buffer coefficient
lambda_precisao   <- 0.10  # Posterior uncertainty prediction error penalty

cat("Resilience Index Structural Coefficients Configuration:\n")
cat(sprintf("  α_comercial   = %.2f (Recurrent genome recovery)\n", alpha_comercial))
cat(sprintf("  β_resistencia = %.2f (CABMV resistance matrix retention)\n", beta_resistencia))
cat(sprintf("  γ_estab       = %.2f (Tie-breaking stability adjustment)\n", gamma_estab))
cat(sprintf("  λ_precisao    = %.2f (Uncertainty prediction error penalty)\n\n", lambda_precisao))

cat("Relative Commercial Weight Matrix (BC1 Strategy Optimization):\n")
print(round(w_comercial, 3))

# =============================================================================
# 3. COMMERCIAL PHENOTYPE COMPONENT — IR_comercial
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("Commercial Phenotype Integration Component (BC1 → P. edulis)\n")
cat(strrep("=", 60), "\n\n")

# Z-score standardization across independent breeding values
blups_std <- blups_est |>
  group_by(trait) |>
  mutate(
    blup_z   = (blup - mean(blup, na.rm=TRUE)) / sd(blup, na.rm=TRUE),
    blup_dir = blup_z * direcao[trait]
  ) |>
  ungroup() |>
  mutate(
    peso_combinado = h2_pond[trait] * w_comercial[trait] * precisao,
    contribuicao   = blup_dir * peso_combinado
  )

# Collapse to individual-level IR_comercial indices
IR_com <- blups_std |>
  group_by(individuo) |>
  summarise(
    IR_comercial    = sum(contribuicao, na.rm=TRUE) /
      sum(peso_combinado, na.rm=TRUE),
    var_pred_medio = mean(var_pred, na.rm=TRUE),
    precisao_media = mean(precisao, na.rm=TRUE),
    .groups = "drop"
  )

cat(sprintf("IR_comercial Summary: Mean=%.3f | SD=%.3f | Range: [%.3f, %.3f]\n\n",
            mean(IR_com$IR_comercial), sd(IR_com$IR_comercial),
            min(IR_com$IR_comercial), max(IR_com$IR_comercial)))

# =============================================================================
# 4. PATHOSYSTEM RESISTANCE COMPONENT — IR_resistencia
# =============================================================================

cat(strrep("=", 60), "\n")
cat("Pathosystem Resistance Retention Component (CABMV Validation)\n")
cat(strrep("=", 60), "\n\n")

# Extract POM-derived 'P_severo_media' — inverted tracking (low values = resilient)
pom_info <- pom_ind |>
  select(individuo, P_severo_media, P_severo_pico, n_aval) |>
  mutate(individuo = as.character(individuo))

# Extract BDSIP-derived 't_surv_obs' — positive tracking (longer survival = optimal)
bdsip_info <- bdsip_std |>
  select(individuo, t_surv_obs, mu_max, area_beta) |>
  mutate(
    individuo    = as.character(individuo),
    t_surv_obs   = as.numeric(t_surv_obs),
    mu_max       = as.numeric(mu_max),
    area_beta    = as.numeric(area_beta)
  )

# Integrate POM and BDSIP metrics to construct 'IR_resistencia'
resist_df <- pom_info |>
  left_join(bdsip_info, by="individuo") |>
  mutate(
    # Z-score standardization: Inverting P_severo (low probability values indicate resilience)
    P_sev_sc    = -scale(P_severo_media)[,1],
    t_surv_sc   =  scale(t_surv_obs)[,1],
    mu_max_sc   = -scale(mu_max)[,1],  # Suppressed epidemic peak is preferred
    
    # Weighted calculation: POM parameters get higher weights due to time-space adjustments
    IR_resistencia = (0.5 * P_sev_sc +
                        0.3 * t_surv_sc +
                        0.2 * mu_max_sc)
  ) |>
  select(individuo, IR_resistencia, P_severo_media,
         t_surv_obs, P_sev_sc, t_surv_sc)

cat(sprintf("IR_resistencia arrays successfully mapped for %d genotypes\n", nrow(resist_df)))
cat(sprintf("IR_resistencia Summary: Mean=%.3f | SD=%.3f | Range: [%.3f, %.3f]\n\n",
            mean(resist_df$IR_resistencia, na.rm=TRUE),
            sd(resist_df$IR_resistencia, na.rm=TRUE),
            min(resist_df$IR_resistencia, na.rm=TRUE),
            max(resist_df$IR_resistencia, na.rm=TRUE)))

# =============================================================================
# 5. ENVIRONMENTAL CLIMATIC BUFFERING COMPONENT — IR_estabilidade
# =============================================================================

cat(strrep("=", 60), "\n")
cat("Longitudinal Environmental Buffering Component (Temporal Stability)\n")
cat(strrep("=", 60), "\n\n")

est_baixa <- c("i18","i19") # Favorable, low stress winter seasons
est_alta  <- c("v19")       # Restrictive, high disease pressure summer season

estab_ind <- blups_std |>
  group_by(individuo, trait) |>
  summarise(
    blup_var_est = var(blup, na.rm=TRUE),
    blup_baixa   = mean(blup[estacao %in% est_baixa], na.rm=TRUE),
    blup_alta    = mean(blup[estacao %in% est_alta],  na.rm=TRUE),
    .groups = "drop"
  ) |>
  mutate(
    IE = ifelse(is.finite(blup_baixa) & blup_baixa != 0,
                blup_alta / blup_baixa, NA_real_)
  ) |>
  group_by(individuo) |>
  summarise(
    estab_media = mean(blup_var_est, na.rm=TRUE),
    IE_media    = mean(IE, na.rm=TRUE),
    .groups = "drop"
  ) |>
  mutate(
    # Invert stability metric: Low tracking variance means high temporal stability buffer
    estab_sc = -scale(estab_media)[,1],
    IE_sc    = ifelse(is.finite(IE_media),
                      scale(IE_media)[,1], 0),
    IR_estabilidade = 0.6 * stab_sc + 0.4 * IE_sc
  )

cat(sprintf("Temporal stability parameters compiled for %d genotypes\n", nrow(estab_ind)))
cat(sprintf("Raw variance limits (estab_media): [%.3f, %.3f]\n",
            min(estab_ind$estab_media), max(estab_ind$estab_media)))
cat(sprintf("Valid stress-tolerance efficiency ratios (IE): %d out of %d genotypes\n\n",
            sum(!is.na(estab_ind$IE_media)), nrow(estab_ind)))

# =============================================================================
# 6. INTEGRATED RESILIENCE INDEX FORMULATION
# =============================================================================

cat(strrep("=", 60), "\n")
cat("INTEGRATED MASTER RESILIENCE INDEX MATRIX\n")
cat(strrep("=", 60), "\n\n")

# Join all isolated selection sub-components
idx_final <- IR_com |>
  mutate(individuo = as.character(individuo)) |>
  left_join(resist_df,  by="individuo") |>
  left_join(estab_ind |>
              mutate(individuo=as.character(individuo)) |>
              select(individuo, estab_media, IE_media,
                     IR_estabilidade),
            by="individuo") |>
  left_join(familia_info, by="individuo")

# Standardize localized metrics prior to deploying the index combination equation
idx_final <- idx_final |>
  mutate(
    IR_com_sc    = scale(IR_comercial)[,1],
    IR_res_sc    = ifelse(is.na(IR_resistencia), 0,
                          scale(IR_resistencia)[,1]),
    IR_est_sc    = ifelse(is.na(IR_estabilidade), 0,
                          IR_estabilidade),
    pen_prec_sc  = scale(var_pred_medio)[,1],
    
    # Integrated Multi-Trait Resilience Equation Calculation
    IR_final = alpha_comercial   * IR_com_sc  +
      beta_resistencia * IR_res_sc  +
      gamma_estab      * IR_est_sc  -
      lambda_precisao  * pen_prec_sc,
    
    rank_comercial   = rank(-IR_comercial),
    rank_resistencia = rank(-IR_res_sc),
    rank_final       = rank(-IR_final)
  )

cat(sprintf("Cross-Component Pearson Correlation Profiles:\n"))
cat(sprintf("  IR_com × IR_res resistencia: r=%.3f\n",
            cor(idx_final$IR_com_sc, idx_final$IR_res_sc,
                use="pairwise.complete.obs")))
cat(sprintf("  IR_com × IR_est estab:        r=%.3f\n",
            cor(idx_final$IR_com_sc, idx_final$IR_est_sc,
                use="pairwise.complete.obs")))
cat(sprintf("  IR_com × IR_final master:     r=%.3f\n",
            cor(idx_final$IR_com_sc, idx_final$IR_final,
                use="pairwise.complete.obs")))
cat(sprintf("  IR_res × IR_final master:     r=%.3f\n\n",
            cor(idx_final$IR_res_sc, idx_final$IR_final,
                use="pairwise.complete.obs")))

# =============================================================================
# 7. MASTER BREEDING SELECTION COHORT RANKING
# =============================================================================

cat(strrep("=", 60), "\n")
cat("FINAL BREEDING COHORT PARENT RANKING\n")
cat(strrep("=", 60), "\n\n")

ranking <- idx_final |>
  arrange(rank_final) |>
  mutate(across(where(is.numeric), ~round(., 3))) |>
  select(individuo, familia, cruzamento,
         IR_comercial, IR_res_sc, IR_estabilidade,
         IR_final, rank_final,
         rank_comercial, rank_resistencia,
         P_severo_media, estab_media, var_pred_medio)

cat("Top 20 Recommended Elite Parents Summary Cohort:\n")
print(head(ranking[, c("individuo","familia",
                       "IR_comercial","IR_res_sc",
                       "IR_final","rank_final",
                       "P_severo_media","estab_media")], 20),
      row.names=FALSE)

cat("\nPerformance Analysis Grouped by Full-Sib Family Structures:\n")
for (fam in c("1","2","3")) {
  sub <- ranking |> filter(as.character(familia) == fam)
  if (nrow(sub) > 0) {
    cat(sprintf("\n  Family Cross Cohort %s (%s): %d genotypes mapped\n",
                fam, unique(sub$cruzamento)[1], nrow(sub)))
    cat(sprintf("     Mean Commercial Component Index:       %.3f\n",
                mean(sub$IR_comercial, na.rm=TRUE)))
    cat(sprintf("     Mean Pathosystem Resistance Component:  %.3f\n",
                mean(sub$IR_res_sc, na.rm=TRUE)))
    cat(sprintf("     Mean Integrated Resilience Score:       %.3f\n",
                mean(sub$IR_final, na.rm=TRUE)))
    cat(sprintf("     Mean Pathosystem Severity Profile:      %.3f\n",
                mean(sub$P_severo_media, na.rm=TRUE)))
    cat(sprintf("     Top 5 Elite Genotypes Checklist: ind %s\n",
                paste(head(sub$individuo, 5), collapse=", ")))
  }
}

# Identify phenotypic antagonistic outliers (High fruit ideotype but highly susceptible, or vice-versa)
cat("\n\nAntagonistic Phenotypic Outliers (Rank Divergence matrix delta > 20):\n")
conflito <- ranking |>
  mutate(delta_rank = abs(rank_comercial - rank_resistencia)) |>
  filter(delta_rank > 20) |>
  arrange(desc(delta_rank)) |>
  select(individuo, familia, rank_comercial, rank_resistencia,
         delta_rank, IR_comercial, IR_res_sc, P_severo_media)
if (nrow(conflito) > 0) {
  print(conflito, row.names=FALSE)
} else {
  cat("  No extreme antagonistic trait linkage outliers detected in this population.\n")
}

# =============================================================================
# 8. STRUCTURAL RELATIVE CONTRIBUTION BREAKDOWN BY TRAIT
# =============================================================================

cat("\nStructural Matrix Relative Contribution Breakdown within IR_comercial:\n")
contrib_trait <- blups_std |>
  group_by(trait) |>
  summarise(
    contrib_media = round(mean(contribuicao, na.rm=TRUE), 5),
    contrib_sd    = round(sd(contribuicao, na.rm=TRUE), 5),
    peso_medio    = round(mean(peso_combinado, na.rm=TRUE), 5),
    .groups = "drop"
  ) |>
  arrange(desc(abs(contrib_media)))
print(contrib_trait, row.names=FALSE)

# =============================================================================
# 9. OUTPUT EXPORT AND CHECKPOINT PREPARATION FOR MULTI-OBJECTIVE MODELING
# =============================================================================

dir.create(here("04_resultados/selecao"), showWarnings=FALSE, recursive=TRUE)

saveRDS(idx_final,
        here("04_resultados/selecao/indice_resiliencia.rds"))
write.csv(ranking,
          here("04_resultados/selecao/ranking_genitores.csv"),
          row.names=FALSE)
saveRDS(blups_std,
        here("04_resultados/selecao/blups_ponderados.rds"))

# Export independent matrix variables optimized for Multi-Objective Pareto Modeling (NSGA-II)
componentes <- idx_final |>
  select(individuo, familia, cruzamento,
         IR_comercial, IR_res_sc, IR_est_sc,
         IR_final, rank_final,
         P_severo_media, estab_media,
         var_pred_medio, precisao_media)
write.csv(componentes,
          here("04_resultados/selecao/componentes_indice.csv"),
          row.names=FALSE)

cat(sprintf("\n✓ Module 05 execution complete. Selection targets established.\n"))
cat("  Exported Cohort Matrix:    ranking_genitores.csv   — Master selection layout table\n")
cat("  Exported Pareto Elements:  componentes_indice.csv  — De-coupled multi-objective indices\n")
cat("  Exported Master RDS:       indice_resiliencia.rds  — Integrated structural index\n")
cat("\nNext pipeline module execution path: 06_nsga2_pareto.R\n")