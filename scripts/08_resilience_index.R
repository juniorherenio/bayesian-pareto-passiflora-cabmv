# =============================================================================
# Índice de Resiliência Bayesiano para Seleção de Genitores
# Contexto: RC1 (edulis × setacea) — recuperar características comerciais
#           de edulis mantendo resistência ao CABMV de setacea
#
# ESTRUTURA DO ÍNDICE:
#   IR_final = 0.60 × IR_comercial   (recuperar edulis)
#            + 0.40 × IR_resistencia (manter resistência setacea)
#            + 0.05 × IR_estabilidade (desempate)
#            - 0.10 × penalidade_precisao
#
# IR_comercial:
#   Traits que distinguem edulis de setacea
#   Pesos maiores para casca (herança mais marcante de setacea)
#   esp_casca_mm invertido (casca FINA = mais próximo de edulis)
#
# IR_resistencia:
#   P_severo_media do POM (invertido — baixa severidade = melhor)
#   + BDSIP t_surv (tempo até nota≥3 — maior = mais resiliente)
#
# IR_estabilidade:
#   Variância entre estações (invertida — baixa = mais estável)
#   IE = eficiência sob estresse (verão/inverno)
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(stringr)

cat("\n", strrep("=", 60), "\n")
cat("05_indice_resiliencia.R — v3\n")
cat("Contexto: RC1 (edulis × setacea)\n")
cat(format(Sys.time(), "[%H:%M:%S]"), "\n")
cat(strrep("=", 60), "\n\n")

# =============================================================================
# 1. CARREGAR DADOS
# =============================================================================

blups_est <- readRDS(here(
  "04_resultados/blups/bayesiano/blups_GxA_POM_por_estacao.rds"))
pom_ind   <- readRDS(here("01_dados/processados/pom_indices.rds"))
bdsip_std <- readRDS(here("01_dados/processados/bdsip_completo.rds"))
dataset   <- as.data.frame(readRDS(here(
  "01_dados/processados/dataset_pom.rds")))

cat(sprintf("BLUPs: %d obs (%d traits × %d ind × %d est)\n\n",
            nrow(blups_est),
            n_distinct(blups_est$trait),
            n_distinct(blups_est$individuo),
            n_distinct(blups_est$estacao)))

blups_est <- blups_est |>
  mutate(var_pred=blup_sd^2, precisao=1/var_pred)

# Info de família
familia_info <- dataset |>
  select(individuo, familia, cruzamento) |>
  distinct() |>
  mutate(individuo=as.character(individuo))

# =============================================================================
# 2. PARÂMETROS — CONTEXTO RC1
# =============================================================================

traits_7 <- c("comp_mm","diam_mm","massa_fruto_g",
              "rend_polpa_pct","esp_casca_mm","sst_brix","massa_polpa_g")

# h²_g ponderados do 04e_v2
h2_pond <- c(
  comp_mm        = 0.372,
  diam_mm        = 0.408,
  massa_fruto_g  = 0.372,
  rend_polpa_pct = 0.378,
  esp_casca_mm   = 0.404,
  sst_brix       = 0.391,
  massa_polpa_g  = 0.364
)

# Pesos comerciais — contexto RC1 → edulis
# esp_casca_mm tem peso maior: casca grossa é herança mais marcante de setacea
# sst_brix e massa_polpa_g têm peso menor: setacea já varia nessas
w_comercial <- c(
  comp_mm        = 0.15,  # fruto pequeno em setacea
  diam_mm        = 0.15,  # idem
  massa_fruto_g  = 0.15,  # fruto leve em setacea
  rend_polpa_pct = 0.20,  # polpa escassa em setacea
  esp_casca_mm   = 0.15,  # casca grossa — herança marcante de setacea
  sst_brix       = 0.10,  # menor peso
  massa_polpa_g  = 0.10   # menor peso
)
w_comercial <- w_comercial / sum(w_comercial)

# Direção desejável (maior = melhor, exceto esp_casca)
direcao <- c(
  comp_mm        =  1,
  diam_mm        =  1,
  massa_fruto_g  =  1,
  rend_polpa_pct =  1,
  esp_casca_mm   = -1,  # casca FINA = mais próximo de edulis
  sst_brix       =  1,
  massa_polpa_g  =  1
)

# Pesos do índice final
alpha_comercial   <- 0.60  # recuperar edulis
beta_resistencia  <- 0.40  # manter resistência setacea
gamma_estab       <- 0.05  # desempate
lambda_precisao   <- 0.10  # penalidade por incerteza

cat("Estrutura do índice:\n")
cat(sprintf("  α_comercial  = %.2f (recuperar edulis)\n", alpha_comercial))
cat(sprintf("  β_resistência = %.2f (resistência CABMV)\n", beta_resistencia))
cat(sprintf("  γ_estab      = %.2f (desempate)\n", gamma_estab))
cat(sprintf("  λ_precisão   = %.2f (penalidade incerteza)\n\n", lambda_precisao))

cat("Pesos comerciais (contexto RC1→edulis):\n")
print(round(w_comercial, 3))

# =============================================================================
# 3. COMPONENTE COMERCIAL — IR_comercial
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("Componente Comercial (RC1→edulis)\n")
cat(strrep("=", 60), "\n\n")

# Padronizar BLUPs por trait
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

# IR_comercial por indivíduo
IR_com <- blups_std |>
  group_by(individuo) |>
  summarise(
    IR_comercial   = sum(contribuicao, na.rm=TRUE) /
      sum(peso_combinado, na.rm=TRUE),
    var_pred_medio = mean(var_pred, na.rm=TRUE),
    precisao_media = mean(precisao, na.rm=TRUE),
    .groups = "drop"
  )

cat(sprintf("IR_comercial: média=%.3f | SD=%.3f | [%.3f, %.3f]\n\n",
            mean(IR_com$IR_comercial), sd(IR_com$IR_comercial),
            min(IR_com$IR_comercial), max(IR_com$IR_comercial)))

# =============================================================================
# 4. COMPONENTE DE RESISTÊNCIA — IR_resistencia
# =============================================================================

cat(strrep("=", 60), "\n")
cat("Componente de Resistência (CABMV)\n")
cat(strrep("=", 60), "\n\n")

# P_severo_media do POM — invertido (baixa severidade = melhor)
pom_info <- pom_ind |>
  select(individuo, P_severo_media, P_severo_pico, n_aval) |>
  mutate(individuo = as.character(individuo))

# BDSIP — t_surv_obs (maior = mais resiliente)
bdsip_info <- bdsip_std |>
  select(individuo, t_surv_obs, mu_max, area_beta) |>
  mutate(
    individuo    = as.character(individuo),
    t_surv_obs   = as.numeric(t_surv_obs),
    mu_max       = as.numeric(mu_max),
    area_beta    = as.numeric(area_beta)
  )

# Combinar POM + BDSIP para IR_resistencia
resist_df <- pom_info |>
  left_join(bdsip_info, by="individuo") |>
  mutate(
    # Padronizar — invertendo P_severo (baixo = resiliente)
    P_sev_sc    = -scale(P_severo_media)[,1],
    t_surv_sc   =  scale(t_surv_obs)[,1],
    mu_max_sc   = -scale(mu_max)[,1],  # pico baixo = melhor
    
    # IR_resistencia = média ponderada
    # POM tem peso maior (controlado para tempo e família)
    IR_resistencia = (0.5 * P_sev_sc +
                        0.3 * t_surv_sc +
                        0.2 * mu_max_sc)
  ) |>
  select(individuo, IR_resistencia, P_severo_media,
         t_surv_obs, P_sev_sc, t_surv_sc)

cat(sprintf("IR_resistencia calculado para %d indivíduos\n", nrow(resist_df)))
cat(sprintf("IR_resist: média=%.3f | SD=%.3f | [%.3f, %.3f]\n\n",
            mean(resist_df$IR_resistencia, na.rm=TRUE),
            sd(resist_df$IR_resistencia, na.rm=TRUE),
            min(resist_df$IR_resistencia, na.rm=TRUE),
            max(resist_df$IR_resistencia, na.rm=TRUE)))

# =============================================================================
# 5. COMPONENTE DE ESTABILIDADE — IR_estabilidade
# =============================================================================

cat(strrep("=", 60), "\n")
cat("Componente de Estabilidade temporal\n")
cat(strrep("=", 60), "\n\n")

est_baixa <- c("i18","i19")
est_alta  <- c("v19")

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
    # Inverter estab — baixa variância = mais estável = melhor
    estab_sc = -scale(estab_media)[,1],
    IE_sc    = ifelse(is.finite(IE_media),
                      scale(IE_media)[,1], 0),
    IR_estabilidade = 0.6 * estab_sc + 0.4 * IE_sc
  )

cat(sprintf("Estabilidade calculada para %d indivíduos\n", nrow(estab_ind)))
cat(sprintf("estab_media: [%.3f, %.3f]\n",
            min(estab_ind$estab_media), max(estab_ind$estab_media)))
cat(sprintf("IE válidos: %d de %d\n\n",
            sum(!is.na(estab_ind$IE_media)), nrow(estab_ind)))

# =============================================================================
# 6. ÍNDICE FINAL
# =============================================================================

cat(strrep("=", 60), "\n")
cat("ÍNDICE FINAL\n")
cat(strrep("=", 60), "\n\n")

# Juntar todos os componentes
idx_final <- IR_com |>
  mutate(individuo = as.character(individuo)) |>
  left_join(resist_df,  by="individuo") |>
  left_join(estab_ind |>
              mutate(individuo=as.character(individuo)) |>
              select(individuo, estab_media, IE_media,
                     IR_estabilidade),
            by="individuo") |>
  left_join(familia_info, by="individuo")

# Padronizar cada componente antes de combinar
idx_final <- idx_final |>
  mutate(
    IR_com_sc    = scale(IR_comercial)[,1],
    IR_res_sc    = ifelse(is.na(IR_resistencia), 0,
                          scale(IR_resistencia)[,1]),
    IR_est_sc    = ifelse(is.na(IR_estabilidade), 0,
                          IR_estabilidade),
    pen_prec_sc  = scale(var_pred_medio)[,1],
    
    # Índice final ponderado
    IR_final = alpha_comercial  * IR_com_sc  +
      beta_resistencia * IR_res_sc  +
      gamma_estab      * IR_est_sc  -
      lambda_precisao  * pen_prec_sc,
    
    rank_comercial   = rank(-IR_comercial),
    rank_resistencia = rank(-IR_res_sc),
    rank_final       = rank(-IR_final)
  )

cat(sprintf("Correlações entre componentes:\n"))
cat(sprintf("  IR_com × IR_resist: r=%.3f\n",
            cor(idx_final$IR_com_sc, idx_final$IR_res_sc,
                use="pairwise.complete.obs")))
cat(sprintf("  IR_com × IR_estab:  r=%.3f\n",
            cor(idx_final$IR_com_sc, idx_final$IR_est_sc,
                use="pairwise.complete.obs")))
cat(sprintf("  IR_com × IR_final:  r=%.3f\n",
            cor(idx_final$IR_com_sc, idx_final$IR_final,
                use="pairwise.complete.obs")))
cat(sprintf("  IR_res × IR_final:  r=%.3f\n\n",
            cor(idx_final$IR_res_sc, idx_final$IR_final,
                use="pairwise.complete.obs")))

# =============================================================================
# 7. RANKING FINAL
# =============================================================================

cat(strrep("=", 60), "\n")
cat("RANKING FINAL DE GENITORES\n")
cat(strrep("=", 60), "\n\n")

ranking <- idx_final |>
  arrange(rank_final) |>
  mutate(across(where(is.numeric), ~round(., 3))) |>
  select(individuo, familia, cruzamento,
         IR_comercial, IR_res_sc, IR_estabilidade,
         IR_final, rank_final,
         rank_comercial, rank_resistencia,
         P_severo_media, estab_media, var_pred_medio)

cat("Top 20 genitores:\n")
print(head(ranking[, c("individuo","familia",
                       "IR_comercial","IR_res_sc",
                       "IR_final","rank_final",
                       "P_severo_media","estab_media")], 20),
      row.names=FALSE)

cat("\nPor família:\n")
for (fam in c("1","2","3")) {
  sub <- ranking |> filter(as.character(familia) == fam)
  if (nrow(sub) > 0) {
    cat(sprintf("\n  Família %s (%s): %d indivíduos\n",
                fam, unique(sub$cruzamento)[1], nrow(sub)))
    cat(sprintf("    IR_comercial médio:   %.3f\n",
                mean(sub$IR_comercial, na.rm=TRUE)))
    cat(sprintf("    IR_resistência médio: %.3f\n",
                mean(sub$IR_res_sc, na.rm=TRUE)))
    cat(sprintf("    IR_final médio:       %.3f\n",
                mean(sub$IR_final, na.rm=TRUE)))
    cat(sprintf("    P_severo médio:       %.3f\n",
                mean(sub$P_severo_media, na.rm=TRUE)))
    cat(sprintf("    Top 5: ind %s\n",
                paste(head(sub$individuo, 5), collapse=", ")))
  }
}

# Indivíduos em conflito (alto comercial, baixa resistência e vice-versa)
cat("\n\nConflito comercial × resistência (diferença de rank > 20):\n")
conflito <- ranking |>
  mutate(delta_rank = abs(rank_comercial - rank_resistencia)) |>
  filter(delta_rank > 20) |>
  arrange(desc(delta_rank)) |>
  select(individuo, familia, rank_comercial, rank_resistencia,
         delta_rank, IR_comercial, IR_res_sc, P_severo_media)
if (nrow(conflito) > 0) {
  print(conflito, row.names=FALSE)
} else {
  cat("  Nenhum conflito expressivo encontrado.\n")
}

# =============================================================================
# 8. CONTRIBUIÇÃO POR TRAIT
# =============================================================================

cat("\nContribuição por trait no IR_comercial:\n")
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
# 9. SALVAR
# =============================================================================

dir.create(here("04_resultados/selecao"), showWarnings=FALSE, recursive=TRUE)

saveRDS(idx_final,
        here("04_resultados/selecao/indice_resiliencia.rds"))
write.csv(ranking,
          here("04_resultados/selecao/ranking_genitores.csv"),
          row.names=FALSE)
saveRDS(blups_std,
        here("04_resultados/selecao/blups_ponderados.rds"))

# Salvar componentes separados para NSGA-II
componentes <- idx_final |>
  select(individuo, familia, cruzamento,
         IR_comercial, IR_res_sc, IR_est_sc,
         IR_final, rank_final,
         P_severo_media, estab_media,
         var_pred_medio, precisao_media)
write.csv(componentes,
          here("04_resultados/selecao/componentes_indice.csv"),
          row.names=FALSE)

cat(sprintf("\n✓ 05 v3 concluído.\n"))
cat("  ranking_genitores.csv  — ranking final\n")
cat("  componentes_indice.csv — componentes separados\n")
cat("  indice_resiliencia.rds — índice completo\n")
cat("\nPróximo: 06_nsga2_pareto.R\n")