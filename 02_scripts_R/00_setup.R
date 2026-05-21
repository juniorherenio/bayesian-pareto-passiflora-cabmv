# =============================================================================
# Global project configuration
# Passion fruit: longitudinal phenomics + asreml-R + AI
# =============================================================================

library(here)
library(readxl)
library(asreml)
library(dplyr)
library(nasapower)
library(EnvRtype)

# Explicitly resolve namespace conflicts
filter <- dplyr::filter
select <- dplyr::select
lag    <- dplyr::lag

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
path_raw    <- here("01_dados/brutos")
path_env    <- here("01_dados/ambientais")
path_proc   <- here("01_dados/processados")
path_results <- here("04_resultados")
path_figures <- here("05_figuras")

# -----------------------------------------------------------------------------
# Dataset loading and preparation
# -----------------------------------------------------------------------------
dataset <- read_excel(here("01_dados/brutos/dataset.xlsx"),
                      col_types = c(rep("text", 10), rep("numeric", 15))) |>
  mutate(data = as.Date(as.numeric(data), origin = "1899-12-30")) |>
  rename(individuo = planta) |>
  mutate(parcela = paste0("F", familia, "B", bloco)) |>
  mutate(across(c(individuo, bloco, parcela, familia,
                  P1, P2, estacao, ano, estacao_ano), as.factor)) |>
  mutate(aacpdm = (aacpd_planta + aacpd_folha) / 2)

# Chronological order for season_year (estacao_ano) — includes summer_2019 (verao_2019)
niveis_ordem <- c("inverno_2018", "primavera_2018", "verao_2018",
                  "verao_2019",   "outono_2019",    "inverno_2019")
niveis_presentes <- niveis_ordem[niveis_ordem %in% levels(dataset$estacao_ano)]
dataset$estacao_ano <- factor(dataset$estacao_ano, levels = niveis_presentes)

# -----------------------------------------------------------------------------
# Response variables and covariates
# -----------------------------------------------------------------------------
vars_continuas   <- c("comp_mm", "diam_mm", "ind_formato",
                      "massa_fruto_g", "massa_polpa_g", "rend_polpa_pct",
                      "esp_casca_mm", "esp_casca_inv", "sst_brix")

vars_contagem    <- c("n_frutos")
vars_binarias    <- c("formato_elipsoide", "casca_amarela", "polpa_amarela")
vars_covariaveis <- c("aacpdm")

# -----------------------------------------------------------------------------
# Environment Confirmation
# -----------------------------------------------------------------------------
message("========================================")
message("Setup OK")
message("asreml  : ", packageVersion("asreml"))
message("dplyr   : ", packageVersion("dplyr"))
message("Dataset : ", nrow(dataset), " obs x ", ncol(dataset), " columns")
message("aacpdm  : ", sum(!is.na(dataset$aacpdm)), " valid obs / ",
        sum(is.na(dataset$aacpdm)), " NAs")
message("Season-year levels: ",
        paste(levels(dataset$estacao_ano), collapse = " | "))
message("ind 240 comp_mm: ",
        paste(round(dataset$comp_mm[dataset$individuo == "240"], 1), collapse = ", "))
message("========================================")