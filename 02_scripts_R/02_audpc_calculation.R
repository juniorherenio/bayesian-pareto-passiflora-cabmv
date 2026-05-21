# =============================================================================
# AUDPC recalculation using the fruit development window (75 days)
#
# STRATEGY:
#   For each fruit harvest (individuo × data_coleta):
#   → Identify disease severity assessments within the preceding 75 days
#   → Compute trapezoidal AUDPC within this specific window
#   → Use whole-plant severity scores (systemic state affecting pericarp)
#
# THREE CALCULATED VERSIONS:
#   audpc_acum    — Accumulated from the beginning (original version)
#   audpc_periodo — Reset at each harvest (period between consecutive harvests)
#   audpc_jan75   — 75-day moving window prior to harvest (RECOMMENDED)
#
# BIOLOGICAL RATIONALE:
#   Passion fruit: ~60-90 days from anthesis (flowering) to fruit maturation
#   Rind thickness (a key symptom of CABMV) develops during this time frame
#   The relevant disease pressure is concentrated within these 75 development days
# =============================================================================

library(here)
library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)

cat("\n", strrep("=", 60), "\n")
cat("AUDPC Recalculation — Fruit Development Window\n")
cat(format(Sys.time(), "[%H:%M:%S]"), "\n")
cat(strrep("=", 60), "\n\n")

# =============================================================================
# 1. LOAD DISEASE SEVERITY WORKSHEET
# =============================================================================

cat("Loading disease severity worksheet...\n")

xlsx_path <- here("01_dados/brutos/Severidade_CABMV.xlsx")

# Check if file exists within the project directory
if (!file.exists(xlsx_path)) {
  # Fallback to local absolute path
  xlsx_path <- "C:/Users/junio/OneDrive/Documentos/Artigos/2026 - Projeto Maracujá/01_dados/brutos/Severidade_CABMV.xlsx"
  cat("  File not found in project directory — using fallback absolute path\n")
}

# Read assessment headers (row 1) and assessment dates (row 2)
header1 <- read_excel(xlsx_path, sheet = "Severidade CABMV",
                      col_names = FALSE, n_max = 1)
header2 <- read_excel(xlsx_path, sheet = "Severidade CABMV",
                      col_names = FALSE, skip = 1, n_max = 1)

# Extract column indices and dates for disease assessments
n_cols  <- ncol(header1)
aval_cols <- which(!is.na(as.character(header1[1,])) &
                     grepl("Avalia", as.character(header1[1,])))

cat(sprintf("  %d assessments found\n", length(aval_cols)))

# Parse assessment dates
datas_aval <- sapply(aval_cols, function(col) {
  data_str <- gsub("Data\\s*:?\\s*", "", as.character(header2[1, col]))
  data_str <- trimws(data_str)
  
  # Try different date formats
  d <- NA
  for (fmt in c("%d/%m/%Y", "%d/%m/%y", "%Y-%m-%d")) {
    d <- tryCatch(as.Date(data_str, fmt), error=function(e) NA)
    if (!is.na(d)) break
  }
  
  # Return numeric format for sapply compatibility
  as.numeric(d)
})
datas_aval <- as.Date(datas_aval, origin="1970-01-01")

# Verify missing values (NAs) prior to typo corrections
cat("NAs in assessment dates:", sum(is.na(datas_aval)), "\n")
cat("Parsed dates:\n")
print(data.frame(col=aval_cols, data=datas_aval))

# Correct typographical year error only where date is not NA
idx_2017 <- which(!is.na(datas_aval) & year(datas_aval) == 2017)
if (length(idx_2017) > 0) {
  datas_aval[idx_2017] <- datas_aval[idx_2017] %m+% years(1)
  cat(sprintf("  %d dates corrected (2017→2018)\n", length(idx_2017)))
}

# Remove assessments with unparseable dates (NA)
ok_aval <- !is.na(datas_aval)
aval_cols  <- aval_cols[ok_aval]
datas_aval <- datas_aval[ok_aval]

cat(sprintf("  %d valid assessments (out of %d identified)\n",
            length(aval_cols), sum(!is.na(datas_aval)) + sum(is.na(datas_aval))))

# Correct remaining layout typos (2017 → 2018)
datas_aval[year(datas_aval) == 2017] <-
  datas_aval[year(datas_aval) == 2017] %m+% years(1)

cat(sprintf("  Experimental period: %s to %s\n",
            format(min(datas_aval, na.rm=TRUE), "%d/%m/%Y"),
            format(max(datas_aval, na.rm=TRUE), "%d/%m/%Y")))

# Load raw data (skipping headers, from row 4 onwards)
dados_raw <- read_excel(xlsx_path, sheet = "Severidade CABMV",
                        col_names = FALSE, skip = 3)

# =============================================================================
# 2. CONVERT TO LONG FORMAT
# =============================================================================

cat("\nConverting dataset to long format...\n")

severity_long <- data.frame()

for (k in seq_along(aval_cols)) {
  col_i <- aval_cols[k]
  data_k <- datas_aval[k]
  
  if (is.na(data_k)) next
  
  # Map columns: individuo (col_i), folha (col_i+1), planta (col_i+2)
  chunk <- dados_raw[, c(col_i, col_i+1, col_i+2)]
  names(chunk) <- c("individuo", "folha", "planta")
  
  chunk <- chunk |>
    mutate(
      data_aval  = data_k,
      avaliacao  = k,
      individuo  = suppressWarnings(as.integer(individuo)),
      folha      = suppressWarnings(as.numeric(
        ifelse(folha == "-" | is.na(folha), NA, folha))),
      planta     = suppressWarnings(as.numeric(
        ifelse(planta == "-" | is.na(planta), NA, planta)))
    ) |>
    filter(!is.na(individuo))
  
  severity_long <- bind_rows(severity_long, chunk)
}

severity_long <- severity_long |>
  arrange(individuo, data_aval) |>
  mutate(
    # Mean score of leaf and whole plant (consistent with the original paper methodology)
    nota_media = rowMeans(cbind(folha, planta), na.rm=TRUE)
  )

cat(sprintf("  Shape: %d obs | %d individuals | %d unique assessment dates\n",
            nrow(severity_long),
            n_distinct(severity_long$individuo),
            n_distinct(severity_long$data_aval)))

# =============================================================================
# 3. AUDPC CALCULATION — TRAPEZOIDAL FUNCTION
# =============================================================================

aacpd_trapezio <- function(notas, datas) {
  # Filter missing values
  ok    <- !is.na(notas) & !is.na(datas)
  notas <- notas[ok]
  datas <- as.numeric(datas[ok])
  
  if (length(notas) < 2) return(NA_real_)
  
  # Chronological sorting
  ord   <- order(datas)
  notas <- notas[ord]
  datas <- datas[ord]
  
  # Trapezoidal summation
  sum(diff(datas) * (head(notas,-1) + tail(notas,-1)) / 2)
}

# =============================================================================
# 4. LOAD FRUIT HARVEST DATASET
# =============================================================================

cat("\nLoading fruit harvest dataset...\n")
dataset <- readRDS(here("01_dados/processados/dataset_completo.rds"))
dataset <- as.data.frame(dataset)

# Convert harvest date to Date class
dataset$data_coleta <- as.Date(dataset$data,
                               origin = "1899-12-30")  # Excel date origin

cat(sprintf("  %d harvests | %d individuals | %d unique harvest dates\n",
            nrow(dataset),
            n_distinct(dataset$individuo),
            n_distinct(dataset$data_coleta)))

cat("\nHarvest dates sample:\n")
print(sort(unique(dataset$data_coleta))[1:10])

# =============================================================================
# 5. COMPUTE THE THREE AUDPC VERSIONS PER HARVEST EVENT
# =============================================================================

cat("\nComputing AUDPC versions per individual × harvest date...\n")

JANELA <- 75  # fruit development window in days

resultado <- dataset |>
  select(individuo, data_coleta) |>
  distinct() |>
  mutate(individuo = as.integer(as.character(individuo)))

resultado$aacpd_acum    <- NA_real_
resultado$aacpd_periodo <- NA_real_
resultado$aacpd_jan75   <- NA_real_

# Experiment onset date (first available disease assessment)
data_inicio <- min(datas_aval, na.rm=TRUE)

for (i in 1:nrow(resultado)) {
  ind    <- resultado$individuo[i]
  d_col  <- resultado$data_coleta[i]
  
  if (is.na(d_col)) next
  
  # Filter disease assessments for the target individual
  sev_ind <- severity_long |>
    filter(individuo == ind) |>
    arrange(data_aval)
  
  if (nrow(sev_ind) == 0) next
  
  # Version A — Accumulated AUDPC from baseline experiment onset
  sev_acum <- sev_ind |>
    filter(data_aval <= d_col)
  resultado$aacpd_acum[i] <- aacpd_trapezio(
    sev_acum$nota_media, sev_acum$data_aval)
  
  # Version B — Period-specific AUDPC since the previous harvest event
  # (identifies the immediate preceding harvest date for the individual)
  coletas_ind <- dataset |>
    filter(as.integer(as.character(individuo)) == ind,
           data_coleta < d_col) |>
    pull(data_coleta) |>
    unique() |>
    sort()
  
  d_prev <- if (length(coletas_ind) > 0) max(coletas_ind) else data_inicio
  
  sev_periodo <- sev_ind |>
    filter(data_aval > d_prev & data_aval <= d_col)
  resultado$aacpd_periodo[i] <- aacpd_trapezio(
    sev_periodo$nota_media, sev_periodo$data_aval)
  
  # Version C — 75-day moving window before harvest (RECOMMENDED)
  d_janela_inicio <- d_col - JANELA
  sev_jan <- sev_ind |>
    filter(data_aval >= d_janela_inicio & data_aval <= d_col)
  resultado$aacpd_jan75[i] <- aacpd_trapezio(
    sev_jan$nota_media, sev_jan$data_aval)
}

cat("Computation completed.\n")

# =============================================================================
# 6. COMPARATIVE SUMMARY STATISTICS
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("Comparative Summary Statistics\n")
cat(strrep("=", 60), "\n\n")

for (v in c("aacpd_acum","aacpd_periodo","aacpd_jan75")) {
  vals <- resultado[[v]]
  cat(sprintf("%s:\n", v))
  cat(sprintf("  Valid N: %d | Missing (NAs): %d\n",
              sum(!is.na(vals)), sum(is.na(vals))))
  cat(sprintf("  Mean=%.1f | SD=%.1f | Min=%.1f | Max=%.1f\n\n",
              mean(vals, na.rm=TRUE), sd(vals, na.rm=TRUE),
              min(vals, na.rm=TRUE), max(vals, na.rm=TRUE)))
}

# Cross-version correlations
cat("Cross-version correlation matrix:\n")
cor_tab <- cor(resultado[, c("aacpd_acum","aacpd_periodo","aacpd_jan75")],
               use="pairwise.complete.obs")
print(round(cor_tab, 3))

# =============================================================================
# 7. MERGE BACK TO DATASET AND STANDARDIZE
# =============================================================================

cat("\nMerging indices into the primary dataset...\n")

resultado$individuo <- as.factor(resultado$individuo)
dataset$individuo   <- as.factor(as.integer(as.character(dataset$individuo)))

dataset_novo <- dataset |>
  left_join(
    resultado |>
      select(individuo, data_coleta,
             aacpd_acum, aacpd_periodo, aacpd_jan75),
    by = c("individuo", "data_coleta")
  )

# Standardize (Z-score transformation) the three indices
dataset_novo$aacpd_acum_sc    <- scale(dataset_novo$aacpd_acum)[,1]
dataset_novo$aacpd_periodo_sc <- scale(dataset_novo$aacpd_periodo)[,1]
dataset_novo$aacpd_jan75_sc   <- scale(dataset_novo$aacpd_jan75)[,1]

# Impute missing records with 0 (representing the population mean after scaling)
for (v in c("aacpd_acum_sc","aacpd_periodo_sc","aacpd_jan75_sc")) {
  dataset_novo[[v]][is.na(dataset_novo[[v]])] <- 0
}

cat(sprintf("  Updated dataset: %d observations\n", nrow(dataset_novo)))
cat(sprintf("  Missing entries (NAs) in aacpd_jan75_sc: %d\n",
            sum(is.na(dataset_novo$aacpd_jan75_sc))))

# =============================================================================
# 8. PRELIMINARY VALIDATION — CORRELATION WITH esp_casca_mm & rend_polpa_pct
# =============================================================================

cat("\nCorrelation between AUDPC indices and rind thickness (esp_casca_mm):\n")
for (v in c("aacpd_acum","aacpd_periodo","aacpd_jan75")) {
  r <- cor(dataset_novo[[v]], dataset_novo$esp_casca_mm,
           use="pairwise.complete.obs")
  cat(sprintf("  %s × esp_casca_mm: r=%.3f\n", v, r))
}

cat("\nCorrelation between AUDPC indices and pulp yield (rend_polpa_pct):\n")
for (v in c("aacpd_acum","aacpd_periodo","aacpd_jan75")) {
  r <- cor(dataset_novo[[v]], dataset_novo$rend_polpa_pct,
           use="pairwise.complete.obs")
  cat(sprintf("  %s × rend_polpa_pct: r=%.3f\n", v, r))
}

# =============================================================================
# 9. EXPORTING OUTPUTS
# =============================================================================

saveRDS(dataset_novo,
        here("01_dados/processados/dataset_aacpd_janela.rds"))

write.csv(resultado,
          here("01_dados/processados/aacpd_tres_versoes.csv"),
          row.names=FALSE)

# Export long-format severity data for future downstream tasks
saveRDS(severity_long,
        here("01_dados/processados/severity_long.rds"))

cat("\n✓ Process completed successfully.\n")
cat("  dataset_aacpd_janela.rds — Main dataset containing 3 distinct AUDPC tracks\n")
cat("  aacpd_tres_versoes.csv    — Comparative table indexed by individual harvest event\n")
cat("  severity_long.rds        — Long-format individual severity data tracking file\n")
cat("\nNext step: Contrast model performance for each version via DIC comparison.\n")