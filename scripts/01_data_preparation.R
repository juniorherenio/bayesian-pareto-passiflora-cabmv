# =============================================================================
# Adds spatial position (row × column) to the dataset via field layout map
# =============================================================================

source(here::here("02_scripts_R/00_setup.R"))

# -----------------------------------------------------------------------------
# 1. FIELD LAYOUT MAP LOADING
# -----------------------------------------------------------------------------
# The field layout map contains 240 genotypes arranged across 4 field rows (rows 11-14
# of the plot, renamed to 1-4). Columns vary by row:
#   Row 1: columns 10-76
#   Rows 2-3: columns 1-74
#   Row 4: columns 49-72
# F1/F2/F3 are family markers (borders/references) and are not genotypated.
# ' -' indicates an empty position.

croqui <- read.csv(here("01_dados/brutos/croqui_processado.csv")) |>
  mutate(
    individuo = as.character(genotipo),
    linha     = as.integer(linha),
    coluna    = as.integer(coluna)
  ) |>
  dplyr::select(individuo, linha, column = coluna) # Renamed to column for consistency

cat("Field layout map loaded:", nrow(croqui), "positions\n")
cat("Field rows:", sort(unique(croqui$linha)), "\n")
cat("Unique genotypes in layout map:", n_distinct(croqui$individuo), "\n\n")

# -----------------------------------------------------------------------------
# 2. MERGE WITH DATASET
# -----------------------------------------------------------------------------
# 'individuo' in the dataset = numeric plant ID (equivalent to 'genotipo' in the layout map)

dataset_espacial <- dataset |>
  mutate(individuo_chr = as.character(as.integer(as.character(individuo)))) |>
  left_join(
    croqui |> rename(individuo_chr = individuo),
    by = "individuo_chr"
  ) |>
  dplyr::select(-individuo_chr)

# Verify spatial coverage
n_com_pos  <- sum(!is.na(dataset_espacial$linha))
n_sem_pos  <- sum(is.na(dataset_espacial$linha))
ind_sem    <- dataset_espacial |>
  filter(is.na(linha)) |>
  distinct(individuo) |>
  pull(individuo)

cat("Observations with spatial coordinates:", n_com_pos, "\n")
cat("Observations without spatial coordinates (NA):", n_sem_pos, "\n")
if (length(ind_sem) > 0) {
  cat("Individuals missing from field layout map:", paste(ind_sem, collapse = ", "), "\n")
}

# -----------------------------------------------------------------------------
# 3. CONVERT row AND column TO ORDERED FACTORS (for asreml ar1())
# -----------------------------------------------------------------------------
dataset_espacial <- dataset_espacial |>
  mutate(
    linha  = factor(linha,  levels = sort(unique(linha[!is.na(linha)]))),
    coluna = factor(coluna, levels = sort(unique(coluna[!is.na(coluna)])))
  )

cat("\nRow levels: ", paste(levels(dataset_espacial$linha),  collapse = " "), "\n")
cat("Number of unique columns:", nlevels(dataset_espacial$coluna), "\n")

# -----------------------------------------------------------------------------
# 4. SAVE COMPLETED PROCESSED DATASET
# -----------------------------------------------------------------------------
saveRDS(dataset_espacial,
        file = here("01_dados/processados/dataset_completo.rds"))

write.csv(dataset_espacial,
          file = here("01_dados/processados/dataset_completo.csv"),
          row.names = FALSE)

# Update dataset object in global environment
dataset <<- dataset_espacial

cat("\n✓ dataset_completo saved to 01_dados/processados/\n")
cat("Final dataset:", nrow(dataset), "obs ×", ncol(dataset), "columns\n")
message("01_data_preparation.R completed.")