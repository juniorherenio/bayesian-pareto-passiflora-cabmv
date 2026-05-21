# =============================================================================
# GMM + UMAP Probabilistic Clustering over BC1 Breeding Populations
#
# MITIGATING DOUBLE SHRINKAGE BIAS:
#   Estimated Breeding Values (EBVs/BLUPs) are pre-shrunk by the animal model.
#   Methodological Solution: Weight each independent BLUP by sqrt(precision_kt) 
#   = 1/blup_sd_kt prior to downstream analysis. This functions as an arithmetic 
#   weighting adjustment, not a nested statistical estimation.
#   Genotypes characterized by high prediction uncertainty are pulled towards 
#   the background population mean, suppressing their weight in cluster formation.
#
# ANALYTICAL PIPELINE:
#   1. Construct an uncertainty-adjusted feature matrix weighted by prediction precision
#   2. Tune UMAP topology hyper-parameters (n_neighbors × min_dist) via silhouette + trustworthiness
#   3. Tune GMM components (n_components × covariance_type) via BIC + silhouette criteria
#   4. Validate empirical cluster partitions via multi-metric screening:
#      - BIC (Bayesian Information Criterion)
#      - Silhouette width score
#      - Calinski-Harabasz variance ratio index
#      - Davies-Bouldin separation index
#      - Shannon Entropy index (cluster assignment purity assessment)
#   5. Project 2D UMAP manifold visualization overlayed with optimal GMM clusters
# =============================================================================

library(here)
library(dplyr)
library(tidyverse)
library(tidyr)
library(mclust)    # GMM with automated BIC screening
library(uwot)      # UMAP manifold learning engine
library(cluster)   # Silhouette widths and cluster distance matrices
library(clusterSim) # Davies-Bouldin and Calinski-Harabasz indexing tools
library(ggplot2)

cat("\n", strrep("=", 60), "\n")
cat("07_clustering.R — GMM + UMAP Manifold Alignment\n")
cat(format(Sys.time(), "[%H:%M:%S]"), "\n")
cat(strrep("=", 60), "\n\n")

# =============================================================================
# 1. FEATURE MATRIX CONSTRUCTION & COMPRESSION
# =============================================================================

blups_est   <- readRDS(here(
  "04_resultados/blups/bayesiano/blups_GxA_POM_por_estacao.rds"))
componentes <- read.csv(here("04_resultados/selecao/componentes_indice.csv"))
pareto_bayes <- readRDS(here("04_resultados/selecao/pareto_bayesiano.rds"))
pom_ind     <- readRDS(here("01_dados/processados/pom_indices.rds"))
bdsip_std    <- readRDS(here("01_dados/processados/bdsip_completo.rds"))

cat("Compiling precision-weighted feature matrix...\n")

# Derive prediction error variances and precision arrays
blups_est <- blups_est |>
  mutate(
    var_pred = blup_sd^2,
    precisao = 1 / var_pred
  )

# =============================================================================
# ARITHMETIC PRECISION WEIGHTING — Bypasses Double Shrinkage Artifacts
# =============================================================================
# BLUP_pond = BLUP × sqrt(precision) = BLUP / blup_sd
# Biometric Rationale: Standardizes estimates by their specific prediction errors.
# Genotypes with high posterior uncertainty (inflated blup_sd) are compressed 
# toward zero. This avoids nested shrinkage inflation because it is an arithmetic 
# normalization step rather than a multi-stage random-effects model.
# =============================================================================

blups_pond <- blups_est |>
  mutate(
    blup_pond = blup * sqrt(precisao)   # Arithmetic precision-weighting transformation
  )

# Cast long array into wide matrix format: Genotype × (Trait_Season structural links)
feat_wide <- blups_pond |>
  mutate(feat_nome = paste0(trait, "_", estacao)) |>
  select(individuo, feat_nome, blup_pond) |>
  pivot_wider(names_from=feat_nome, values_from=blup_pond) |>
  column_to_rownames("individuo")

# Enforce matrix cleaning by dropping features containing missing inputs (NAs)
feat_wide <- feat_wide[, colSums(is.na(feat_wide)) == 0]

cat(sprintf("  Clean Feature Matrix Structure: %d genotypes × %d longitudinal traits mapped\n",
            nrow(feat_wide), ncol(feat_wide)))

# Standardize columns via Z-score scaling matrix
# Necessary prerequisite for UMAP and GMM algorithms — acts as scaling, not random-effects shrinkage
feat_sc <- scale(feat_wide)

# Merge back cross-metadata, validation indexes and target indices
resist_feat <- componentes |>
  mutate(individuo = as.character(individuo)) |>
  select(individuo, IR_comercial, IR_res_sc,
         IR_est_sc, IR_final, var_pred_medio) |>
  filter(!is.na(IR_res_sc))

# Identify intersection cohort containing completed profiles
ind_completos <- intersect(rownames(feat_sc),
                           resist_feat$individuo)
feat_sc_sub <- feat_sc[ind_completos, ]

cat(sprintf("Genotypes with full feature matrix coverage: %d lines\n\n",
            length(ind_completos)))

# =============================================================================
# 2. UMAP TOPOLOGICAL MANIFOLD HYPER-PARAMETER TUNING
# =============================================================================

cat(strrep("=", 60), "\n")
cat("UMAP Manifold Parameter Optimization\n")
cat(strrep("=", 60), "\n\n")

# Hyper-parameter search space optimization grid
n_neighbors_grid <- c(5, 10, 15, 20)
min_dist_grid    <- c(0.01, 0.05, 0.1, 0.3)
n_components_umap <- 2

# Primary validation metric: Trustworthiness index
# Quantifies how effectively local structure neighborhoods are preserved during projection
# from high-dimensional matrix spaces down to low-dimensional embedding coordinates.

trustworthiness <- function(X_orig, X_emb, k=10) {
  n <- nrow(X_orig)
  d_orig <- as.matrix(dist(X_orig))
  d_emb  <- as.matrix(dist(X_emb))
  
  rank_orig <- apply(d_orig, 1, rank)
  rank_emb  <- apply(d_emb,  1, rank)
  
  t_sum <- 0
  for (i in 1:n) {
    # Isolate k-nearest neighbors in low-dimensional space coordinates
    nn_emb <- order(d_emb[i,])[-1][1:k]
    for (j in nn_emb) {
      r_ij <- rank_orig[j, i] - 1
      if (r_ij > k) t_sum <- t_sum + (r_ij - k)
    }
  }
  1 - (2 / (n * k * (2*n - 3*k - 1))) * t_sum
}

cat("UMAP Hyper-parameter Tuning Matrix (n_neighbors × min_dist):\n")
cat(sprintf("  %d × %d = %d search space configurations explored\n\n",
            length(n_neighbors_grid), length(min_dist_grid),
            length(n_neighbors_grid) * length(min_dist_grid)))

umap_results <- list()
umap_grid    <- expand.grid(n_neighbors=n_neighbors_grid,
                            min_dist=min_dist_grid)

for (i in 1:nrow(umap_grid)) {
  nn <- umap_grid$n_neighbors[i]
  md <- umap_grid$min_dist[i]
  
  set.seed(42)
  emb <- tryCatch(
    umap(feat_sc_sub, n_neighbors=nn, min_dist=md,
         n_components=n_components_umap, verbose=FALSE),
    error = function(e) NULL
  )
  if (is.null(emb)) next
  
  # Calculate topological neighborhood preservation
  tw <- trustworthiness(feat_sc_sub, emb, k=min(nn, 10))
  
  # Evaluate cluster separation over low-dimensional space via an intermediate GMM model matrix (k=3 baseline)
  set.seed(42)
  gmm_tmp <- Mclust(emb, G=3, verbose=FALSE)
  sil_tmp <- if (!is.null(gmm_tmp)) {
    mean(silhouette(gmm_tmp$classification, dist(emb))[,3])
  } else NA
  
  umap_results[[i]] <- list(
    n_neighbors = nn,
    min_dist    = md,
    trust       = round(tw, 4),
    silhouette  = round(sil_tmp, 4),
    embedding   = emb
  )
  
  cat(sprintf("  nn=%2d | md=%.2f | trustworthiness=%.4f | silhouette width=%.4f\n",
              nn, md, tw, sil_tmp))
}

# Select optimal manifold embedding coordinates
# Fitness function maximizing: trustworthiness × silhouette width (structural balance)
scores_umap <- sapply(umap_results, function(r) {
  if (is.null(r)) return(-Inf)
  tw  <- ifelse(is.na(r$trust), 0, r$trust)
  sil <- ifelse(is.na(r$silhouette), 0, r$silhouette)
  tw * sil
})

best_umap_idx <- which.max(scores_umap)
best_umap     <- umap_results[[best_umap_idx]]

cat(sprintf("\n✓ Optimal Manifold Model Reached: n_neighbors=%d | min_dist=%.2f\n",
            best_umap$n_neighbors, best_umap$min_dist))
cat(sprintf("  trustworthiness parameter=%.4f | average silhouette width=%.4f\n\n",
            best_umap$trust, best_umap$silhouette))

emb_best <- best_umap$embedding

# =============================================================================
# 3. GMM COVARIANCE STRUCTURE TYPING & CLUSTER COMPONENT SELECTION
# =============================================================================

cat(strrep("=", 60), "\n")
cat("Gaussian Mixture Model (GMM) Component Optimization\n")
cat(strrep("=", 60), "\n\n")

# mclust executes automated grid screening across multiple geometric covariance matrix types
# Geometric parameterizations evaluated via Bayesian Information Criterion (BIC):
#   Spherical / Diagonal: EII, VII, EEI, VEI, EVI, VVI
#   Ellipsoidal / Unstructured: EEE, EVE, VEE, VVE, EEV, VEV, EVV, VVV

cat("Fitting GMM parameters with automated structural selection via BIC optimization...\n")
cat("  Screening components G=2..8 | Covariance configurations: EII,VII,EEI,VEI,EVI,VVI,EEE,VVV\n\n")

set.seed(42)
gmm_bic <- Mclust(emb_best,
                  G = 2:8,
                  modelNames = c("EII","VII","EEI","VEI",
                                 "EVI","VVI","EEE","VVV"),
                  verbose = FALSE)

cat("BIC Diagnostic Selection Matrix Summary across components:\n")
print(summary(gmm_bic$BIC))

cat(sprintf("\n✓ Optimal Covariance Model: %s parameterization | Optimal K=%d components | BIC score=%.2f\n\n",
            gmm_bic$modelName, gmm_bic$G, gmm_bic$BIC))

# Compile comparative diagnostics block for fixed cluster arrays (G=2..6 manually verified)
cat("Comparative Metrics Matrix for Fixed Partition Ranges G=2..6:\n")
cat(sprintf("  %-8s %-10s %-10s %-10s %-10s\n",
            "K-groups", "BIC", "Silhouette", "CH-Index", "DB-Index"))

gmm_manual <- list()
for (G in 2:6) {
  set.seed(42)
  m <- tryCatch(
    Mclust(emb_best, G=G,
           modelNames=gmm_bic$modelName,
           verbose=FALSE),
    error=function(e) NULL
  )
  if (is.null(m)) next
  
  cls   <- m$classification
  d_mat <- dist(emb_best)
  
  # Mean silhouette width
  sil <- mean(silhouette(cls, d_mat)[,3])
  
  # Calinski-Harabasz variance ratio index
  ch  <- tryCatch(
    index.G1(emb_best, cls),
    error=function(e) NA
  )
  
  # Davies-Bouldin cluster separation matrix parameter
  db  <- tryCatch(
    index.DB(emb_best, cls)$DB,
    error=function(e) NA
  )
  
  gmm_manual[[G]] <- list(
    G=G, BIC=m$BIC, sil=sil, ch=ch, db=db, model=m
  )
  
  cat(sprintf("  G=%-6d %-10.2f %-10.4f %-10.4f %-10.4f\n",
              G, m$BIC, sil, ch, db))
}

# =============================================================================
# 4. FINAL INTEGRATED SELECTION SCORE FOR OPTIMAL COMPONENT K
# =============================================================================

cat("\n", strrep("=", 60), "\n")
cat("Consensus Identification of the Optimal Number of Clusters\n")
cat(strrep("=", 60), "\n\n")

# Mathematical integration of criteria bounds:
#   BIC (Bayesian Information Criterion): Maximize score
#   Silhouette average width: Maximize score within [0,1]
#   Calinski-Harabasz (CH-Index): Maximize variance ratio
#   Davies-Bouldin (DB-Index): Minimize distance ratio

metricas <- do.call(rbind, lapply(gmm_manual, function(r) {
  if (is.null(r)) return(NULL)
  data.frame(G=r$G, BIC=r$BIC, sil=r$sil,
             ch=ifelse(is.na(r$ch),0,r$ch),
             db=ifelse(is.na(r$db),Inf,r$db))
}))

if (nrow(metricas) > 0) {
  metricas <- metricas |>
    mutate(
      BIC_norm = (BIC - min(BIC)) / (max(BIC) - min(BIC)),
      sil_norm = (sil - min(sil)) / (max(sil) - min(sil)),
      ch_norm  = (ch  - min(ch))  / (max(ch)  - min(ch)),
      db_norm  = 1 - (db - min(db)) / (max(db) - min(db)),
      score    = (BIC_norm + sil_norm + ch_norm + db_norm) / 4
    ) |>
    arrange(desc(score))
  
  cat("Consensus Integrated Selection Score Matrix (Normalized 0-1, higher = optimal):\n")
  print(metricas[, c("G","BIC","sil","ch","db","score")],
        row.names=FALSE)
  
  G_otimo <- metricas$G[1]
  cat(sprintf("\n✓ Optimal partition defined by Integrated Selection Score: K=%d clusters\n\n", G_otimo))
} else {
  G_otimo <- gmm_bic$G
  cat(sprintf("  Fallback mechanism triggered: Deploying automated BIC optimal K=%d components\n\n", G_otimo))
}

# =============================================================================
# 5. FINAL GAUSSIAN MIXTURE MODEL PROBABILISTIC EXECUTION
# =============================================================================

set.seed(42)
gmm_final <- Mclust(emb_best, G=G_otimo,
                    modelNames=gmm_bic$modelName,
                    verbose=FALSE)

cat(sprintf("Final Model Specification: %s covariance structure | K=%d components mapped\n",
            gmm_final$modelName, gmm_final$G))
cat(sprintf("BIC checkpoint=%.2f | Maximized logLikelihood=%.2f\n\n",
            gmm_final$BIC, gmm_final$loglik))

# Extract soft-clustering posterior probability vectors
probs_cluster <- gmm_final$z
colnames(probs_cluster) <- paste0("P_cluster_", 1:G_otimo)
cluster_assign <- gmm_final$classification

cat("Genotype assignment distribution per cluster matrix:\n")
print(table(cluster_assign))

# Compute posterior assignment Shannon Entropy profiles per individual genotype
entropia <- apply(probs_cluster, 1, function(p) {
  p <- p[p > 0]; -sum(p * log(p))
})
cat(sprintf("\nPopulation Mean Entropy: %.3f (Theoretical absolute uncertainty max: %.3f)\n\n",
            mean(entropia), log(G_otimo)))

# Validation summary matrix
d_mat     <- dist(emb_best)
sil_final <- silhouette(cluster_assign, d_mat)
sil_medio <- mean(sil_final[,3])
ch_final  <- tryCatch(index.G1(emb_best, cluster_assign), error=function(e) NA)
db_final  <- tryCatch(index.DB(emb_best, cluster_assign)$DB, error=function(e) NA)

cat(sprintf("Final Cluster Validity Metrics:\n"))
cat(sprintf("  Average Silhouette Width: %.4f\n", sil_medio))
cat(sprintf("  Calinski-Harabasz Index:  %.2f\n", ch_final))
cat(sprintf("  Davies-Bouldin Index:     %.4f\n\n", db_final))

cat("Internal Silhouette widths breakdown per component block:\n")
for (k in 1:G_otimo) {
  sil_k <- sil_final[cluster_assign == k, 3]
  cat(sprintf("  Component K=%d: average width=%.3f | allocation size=%d lines | lower bound=%.3f\n",
              k, mean(sil_k), length(sil_k), min(sil_k)))
}

# =============================================================================
# 6. STRUCTURAL BIOMETRIC & ECO-EPIDEMIOLOGICAL INTERPRETATION
# =============================================================================

cluster_df <- data.frame(
  individuo = ind_completos,
  cluster   = cluster_assign,
  entropia  = entropia,
  umap1     = emb_best[,1],
  umap2     = emb_best[,2]
) |>
  left_join(
    componentes |>
      mutate(individuo=as.character(individuo)) |>
      select(individuo, familia, cruzamento,
             IR_comercial, IR_res_sc, IR_final,
             P_severo_media, var_pred_medio),
    by="individuo"
  ) |>
  left_join(
    pareto_bayes |> select(individuo, P_pareto),
    by="individuo"
  )

cat("\nSub-component breeding parameter means grouped by cluster partition:\n")
resumo_cluster <- cluster_df |>
  group_by(cluster) |>
  summarise(
    n              = n(),
    IR_com_medio   = round(mean(IR_comercial, na.rm=TRUE), 3),
    IR_res_medio   = round(mean(IR_res_sc, na.rm=TRUE), 3),
    IR_final_medio = round(mean(IR_final, na.rm=TRUE), 3),
    P_severo_medio = round(mean(P_severo_media, na.rm=TRUE), 3),
    P_pareto_medio = round(mean(P_pareto, na.rm=TRUE), 3),
    entropia_medio = round(mean(entropia), 3),
    .groups="drop"
  )
print(resumo_cluster, row.names=FALSE)

cat("\nCross-tabulation array: Full-Sib Family Structures × Cluster allocation:\n")
print(table(cluster_df$familia, cluster_df$cluster))

cat("\nBiometric Ideotype Assignment Characterization:\n")
for (k in 1:G_otimo) {
  r <- resumo_cluster |> filter(cluster == k)
  perfil <- case_when(
    r$IR_com_medio > 0.1  & r$IR_res_medio > 0.1  ~ "Balanced Ideotype (Recovered commercial traits + Introgressed resistance)",
    r$IR_com_medio > 0.1  & r$IR_res_medio <= 0.1 ~ "P. edulis Phenotypic Type (High commercial merit / Susceptible)",
    r$IR_com_medio <= 0.1 & r$IR_res_medio > 0.1  ~ "P. setacea Phenotypic Type (Low commercial merit / Highly resistant)",
    TRUE                                          ~ "Compromised / Transgressive Intermediate Variant"
  )
  cat(sprintf("  Cluster K=%d (Size count=%d lines): Ideotype Profile Class: %s\n", k, r$n, perfil))
  cat(sprintf("     Agronomic Index=%.3f | Resistance Index=%.3f | Chronic Pathosystem Probability=%.3f\n",
              r$IR_com_medio, r$IR_res_medio, r$P_severo_media))
}

# Identify allocation clusters for Consensus Elite Parents
consenso_ids <- c("97","9","131","2","83","223")
cat("\nConsensus Recommended Elite Parents (NSGA-II Frontiers) Mapping across Clusters:\n")
cluster_df |>
  filter(individuo %in% consenso_ids) |>
  select(individuo, cluster, IR_comercial, IR_res_sc, P_pareto) |>
  arrange(cluster, desc(P_pareto)) |>
  print(row.names=FALSE)

# =============================================================================
# 7. OUTPUT MATRIX GENERATION & DATA EXPORTING
# =============================================================================

dir.create(here("04_resultados/selecao"), showWarnings=FALSE, recursive=TRUE)
saveRDS(list(emb=emb_best, gmm=gmm_final, clusters=cluster_df,
             metricas=metricas,
             umap_params=list(n_neighbors=best_umap$n_neighbors,
                              min_dist=best_umap$min_dist)),
        here("04_resultados/selecao/clustering_resultado.rds"))
write.csv(cluster_df,
          here("04_resultados/selecao/clusters_genitores.csv"),
          row.names=FALSE)
write.csv(resumo_cluster,
          here("04_resultados/selecao/resumo_clusters.csv"),
          row.names=FALSE)

cat("\n✓ Module 07 mathematical execution successfully completed.\n")
cat("  Exported Cluster Array:    clusters_genitores.csv  — Allocation layout matrix per genotype\n")
cat("  Exported Parameter Means:   resumo_clusters.csv     — Component-specific biometric summary\n")
cat("  Exported Checkpoint RDS:   clustering_resultado.rds — Integrated structural model state\n")
cat("\nNext downstream pipeline execution module path: 08_dashboard.R (AI + Shiny deployment interface)\n")