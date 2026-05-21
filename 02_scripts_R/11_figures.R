# =============================================================================
# Master plotting script for multi-trait selection and genetic parameters
# =============================================================================

library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)
library(here)

FONT_BASE   <- "serif"
PAL_CLUSTER <- c("1"="#1B4F72","2"="#C0392B","3"="#1E8449")
PAL_FAM     <- c("1"="#1B4F72","2"="#D4AC0D","3"="#C0392B")

theme_q1 <- function(base_size=11) {
  theme_classic(base_family=FONT_BASE, base_size=base_size) +
    theme(
      axis.line        = element_line(color="grey25", linewidth=0.45),
      axis.ticks       = element_line(color="grey25", linewidth=0.3),
      axis.text        = element_text(color="grey20", size=base_size-2),
      axis.title       = element_text(color="grey10", size=base_size-1,
                                      face="bold"),
      panel.grid.major = element_line(color="grey93", linewidth=0.25),
      panel.grid.minor = element_blank(),
      legend.title     = element_text(size=base_size-2, face="bold"),
      legend.text      = element_text(size=base_size-3),
      legend.key.size  = unit(0.38,"cm"),
      legend.background= element_rect(fill=alpha("white",0.85), color=NA),
      plot.margin      = margin(10,10,10,10),
      plot.background  = element_rect(fill="white", color=NA),
      panel.background = element_rect(fill="white", color=NA)
    )
}

clusters <- read.csv(here("04_resultados/selecao/clusters_genitores.csv")) |>
  mutate(
    cluster_f = factor(cluster, levels=c(1,2,3)),
    familia_f = factor(familia, levels=c(1,2,3)),
    consenso  = as.character(individuo) %in% c("97","9","131","2","83","223"),
    individuo = as.character(individuo)
  )

pareto <- read.csv(here("04_resultados/selecao/pareto_bayesiano.csv")) |>
  mutate(
    individuo = as.character(individuo),
    familia_f = factor(familia, levels=c(1,2,3)),
    consenso  = individuo %in% c("97","9","131","2","83","223")
  )

rg <- read.csv(here(
  "04_resultados/parametros_geneticos/bayesiano/rg_bivariado.csv"))

dir.create(here("06_manuscrito/figuras"), showWarnings=FALSE, recursive=TRUE)


# =============================================================================
# FIGURE 1 — PROBABILISTIC MANIFOLD CLUSTERING (UMAP + GMM)
# =============================================================================

FONT_BASE   <- "serif"
PAL_CLUSTER <- c("1"="#1B4F72","2"="#C0392B","3"="#1E8449")
PAL_SHAPE   <- c("1"=16,"2"=17,"3"=15)

# Biologically validated cluster definitions
CLUSTER_LABELS <- c(
  "1" = "Cluster I — Balanced ideotype (293 × 355)",
  "2" = "Cluster II — Poor allelic combination (516 × 355)",
  "3" = "Cluster III — Low commercial recovery (501 × 355)"
)

theme_q1 <- function(base_size=10) {
  theme_classic(base_family=FONT_BASE, base_size=base_size) +
    theme(
      axis.line        = element_line(color="grey25", linewidth=0.4),
      axis.ticks       = element_line(color="grey25", linewidth=0.3),
      axis.text        = element_text(color="grey25", size=base_size-2),
      axis.title       = element_text(color="grey10", size=base_size-1,
                                      face="bold"),
      panel.grid.major = element_line(color="grey93", linewidth=0.22),
      panel.grid.minor = element_blank(),
      legend.title     = element_text(size=base_size-1, face="bold"),
      legend.text      = element_text(size=base_size-2),
      legend.key.size  = unit(0.35,"cm"),
      plot.background  = element_rect(fill="white", color=NA),
      panel.background = element_rect(fill="white", color=NA),
      plot.margin      = margin(5,5,5,5),
      plot.tag         = element_text(face="bold", size=base_size+2,
                                      family=FONT_BASE)
    )
}

clusters <- read.csv(here("04_resultados/selecao/clusters_genitores.csv")) |>
  mutate(
    cluster_f = factor(cluster, levels=c(1,2,3)),
    familia_f = factor(familia, levels=c(1,2,3)),
    consenso  = as.character(individuo) %in%
      c("97","9","131","2","83","223"),
    individuo = as.character(individuo)
  )

buf <- function(x, pad=0.4) range(x) + c(-pad, pad)

cent_all <- clusters |>
  group_by(cluster_f) |>
  summarise(cx=mean(umap1), cy=mean(umap2), .groups="drop") |>
  mutate(label=c("I","II","III"))

# ── Panel A — Global Manifold Projection ─────────────────────────────────────
pA <- ggplot(clusters, aes(umap1, umap2)) +
  stat_ellipse(aes(fill=cluster_f, color=cluster_f),
               geom="polygon", level=0.95, alpha=0.10,
               linewidth=NA) +
  stat_ellipse(aes(color=cluster_f), level=0.90, linewidth=0.5) +
  geom_point(aes(color=cluster_f, shape=familia_f),
             size=1.5, alpha=0.65, stroke=0.2) +
  geom_point(data=filter(clusters, consenso),
             shape=21, size=3.0, fill="#F4D03F",
             color="#7D6608", stroke=1.1) +
  geom_text(data=cent_all,
            aes(cx, cy+1.5, label=label, color=cluster_f),
            size=4.5, fontface="bold", family=FONT_BASE,
            show.legend=FALSE) +
  scale_color_manual(
    values  = PAL_CLUSTER,
    labels  = CLUSTER_LABELS,
    name    = "Cluster"
  ) +
  scale_fill_manual(values=PAL_CLUSTER, guide="none") +
  scale_shape_manual(
    values = PAL_SHAPE,
    labels = c("1"="FBC₁ (293 × 355)",
               "2"="FBC₁ (501 × 355)",
               "3"="FBC₁ (516 × 355)"),
    name   = "Full-sib family"
  ) +
  guides(
    color = guide_legend(order=1, nrow=3, 
                         override.aes=list(size=3, shape=16, linetype=0)),
    shape = guide_legend(order=2, nrow=3,
                         override.aes=list(size=3, color="grey35"))
  ) +
  labs(x="UMAP 1", y="UMAP 2", tag="A") +
  theme_q1()

# ── Panel B — Magnified Target Cohort (Cluster I) ──────────────────────────
cl1 <- filter(clusters, cluster==1)

pB <- ggplot(cl1, aes(umap1, umap2)) +
  stat_ellipse(aes(fill=cluster_f, color=cluster_f),
               geom="polygon", level=0.95, alpha=0.10,
               linewidth=NA) +
  stat_ellipse(aes(color=cluster_f), level=0.90, linewidth=0.5) +
  geom_point(data=filter(cl1, !consenso),
             color=PAL_CLUSTER["1"], shape=16,
             size=2.2, alpha=0.60, stroke=0.2) +
  geom_point(data=filter(cl1, consenso),
             shape=21, size=5.0, fill="#F4D03F",
             color="#7D6608", stroke=1.4) +
  geom_label_repel(
    data=filter(cl1, consenso),
    aes(label=individuo),
    size=3.2, fontface="bold", family=FONT_BASE,
    color="grey10",
    label.padding=unit(0.13,"lines"), label.size=0.22,
    box.padding=0.8, point.padding=0.5,
    force=10, min.segment.length=0.2,
    segment.color="grey40", segment.size=0.4,
    fill=alpha("white",0.93),
    show.legend=FALSE, max.overlaps=Inf, seed=42
  ) +
  scale_color_manual(values=PAL_CLUSTER, guide="none") +
  scale_fill_manual(values=PAL_CLUSTER, guide="none") +
  coord_cartesian(xlim=buf(cl1$umap1, 0.5),
                  ylim=buf(cl1$umap2, 0.5)) +
  labs(x="UMAP 1", y="UMAP 2", tag="B") +
  theme_q1() +
  theme(legend.position="none")

# ── Panel C — Magnified Outlier Cohort (Cluster II) ─────────────────────────
cl2 <- filter(clusters, cluster==2)

pC <- ggplot(cl2, aes(umap1, umap2)) +
  stat_ellipse(aes(fill=cluster_f, color=cluster_f),
               geom="polygon", level=0.95, alpha=0.12,
               linewidth=NA) +
  stat_ellipse(aes(color=cluster_f), level=0.90, linewidth=0.5) +
  geom_point(color=PAL_CLUSTER["2"], shape=15,
             size=2.5, alpha=0.70, stroke=0.2) +
  scale_color_manual(values=PAL_CLUSTER, guide="none") +
  scale_fill_manual(values=PAL_CLUSTER, guide="none") +
  coord_cartesian(xlim=buf(cl2$umap1, 0.3),
                  ylim=buf(cl2$umap2, 0.3)) +
  labs(x="UMAP 1", y="UMAP 2", tag="C") +
  theme_q1() +
  theme(legend.position="none")

# ── Panel D — Magnified Outlier Cohort (Cluster III) ────────────────────────
cl3 <- filter(clusters, cluster==3)

pD <- ggplot(cl3, aes(umap1, umap2)) +
  stat_ellipse(aes(fill=cluster_f, color=cluster_f),
               geom="polygon", level=0.95, alpha=0.12,
               linewidth=NA) +
  stat_ellipse(aes(color=cluster_f), level=0.90, linewidth=0.5) +
  geom_point(color=PAL_CLUSTER["3"], shape=17,
             size=2.5, alpha=0.70, stroke=0.2) +
  scale_color_manual(values=PAL_CLUSTER, guide="none") +
  scale_fill_manual(values=PAL_CLUSTER, guide="none") +
  coord_cartesian(xlim=buf(cl3$umap1, 0.3),
                  ylim=buf(cl3$umap2, 0.3)) +
  labs(x="UMAP 1", y="UMAP 2", tag="D") +
  theme_q1() +
  theme(legend.position="none")

# ── Grid Compilation and Layout Optimization ─────────────────────────────────
fig1 <- (pA | (pB / (pC | pD))) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",          
    legend.justification = "center",    
    legend.margin = margin(t = 10, b = 15) 
  )

ggsave(here("06_manuscrito/figuras/fig1_umap_final.pdf"),
       fig1, width=24, height=14.5, units="cm", dpi=300) 
ggsave(here("06_manuscrito/figuras/fig1_umap_final.png"),
       fig1, width=24, height=14.5, units="cm", dpi=300, bg="white")
cat("✓ Figure 1 exported successfully with integrated layouts and aligned labels.\n")

# =============================================================================
# FIGURE 2 — LONGITUDINAL HERITABILITY (h²_g Heatmap matrix)
# =============================================================================

h2_plot <- h2_est |>
  mutate(
    trait_f   = factor(trait, levels=rev(names(trait_labels)),
                       labels=rev(trait_labels)),
    estacao_f = factor(estacao,
                       levels=c("i18","p18","v19","o19","i19"),
                       labels=estacao_labels)
  )

fig2 <- ggplot(h2_plot, aes(estacao_f, trait_f, fill=h2_media)) +
  
  geom_tile(color="white", linewidth=0.6) +
  
  # Structural values annotations
  geom_text(aes(label=sprintf("%.2f", h2_media)),
            size=2.8, family=FONT_BASE, fontface="bold",
            color=ifelse(h2_plot$h2_media > 0.38, "white", "grey15")) +
  
  # 95% Bayesian HPD intervals mapping tracking
  geom_errorbar(
    aes(xmin=as.numeric(estacao_f) - 0.45,
        xmax=as.numeric(estacao_f) + 0.45,
        ymin=as.numeric(trait_f) + 0.28,
        ymax=as.numeric(trait_f) + 0.28),
    linewidth=0
  ) +
  
  scale_fill_gradientn(
    colors=c("#EBF5FB","#2471A3","#1B2631"),
    values=c(0, 0.5, 1),
    limits=c(0.2, 0.56),
    name=expression(italic(h)^2),
    guide=guide_colorbar(barwidth=0.5, barheight=5,
                         ticks.linewidth=0.3)
  ) +
  
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(x=NULL, y=NULL) +
  theme_q1() +
  theme(
    axis.text.x    = element_text(angle=30, hjust=1, size=8.5),
    axis.text.y    = element_text(size=9),
    panel.grid     = element_blank(),
    legend.position= "right"
  )

ggsave(here("06_manuscrito/figuras/fig2_h2_heatmap.pdf"),
       fig2, width=14, height=10, units="cm", dpi=300)
ggsave(here("06_manuscrito/figuras/fig2_h2_heatmap.png"),
       fig2, width=14, height=10, units="cm", dpi=300)
cat("✓ Figure 2 longitudinal heritability matrix heatmap generated successfully.\n")

# =============================================================================
# FIGURE 3 — UNCERTAINTY-AWARE POSTERIOR PARETO OPTIMIZATION FRONTIER
# =============================================================================

pareto_norm <- filter(pareto, !consenso)
pareto_dest <- filter(pareto, consenso)

fig3 <- ggplot() +
  
  # Multi-objective selection quadrant partitioning overlays
  annotate("rect", xmin=-Inf,xmax=0, ymin=0,ymax=Inf,
           fill="#D6EAF8", alpha=0.40) +
  annotate("rect", xmin=0,xmax=Inf, ymin=0,ymax=Inf,
           fill="#D5F5E3", alpha=0.40) +
  annotate("rect", xmin=-Inf,xmax=0, ymin=-Inf,ymax=0,
           fill="#FADBD8", alpha=0.40) +
  annotate("rect", xmin=0,xmax=Inf, ymin=-Inf,ymax=0,
           fill="#FEF9E7", alpha=0.40) +
  
  # Ideotype boundary texts placement anchoring
  annotate("text", x=-Inf, y=-Inf,
           label="Low resistance\nLow commercial",
           size=3, color="grey50", family=FONT_BASE,
           hjust=-0.1, vjust=-0.5, lineheight=0.9) +
  
  annotate("text", x=Inf, y=-Inf,
           label="Low resistance\nHigh commercial",
           size=3, color="grey50", family=FONT_BASE,
           hjust=1.1, vjust=-0.5, lineheight=0.9) +
  
  annotate("text", x=-Inf, y=Inf,
           label="High resistance\nLow commercial",
           size=3, color="grey50", family=FONT_BASE,
           hjust=-0.1, vjust=1.5, lineheight=0.9) +
  
  annotate("text", x=Inf, y=Inf,
           label="High resistance\nHigh commercial\nIdeal region",
           size=3.2, color="#1A5E20", fontface="bold",
           family=FONT_BASE, hjust=1.1, vjust=1.2, lineheight=0.9) +
  
  geom_hline(yintercept=0, color="grey45", linewidth=0.35, linetype="dashed") +
  geom_vline(xintercept=0, color="grey45", linewidth=0.35, linetype="dashed") +
  
  # Target tracking data distributions array
  geom_point(data=pareto_norm,
             aes(IR_comercial, IR_res_sc,
                 size=P_pareto, color=familia_f),
             alpha=0.60, shape=16) +
  
  # Consensus Elite Targets — Fixed layout diamond tracking points
  geom_point(data=pareto_dest,
             aes(IR_comercial, IR_res_sc, fill=familia_f),
             shape=23, size=5.5, color="black", 
             stroke=0.8, alpha=1) +
  
  # Labels overlay mappings using standard halo mechanisms
  geom_text_repel(
    data=pareto_dest,
    aes(IR_comercial, IR_res_sc, label=individuo),
    size=3.5, fontface="bold", family=FONT_BASE,
    color="black", bg.color="white", bg.r=0.15,
    box.padding=0.8, point.padding=0.4,
    force=8, min.segment.length=0,
    segment.color="grey30", segment.size=0.4,
    show.legend=FALSE, max.overlaps=Inf, seed=42
  ) +
  
  scale_color_manual(values=PAL_FAM,
                     labels=c("1"="293 × 355",
                              "2"="501 × 355",
                              "3"="516 × 355"),
                     name="Full-sib cross") +
  scale_fill_manual(values=PAL_FAM, guide="none") +
  scale_size_continuous(range=c(1.2, 5.8),
                        name=expression(italic(P)[Pareto]),
                        breaks=c(0.05,0.15,0.25,0.30)) +
  
  guides(
    color=guide_legend(order=1,
                       override.aes=list(size=3.5, alpha=1)),
    size =guide_legend(order=2,
                       override.aes=list(color="grey40", shape=16))
  ) +
  
  scale_x_continuous(expand=expansion(mult=c(0.05, 0.22))) +
  scale_y_continuous(expand=expansion(mult=c(0.05, 0.15))) +
  
  labs(x=expression(IR[commercial]),
       y=expression(IR[resistance])) +
  theme_q1() +
  theme(legend.position="right")

ggsave(here("06_manuscrito/figuras/fig3_pareto_v3.pdf"),
       fig3, width=17, height=14, units="cm", dpi=300)
ggsave(here("06_manuscrito/figuras/fig3_pareto_v3.png"),
       fig3, width=17, height=14, units="cm", dpi=300, bg="white")
cat("✓ Figure 3 Pareto optimization frontier chart updated successfully.\n")

# =============================================================================
# FIGURE 4 — GENETIC CORRELATIONS NETWORK MATRIX
# =============================================================================

FONT_BASE <- "serif"

# Node structures configuration definitions
trait_labels <- c(
  comp_mm        = "Length\n(mm)",
  diam_mm        = "Diameter\n(mm)",
  massa_fruto_g  = "Fruit mass\n(g)",
  rend_polpa_pct = "Pulp yield\n(%)",
  esp_casca_mm   = "Rind\nthickness (mm)",
  sst_brix       = "TSS\n(°Brix)",
  massa_polpa_g  = "Pulp mass\n(g)"
)

nodes <- data.frame(
  trait = names(trait_labels),
  label = unname(trait_labels),
  group = c("Morphology","Morphology","Morphology",
            "Quality","Disease","Quality","Quality"),
  stringsAsFactors=FALSE
)

# Circular layout mapping generation matrix
n  <- nrow(nodes)
th <- seq(pi/2, pi/2 - 2*pi*(n-1)/n, length.out=n)
nodes$x <- cos(th)
nodes$y <- sin(th)

# Local coordinate boundary text adjustments
r_lab <- 1.30
nodes$lx    <- cos(th) * r_lab
nodes$ly    <- sin(th) * r_lab
nodes$hjust <- ifelse(cos(th) >  0.15, 0,
                      ifelse(cos(th) < -0.15, 1, 0.5))
nodes$vjust <- ifelse(sin(th) >  0.15, 0,
                      ifelse(sin(th) < -0.15, 1, 0.5))

edges <- rg |>
  left_join(nodes |> select(trait,x,y) |>
              rename(trait1=trait,x1=x,y1=y), by="trait1") |>
  left_join(nodes |> select(trait,x,y) |>
              rename(trait2=trait,x2=x,y2=y), by="trait2") |>
  mutate(
    sinal  = ifelse(rg >= 0, "Positive","Negative"),
    abs_rg = abs(rg),
    xm     = (x1+x2)/2 * 0.75,   
    ym     = (y1+y2)/2 * 0.75
  )

pal_group <- c("Morphology"="#2471A3",
               "Quality"    ="#6C3483",
               "Disease"    ="#C0392B")
pal_rg    <- c("Positive"="#1E8449","Negative"="#C0392B")

fig4 <- ggplot() +
  
  # Background structural matrices links (negligible lines)
  geom_segment(
    data=filter(edges, abs_rg <= 0.10),
    aes(x=x1,y=y1,xend=x2,yend=y2),
    color="grey85", linewidth=0.3, alpha=0.6
  ) +
  
  # Robust link parameters (significant phenotypic lines mapped by sign)
  geom_segment(
    data=filter(edges, abs_rg > 0.10),
    aes(x=x1,y=y1,xend=x2,yend=y2,
        color=sinal, linewidth=abs_rg),
    alpha=0.85, lineend="round"
  ) +
  
  # Marginal parameter estimates overlay mappings
  geom_label(
    data=filter(edges, abs_rg > 0.10),
    aes(x=xm, y=ym,
        label=sprintf("%+.2f", rg),
        color=sinal),
    size=3.0, family=FONT_BASE, fontface="bold",
    label.padding=unit(0.15,"lines"),
    label.size=NA, 
    fill=alpha("white", 0.85),
    show.legend=FALSE
  ) +
  
  # Core nodes plots mapping
  geom_point(
    data=nodes,
    aes(x, y, fill=group),
    shape=21, size=13, color="white", stroke=2.0
  ) +
  
  # Core text arrays labeling
  geom_text(
    data=nodes,
    aes(lx, ly, label=label,
        hjust=hjust, vjust=vjust),
    size=3.2, family=FONT_BASE, fontface="bold",
    color="grey15", lineheight=0.85
  ) +
  
  scale_color_manual(
    values=pal_rg,
    name="Genetic\ncorrelation",
    labels=c("Negative"="Negative  (rg < 0)",
             "Positive" ="Positive  (rg > 0)")
  ) +
  scale_fill_manual(
    values=pal_group,
    name="Trait group"
  ) +
  scale_linewidth_continuous(range=c(0.5,3.5), guide="none") +
  
  coord_fixed(xlim=c(-2.2, 3.2), ylim=c(-2.0, 2.0)) +
  theme_void(base_family=FONT_BASE) +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size=9.5, face="bold", color="grey10"),
    legend.text      = element_text(size=8.5, color="grey25"),
    legend.key.size  = unit(0.4,"cm"),
    legend.spacing.y = unit(0.3,"cm"),
    legend.box.margin= margin(l = 25), 
    plot.background  = element_rect(fill="white", color=NA),
    panel.background = element_rect(fill="white", color=NA),
    plot.margin      = margin(10,10,10,10)
  ) +
  guides(
    color=guide_legend(order=1,
                       override.aes=list(linewidth=1.8)),
    fill =guide_legend(order=2,
                       override.aes=list(size=6)) 
  )

dir.create(here("06_manuscrito/figuras"), showWarnings=FALSE, recursive=TRUE)

ggsave(here("06_manuscrito/figuras/fig4_rg_v3.pdf"),
       fig4, width=20, height=14, units="cm", dpi=300)
ggsave(here("06_manuscrito/figuras/fig4_rg_v3.png"),
       fig4, width=20, height=14, units="cm", dpi=300, bg="white")

cat("✓ Figure 4 network correlation chart exported successfully.\n")

# =============================================================================
# FIGURE 5 — ELITE COHORT SELECTION RANKING BAR CHART (IR_final)
# =============================================================================

ranking_plot <- ranking |>
  mutate(
    individuo = as.character(individuo),
    familia_f = factor(familia, levels=c(1,2,3)),
    consenso  = individuo %in% consenso_ids
  ) |>
  arrange(desc(IR_final)) |>
  mutate(rank_plot = row_number())

# Isolate Top 30 cohort for graphical legibility
top30 <- ranking_plot |> head(30)

fig5 <- ggplot(top30,
               aes(x=reorder(individuo, IR_final),
                   y=IR_final)) +
  
  # Error bars mapping using prediction variance parameters as proxy links
  geom_linerange(aes(
    ymin = IR_final - 0.5*sqrt(var_pred_medio/mean(var_pred_medio)),
    ymax = IR_final + 0.5*sqrt(var_pred_medio/mean(var_pred_medio)),
    color = familia_f
  ), linewidth=0.8, alpha=0.4) +
  
  # Primary bars mapping
  geom_col(aes(fill=familia_f, alpha=consenso),
           width=0.72, color=NA) +
  
  # Recommended candidates checklist markers mapping
  geom_point(data=filter(top30, consenso),
             aes(color=familia_f),
             shape=8, size=2.5, stroke=0.8,
             position=position_nudge(y=0.05)) +
  
  scale_fill_manual(values=PAL_FAM,
                    labels=c("1"="293 × 355","2"="501 × 355","3"="516 × 355"),
                    name="Cross") +
  scale_color_manual(values=PAL_FAM, guide="none") +
  scale_alpha_manual(values=c("TRUE"=1,"FALSE"=0.65), guide="none") +
  
  geom_hline(yintercept=0, color="grey30", linewidth=0.35) +
  
  coord_flip() +
  labs(x=NULL, y=expression(IR[final])) +
  theme_q1() +
  theme(axis.text.y=element_text(size=7.5),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_line(color="grey90"))

ggsave(here("06_manuscrito/figuras/fig5_ranking.pdf"),
       fig5, width=12, height=16, units="cm", dpi=300)
ggsave(here("06_manuscrito/figuras/fig5_ranking.png"),
       fig5, width=12, height=16, units="cm", dpi=300)
cat("✓ Figure 5 final resilience index ranking bar chart generated.\n")

# =============================================================================
# FIGURE 6 — INTEGRATED COMPONENT PANEL: β_POM SLOPE VS WEIGHTED HERITABILITY
# =============================================================================

params_plot <- params |>
  mutate(
    trait_f    = factor(trait, levels=names(trait_labels),
                        labels=trait_labels),
    sig        = case_when(
      P_beta_neg > 0.95 ~ "P(β<0) > 0.95",
      P_beta_neg < 0.05 ~ "P(β>0) > 0.95",
      TRUE              ~ "Inconclusive"
    ),
    beta_sc    = beta_POM / max(abs(beta_POM), na.rm=TRUE)
  )

# ── Panel 6A — Fixed Disease Covariate Effects (β_POM) ──────────────────────
p6a <- ggplot(params_plot,
              aes(x=reorder(trait_f, beta_POM), y=beta_POM)) +
  
  geom_col(aes(fill=sig), width=0.50, color=NA) +
  geom_hline(yintercept=0, color="grey30", linewidth=0.4) +
  
  scale_fill_manual(
    values=c(
      "P(β<0) > 0.95" = "#C0392B",
      "P(β>0) > 0.95" = "#1E8449",
      "Inconclusive"  = "#AAB7B8"
    ),
    name=NULL
  ) +
  
  coord_flip() +
  labs(x=NULL, y=expression(beta[POM])) +
  theme_q1() +
  theme(grid.major.y=element_blank(),
        legend.position="bottom",
        legend.text=element_text(size=7.5))

# ── Panel 6B — Marginal Integrated Heritabilities (Weighted Mean h²) ────────
p6b <- ggplot(params_plot,
              aes(x=reorder(trait_f, h2_pond), y=h2_pond)) +
  
  geom_col(fill="#2471A3", alpha=0.85, width=0.50) +
  
  # Arbitrary reference line baseline marker mapping
  geom_hline(yintercept=0.3, color="#E74C3C", linewidth=0.5,
             linetype="dashed") +
  
  annotate("text", x=0.6, y=0.31, label="h² = 0.30",
           size=2.5, color="#E74C3C", family=FONT_BASE) +
  
  scale_y_continuous(limits=c(0, 0.55), expand=c(0,0)) +
  coord_flip() +
  labs(x=NULL, y=expression(italic(h)^2~(weighted~mean))) +
  theme_q1() +
  theme(panel.grid.major.y=element_blank(),
        axis.text.y=element_blank())

fig6 <- p6a + p6b +
  plot_layout(widths=c(1.2, 1), guides="collect") &
  theme(legend.position="bottom")

ggsave(here("06_manuscrito/figuras/fig6_beta_h2.pdf"),
       fig6, width=16, height=8, units="cm", dpi=300)
ggsave(here("06_manuscrito/figuras/fig6_beta_h2.png"),
       fig6, width=16, height=8, units="cm", dpi=300)
cat("✓ Figure 6 integrated panel (POM slope vs heritability means) generated successfully.\n")


# =============================================================================
# APPENDIX PLOT A — EPIDEMIOLOGICAL TIMING ANALYSIS (Kaplan-Meier Curves)
# Target event: Time elapsed to reach critical CABMV score ≥ 3
# =============================================================================

library(survival)
library(ggplot2)
library(dplyr)
library(here)

FONT_BASE   <- "serif"
PAL_FAM     <- c("1"="#1B4F72","2"="#D4AC0D","3"="#C0392B")
FAM_LABELS  = c("1"="FBC₁ (293 \u00d7 355)",
                "2"="FBC₁ (501 \u00d7 355)",
                "3"="FBC₁ (516 \u00d7 355)")

theme_q1 <- function(base_size=11) {
  theme_classic(base_family=FONT_BASE, base_size=base_size) +
    theme(
      axis.line        = element_line(color="grey25", linewidth=0.45),
      axis.ticks       = element_line(color="grey25", linewidth=0.3),
      axis.text        = element_text(color="grey20", size=base_size-1),
      axis.title       = element_text(color="grey10", size=base_size,
                                      face="bold"),
      panel.grid.major = element_line(color="grey93", linewidth=0.25),
      panel.grid.minor = element_blank(),
      legend.title     = element_text(size=base_size-1, face="bold"),
      legend.text      = element_text(size=base_size-1),
      legend.key.size  = unit(0.45,"cm"),
      legend.position  = "right",
      plot.background  = element_rect(fill="white", color=NA),
      panel.background = element_rect(fill="white", color=NA),
      plot.margin      = margin(8,8,8,8)
    )
}

# ── Load and align survival datasets frames ──────────────────────────────────
bdsip_std   <- readRDS(here("01_dados/processados/bdsip_completo.rds"))
dataset     <- as.data.frame(readRDS(here("01_dados/processados/dataset_pom.rds")))

familia_info <- dataset |>
  select(individuo, familia) |>
  distinct() |>
  mutate(individuo = as.integer(as.character(individuo)))

surv_df <- bdsip_std |>
  select(-familia) |> 
  mutate(
    individuo  = as.integer(as.character(individuo)),
    t_surv_obs = as.numeric(t_surv_obs),
    evento     = as.integer(evento)
  ) |>
  left_join(familia_info, by = "individuo") |>
  filter(!is.na(familia)) |>
  mutate(familia_f = factor(familia, levels = c(1, 2, 3)))

cat(sprintf("Survival Sample Profile: %d lines | Events reached: %d | Right-censored data: %d\n",
            nrow(surv_df), sum(surv_df$evento), sum(surv_df$evento==0)))

# ── Fit Kaplan-Meier Non-Parametric Estimators ───────────────────────────────
km_fit <- survfit(Surv(t_surv_obs, evento) ~ familia_f,
                  data    = surv_df,
                  conf.type = "log") 

# Execute log-rank non-parametric significance test
lr_test <- survdiff(Surv(t_surv_obs, evento) ~ familia_f, data=surv_df)
chi2 <- lr_test$chisq
pval <- 1 - pchisq(chi2, df=length(unique(surv_df$familia_f))-1)
cat(sprintf("\nLog-rank Test Significance: χ²=%.3f | p-value=%.3f\n", chi2, pval))

# ── Extract tracking vectors manually to format step structures ──────────────
km_data <- do.call(rbind, lapply(1:3, function(fam) {
  idx   <- which(names(km_fit$strata) == paste0("familia_f=", fam))
  start <- if(idx==1) 1 else sum(km_fit$strata[1:(idx-1)])+1
  end   <- sum(km_fit$strata[1:idx])
  data.frame(
    familia_f = factor(fam, levels=c(1,2,3)),
    time      = c(0, km_fit$time[start:end]),
    surv      = c(1, km_fit$surv[start:end]),
    lower     = c(1, km_fit$lower[start:end]),
    upper     = c(1, km_fit$upper[start:end])
  )
}))

# Align markers for right-censored observations (+)
cens_df <- surv_df |>
  filter(evento == 0) |>
  mutate(
    surv_at = mapply(function(t, fam) {
      sub <- km_data |> filter(familia_f == fam)
      approx(x = sub$time, y = sub$surv, xout = t,
             method = "constant", rule = 2)$y
    }, t_surv_obs, familia_f)
  )

mediana_fic3 <- 280

# ── Plotting Pipeline Execution ──────────────────────────────────────────────
fig_km <- ggplot() +
  
  geom_ribbon(data=km_data,
              aes(x=time, ymin=lower, ymax=upper,
                  fill=familia_f),
              alpha=0.15) +
  
  geom_step(data=km_data,
            aes(x=time, y=surv, color=familia_f),
            linewidth=0.75) +
  
  geom_point(data=cens_df, 
             aes(x=t_surv_obs, y=surv_at, color=familia_f),
             shape=3, size=2.0, stroke=0.8, alpha=0.75) +
  
  geom_vline(xintercept=mediana_fic3,
             color=PAL_FAM["3"], linewidth=0.45,
             linetype="dashed") +
  
  annotate("text",
           x=mediana_fic3 + 12, y=0.62,
           skinny_label=sprintf("%d days", mediana_fic3),
           color=PAL_FAM["3"], size=3.1,
           family=FONT_BASE, hjust=0, fontface="italic") +
  
  annotate("text",
           x=460, y=0.95, 
           label=sprintf("Log-rank\nχ² = %.3f\np = %.3f", chi2, pval),
           size=2.9, family=FONT_BASE, hjust=1,
           color="grey30", lineheight=1.1) +
  
  scale_color_manual(values=PAL_FAM, labels=FAM_LABELS,
                     name="Full-sib family") +
  scale_fill_manual (values=PAL_FAM, labels=FAM_LABELS,
                     name="Full-sib family") +
  scale_x_continuous(limits=c(0, 480),
                     breaks=seq(0, 450, 90),
                     expand=c(0.01, 0)) +
  scale_y_continuous(limits=c(0, 1.02),
                     breaks=seq(0, 1, 0.2),
                     labels=scales::percent_format(accuracy=1),
                     expand=c(0.01, 0)) +
  
  labs(x="Time (days after transplanting)",
       y="Proportion below severity threshold") +
  
  theme_q1() +
  theme(
    legend.position   = c(0.25, 0.25), 
    legend.background = element_rect(fill=alpha("white", 0.85),
                                     color="grey80", linewidth=0.3)
  )

dir.create(here("06_manuscrito/figuras"), showWarnings=FALSE, recursive=TRUE)

ggsave(here("06_manuscrito/figuras/fig_kaplan_meier.pdf"),
       fig_km, width=14, height=10, units="cm", dpi=300, 
       device = cairo_pdf) 

ggsave(here("06_manuscrito/figuras/fig_kaplan_meier.png"),
       fig_km, width=14, height=10, units="cm", dpi=300, bg="white")

cat("✓ Appendix Figure A (Kaplan-Meier survival plot) exported successfully.\n")


# =============================================================================
# APPENDIX PLOT B — INDIVIDUAL CHRONIC PATHOSYSTEM PROBABILITIES (POM Summary)
# Individual-level posterior tracking partitioned by cross structures
# =============================================================================

library(ggplot2)
library(dplyr)
library(ggrepel)
library(here)

FONT_BASE  <- "serif"
PAL_FAM    <- c("1"="#1B4F72","2"="#D4AC0D","3"="#C0392B")
FAM_LABELS <- c("1"="FBC₁ (293 \u00d7 355)",
                "2"="FBC₁ (501 \u00d7 355)",
                "3"="FBC₁ (516 \u00d7 355)")

CONSENSO_IDS <- c("97","9","131","2","83","223")

theme_q1 <- function(base_size=11) {
  theme_classic(base_family=FONT_BASE, base_size=base_size) +
    theme(
      axis.line        = element_line(color="grey25", linewidth=0.45),
      axis.ticks       = element_line(color="grey25", linewidth=0.3),
      axis.text        = element_text(color="grey20", size=base_size-2),
      axis.title       = element_text(color="grey10", size=base_size-1,
                                      face="bold"),
      panel.grid.major = element_line(color="grey93", linewidth=0.25),
      panel.grid.minor = element_blank(),
      legend.title     = element_text(size=base_size-1, face="bold"),
      legend.text      = element_text(size=base_size-1),
      legend.key.size  = unit(0.40,"cm"),
      plot.background  = element_rect(fill="white", color=NA),
      panel.background = element_rect(fill="white", color=NA),
      plot.margin      = margin(8,10,8,8)
    )
}

# ── Load and compile dataset matrices ────────────────────────────────────────
pom_ind  <- readRDS(here("01_dados/processados/pom_indices.rds"))
dataset  <- as.data.frame(readRDS(here("01_dados/processados/dataset_pom.rds")))

familia_info <- dataset |>
  select(individuo, familia, cruzamento) |>
  distinct() |>
  mutate(individuo = as.character(individuo))

pom_plot <- pom_ind |>
  select(-familia) |> 
  mutate(individuo = as.character(individuo)) |>
  left_join(familia_info, by="individuo") |>
  filter(!is.na(familia)) |>
  mutate(
    familia_f  = factor(familia, levels=c(1,2,3)),
    consenso   = individuo %in% CONSENSO_IDS,
    rank_within = ave(P_severo_media, familia_f,
                      FUN=function(x) rank(x, ties.method="first"))
  ) |>
  arrange(familia_f, rank_within) |>
  mutate(x_pos = row_number())   

media_geral <- mean(pom_plot$P_severo_media, na.rm=TRUE)

# Isolate cross partitions intervals coordinates for grid dividers
fam_breaks <- pom_plot |>
  group_by(familia_f) |>
  summarise(xmax=max(x_pos), xmin=min(x_pos), xmid=mean(x_pos),
            .groups="drop")

# ── Plotting Pipeline Execution ──────────────────────────────────────────────
fig_pom <- ggplot(pom_plot, aes(x=x_pos, y=P_severo_media)) +
  
  # Structural segment dividers mapping
  geom_vline(data=fam_breaks[-nrow(fam_breaks),],
             aes(xintercept=xmax + 0.5),
             color="grey70", linewidth=0.4, linetype="solid") +
  
  # Color palette shading boxes per tracking family
  annotate("rect",
           xmin=fam_breaks$xmin[1]-0.5, xmax=fam_breaks$xmax[1]+0.5,
           ymin=-Inf, ymax=Inf,
           fill=PAL_FAM["1"], alpha=0.04) +
  annotate("rect",
           xmin=fam_breaks$xmin[2]-0.5, xmax=fam_breaks$xmax[2]+0.5,
           ymin=-Inf, ymax=Inf,
           fill=PAL_FAM["2"], alpha=0.04) +
  annotate("rect",
           xmin=fam_breaks$xmin[3]-0.5, xmax=fam_breaks$xmax[3]+0.5,
           ymin=-Inf, ymax=Inf,
           fill=PAL_FAM["3"], alpha=0.04) +
  
  # Population mean horizontal tracker marker mapping
  geom_hline(yintercept=media_geral,
             color="grey40", linewidth=0.45, linetype="dashed") +
  annotate("text",
           x=1, y=media_geral + 0.004,
           label=sprintf("Population mean = %.3f", media_geral),
           size=2.8, family=FONT_BASE, color="grey35",
           hjust=0, fontface="italic") +
  
  # Phenotypic data distributions array
  geom_point(data=filter(pom_plot, !consenso),
             aes(color=familia_f),
             size=2.0, alpha=0.70, shape=16) +
  
  # Recommended Elite Consensus Targets mapping (Golden star markers)
  geom_point(data=filter(pom_plot, consenso),
             aes(color=familia_f),
             shape=8, size=3.5, stroke=1.0,
             color="#B7950B", alpha=1) +
  
  # Labels overlay mappings via repelling algorithms
  geom_label_repel(
    data=filter(pom_plot, consenso),
    aes(label=individuo, color=familia_f),
    size=2.8, fontface="bold", family=FONT_BASE,
    label.padding=unit(0.12,"lines"), label.size=0.20,
    box.padding=0.5, point.padding=0.3,
    force=5, min.segment.length=0.2,
    segment.color="grey45", segment.size=0.35,
    fill=alpha("white",0.90),
    show.legend=FALSE, max.overlaps=Inf, seed=42,
    direction="y"
  ) +
  
  # Family specific sub-labels under the primary X-axis alignment
  annotate("text",
           x=fam_breaks$xmid, y=-0.012,
           label=c("FBC₁","FBC₂","FBC₃"),
           color=PAL_FAM[as.character(1:3)],
           size=3.2, family=FONT_BASE, fontface="bold") +
  
  scale_color_manual(values=PAL_FAM, labels=FAM_LABELS,
                     name="Full-sib family") +
  scale_y_continuous(limits=c(-0.02, max(pom_plot$P_severo_media)*1.15),
                     breaks=seq(0, 0.25, 0.05),
                     expand=c(0,0)) +
  scale_x_continuous(breaks=NULL, expand=c(0.01,0)) +
  
  labs(x=NULL,
       y=expression(bar(italic(P))(score >= 3))) +
  
  theme_q1() +
  theme(
    axis.line.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    legend.position  = "right"
  )

# ── Save Checkpoint Verification ─────────────────────────────────────────────
dir.create(here("06_manuscrito/figuras"), showWarnings=FALSE, recursive=TRUE)

ggsave(here("06_manuscrito/figuras/fig_pom_individuos.pdf"),
       fig_pom, width=18, height=10, units="cm", dpi=300)
ggsave(here("06_manuscrito/figuras/fig_pom_individuos.png"),
       fig_pom, width=18, height=10, units="cm", dpi=300, bg="white")

cat("✓ Appendix Figure B (Individual POM tracking curves) exported successfully.\n")
cat("All figures are compiled inside: 06_manuscrito/figuras/\n")