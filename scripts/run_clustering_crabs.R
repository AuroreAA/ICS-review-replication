# *********************************************
# Authors: Aurore Archimbaud
#          TBS Business School
# 20/10/2026
# *********************************************
# Reproduce Figure 4

# Libraries --------------------------------------------------------------
# For ICS
library(ICSOutlier)
library(ICSClust)

# For PCA
library(rrcov)

# For plotting
library(ggplot2)
library(magrittr)
library(ggthemes)
library(patchwork)

# Data  --------------------------------------------------------------------
data("crabs", package = "MASS")

# Log transformation, no standardization
data_crabs <- apply(crabs[, 4:8], 2, log)

# Define clusters based on sex and species
clusters <- factor(paste(crabs$sp, crabs$sex, sep = "_"))

# Common colour palette (colour-blind friendly)
colors <- ggthemes::colorblind_pal()(4)[c(3, 1, 2, 4)]


# ICS with TCOV-COV ---------------------------------------------------------
ics_res <- ICS(data_crabs, S1 = ICS_tcov, S2 = ICS_cov)

df_ics <- data.frame(components(ics_res), clusters = clusters)

# Run kmeans clustering method with 4 clusters on IC.1 and IC.2
ics_clusters <- kmeans_clust(df_ics[, 1:2], k = 4, clusters_only = TRUE)
mclust::adjustedRandIndex(clusters, ics_clusters)


# PCA ---------------------------------------------------------------------
PCA <- PcaClassic(data_crabs, scale = TRUE)

# Name the components PC.1, PC.2, ... (same style as IC.1, IC.2, ...)
pca_scores <- PCA@scores
colnames(pca_scores) <- paste0("PC.", seq_len(ncol(pca_scores)))

df_pca <- data.frame(pca_scores, clusters = clusters)

# kmeans on PC.1 only
pca_clusters <- kmeans_clust(df_pca[, 1], k = 4, clusters_only = TRUE)
mclust::adjustedRandIndex(clusters, pca_clusters)

# kmeans on PC.2 and PC.3
pca_clusters <- kmeans_clust(df_pca[, 2:3], k = 4, clusters_only = TRUE)
mclust::adjustedRandIndex(clusters, pca_clusters)

# kmeans on PC.1, PC.2 and PC.3
pca_clusters <- kmeans_clust(df_pca[, 1:3], k = 4, clusters_only = TRUE)
mclust::adjustedRandIndex(clusters, pca_clusters)


# Visualisation ----------------------------------------------------------

# Common theme for all panels: centred title (method, top row only)
# and centred subtitle (what each panel shows)
theme_paper <- theme_minimal(base_size = 16) +
  theme(plot.title    = element_text(hjust = 0.5,  size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 15),
        legend.position = "bottom")

# 1. Screeplots: same y scale (0-5); each method keeps its own y label
make_screeplot <- function(component, value, title, ylab) {
  df <- data.frame(component = factor(component, levels = component),
                   value     = value)
  ggplot(df, aes(x = component, y = value, group = 1)) +
    geom_point() +
    geom_line() +
    scale_y_continuous(limits = c(0, 5)) +
    labs(title = title, x = "", y = ylab) +
    theme_paper
}

plot_eig_ics <- make_screeplot(colnames(components(ics_res)),
                               ics_res$gen_kurtosis,
                               title    = expression("ICS (TCOV"[2]*"-COV)"),
                               ylab     = "Generalized eigenvalues")

plot_eig_pca <- make_screeplot(colnames(pca_scores),
                               PCA@eigenvalues,
                               title    = "PCA",
                               ylab     = "Eigenvalues")

# 2. Scatterplots of the selected components (no title: the method is
#    already given at the top of each column)
plot_ics <- ggplot(df_ics, aes(x = IC.1, y = IC.2, color = clusters)) +
  geom_point() +
  scale_color_manual(values = colors) +
  labs(color = "Group") +
  theme_paper

plot_pca <- ggplot(df_pca, aes(x = PC.2, y = PC.3, color = clusters)) +
  geom_point() +
  scale_color_manual(values = colors) +
  labs(color = "Group") +
  theme_paper

# 3. One figure: ICS in the left column, PCA in the right column,
#    screeplots on top, scatterplots below, one legend
fig_crabs <- plot_eig_ics + plot_eig_pca +
  plot_ics     + plot_pca +
  plot_layout(ncol = 2, guides = "collect", heights = c(1, 1.3)) &
  theme(legend.position = "bottom")
fig_crabs

ggsave("figures/crabs_ics_pca.pdf", plot = fig_crabs, width = 13, height = 11)