# *********************************************
# Authors: Aurore Archimbaud 
#          TBS Business School
# 29/05/2026
# *********************************************


# Libraries --------------------------------------------------------------
# For crabs data
library(cellWise)

# For ICS
library(ICSOutlier)
library(ICSClust)
library(ICSOutlier)

# For PCA
library(rrcov)

# For plotting
library(ggplot2)
library(ggrepel)
library(magrittr)
library(ggthemes)

# Data  --------------------------------------------------------------------
data("crabs", package = "MASS")

# Log transformation, no standardization
data_crabs <- apply(crabs[,4:8], 2, log)

# Define clusters based on sex and species
clusters <-  factor(paste(crabs$sp, crabs$sex, sep = "_"))



# ICS with TCOV-COV ---------------------------------------------------------------------
ics_res <- ICS(data_crabs, S1 = ICS_tcov, S2 = ICS_cov)

# Selection of components based on screeplot
plot_ics_eigenvalues <- select_plot(ics_res, type = "lines") +
  theme_minimal() +
  ylim(0,4)+
  labs(x = "ICS-Components", fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position = "bottom") 


pdf("figures/crabs_ics_eigenvalues.pdf",  width = 6.5, height = 5.75)
print(plot_ics_eigenvalues)
dev.off()


# Plot the selected components: IC1 and IC2
df_ics <- data.frame(components(ics_res), clusters = clusters)
colors <- ggthemes::colorblind_pal()(4)[c(3, 1, 2, 4)]

plot_ics <- df_ics %>% 
  ggplot(aes(x= IC.1, y = IC.2, color = clusters)) +
  geom_point() +
  scale_color_manual(values = colors) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position ="bottom") 
plot_ics
pdf("figures/crabs_clus_ics.pdf",  width = 6.5, height = 5.75)
print(plot_ics)
dev.off()

# Run kmeans clustering method with 4 clusters
ics_clusters <- kmeans_clust(df_ics[,1:2], k = 4, clusters_only = TRUE)

# Compute Adjusted Rand Index
mclust::adjustedRandIndex(clusters, ics_clusters)



# PCA ---------------------------------------------------------------------
PCA  <- PcaClassic(data_crabs, scale = TRUE)

# Selection of components based on screeplot
df_PCA_eigenvalues <- data.frame(eigenvalues = PCA@eigenvalues, 
                                 PC = colnames(PCA@scores))

plot_pca_eigenvalues <- df_PCA_eigenvalues %>% 
  ggplot(aes(x = PC, y = eigenvalues, group = 0)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(x = "PCA-Components", y = "Eigenvalues", fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position="bottom") 

pdf("figures/crabs_pca_eigenvalues.pdf",  width = 6.5, height = 5.75)
print(plot_pca_eigenvalues)
dev.off()

# Keeping only PC1
df_pca <- data.frame(data.frame(PCA@scores), clusters = clusters)

# Run kmeans clustering method with 4 clusters
pca_clusters <- kmeans_clust(df_pca[,1], k = 4, 
                             clusters_only = TRUE)

# Compute Adjusted Rand Index
mclust::adjustedRandIndex(clusters, pca_clusters)


# Keeping PC1 and PC2
# Plot the selected components
plot_pca <- df_pca %>% 
  ggplot(aes(x = PC2, y = PC3, color = clusters)) +
  geom_point() +
  scale_color_manual(values = colors) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position="bottom") 

pdf("figures/crabs_clus_pca.pdf", width = 6.5, height = 5.75)
print(plot_pca)
dev.off()

# Run kmeans clustering method with 4 clusters
pca_clusters <- kmeans_clust(df_pca[,2:3], k = 4, clusters_only = TRUE)

# Compute Adjusted Rand Index
mclust::adjustedRandIndex(clusters, pca_clusters)


# Run kmeans clustering method with 4 clusters on PC1, PC2 and PC3
pca_clusters <- kmeans_clust(df_pca[,1:3], k = 4, clusters_only = TRUE)

# Compute Adjusted Rand Index
mclust::adjustedRandIndex(clusters, pca_clusters)

