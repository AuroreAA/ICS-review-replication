# *********************************************
# Authors: Aurore Archimbaud 
#          TBS Business School
# 05/06/2026
# *********************************************


# Load required libraries -------
# For statistical analysis
library(MASS)       # For mvrnorm (multivariate normal distribution)
library(stats)      # For cov, prcomp, lda
library(robustbase) # For covMcd (robust scatter matrix)
library(ICS)        # For ICS 
library(ICSClust)

# For visualisation
library(ggplot2) 
library(gridExtra) 
library(ggpubr)     

# Function to calculate an ellipse (from the tolEllipsePlot function)
# Arguments:
#   - loc: location estimate
#   - cov: scatter estimate
ellips <- function(loc, cov) {
  dist <- sqrt(qchisq(0.975, 2))
  A <- solve(cov)
  eA <- eigen(A, symmetric = TRUE)
  ev <- eA$values
  lambda1 <- max(ev)
  lambda2 <- min(ev)
  eigvect <- eA$vectors[, order(ev, decreasing = TRUE)[2]]
  z <- seq(0, 2 * pi, by = 0.01)
  z1 <- dist/sqrt(lambda1) * cos(z)
  z2 <- dist/sqrt(lambda2) * sin(z)
  alfa <- atan2(eigvect[2], eigvect[1])
  if (alfa < 0) alfa <- alfa + 2 * pi
  r <- matrix(c(cos(alfa), -sin(alfa),
                sin(alfa), cos(alfa)),
              ncol = 2)
  result <- loc + cbind(z1, z2) %*% r
  return(result)
}


# Data generation -------
set.seed(1234)  # For reproducibility

# Simulation of a 2-mixture of Gaussian distributions
n_group1 <- 980
n_group2 <- 20
n_variables <- 2

# First group
df_group1 <- data.frame(cluster = "Group 1", 
                        mvtnorm::rmvnorm(n = n_group1,
                                         mean = c(0, 0), 
                                         sigma = diag(c(1, 5))))
                     
# Second group
df_group2 <- data.frame(cluster = "Group 2", 
                        mvtnorm::rmvnorm(n = n_group2,
                                         mean = c(8, -2), 
                                         sigma = diag(c(1,5))))

df <- rbind(df_group1, df_group2)

# Scatter matrices -----

# Calculate scatter matrices and location vectors
cov_X <- ICS_cov(df[, 2:3])  
cov4_X <- ICS_cov4(df[, 2:3])  

# Calculate ellipses
df_ellipse_cov <- data.frame(ellips(cov_X$location, cov_X$scatter))
#df_ellipse_cov <- df_ellipse_cov[seq(1, nrow(df_ellipse_cov), by = 2), ]
df_ellipse_cov4 <- data.frame(ellips(cov4_X$location, cov4_X$scatter))
#df_ellipse_cov4 <- df_ellipse_cov4[seq(1, nrow(df_ellipse_cov4), by = 2), ]

library(ellipse)
df_ellipse_cov <- data.frame(ellipse(cov_X$scatter, cov_X$location))



# Create the plot
colors <- ggthemes::colorblind_pal()(5)[c(1, 3,  2, 4, 5)]

p <- ggplot(df, aes(x = X1, y = X2)) +
  geom_point(aes(color = cluster)) +
  geom_path(data = df_ellipse_cov, aes(x = X1, y = X2), color = "purple", linewidth = 0.1) +
  geom_path(data = df_ellipse_cov4, aes(x = X1, y = X2), color = "red", linewidth = 0.1) +
  #geom_path(data = df_ellipse_cov4, aes(x = X1, y = X2), color = "purple", linewidth = 0.1) +
  # geom_segment(data = df_lda, aes(x = X1, y = X2, xend = X1[2], yend = X2[2]),
  #              color = "cyan4", linewidth = 1) +
  # geom_segment(data = df_pca, aes(x = X1, y = X2, xend = X1[2], yend = X2[2]),
  #              color = "hotpink", linewidth = 1) +
  # geom_segment(data = df_ics, aes(x = X1, y = X2, xend = X1[2], yend = X2[2]),
  #              color = "red3", linewidth = 1) +
  # scale_color_manual(values = c("Group_1" = "dodgerblue3",
  #                               "Group_2" = "darkorange",
  #                               "COV" = "darkgreen",
  #                               "COV4" = "purple",
  #                               "LDA" = "cyan4",
  #                               "PC1" = "hotpink",
  #                               "ICS" = "red3")) +
  scale_color_manual(values = colors[1:2]) +
  theme_minimal() 
p
# +
#   theme(legend.position = "none",
#         plot.margin = margin(0, 0, 0, 0),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   coord_fixed(ratio = 1)



# Convert ellipses to data frames
df_ellipse_cov <- as.data.frame(ellipse_cov)
df_ellipse_cov$Group <- "COV"
colnames(df_ellipse_cov) <- c("X1", "X2", "Group")

df_ellipse_cov4 <- as.data.frame(ellipse_cov4)
df_ellipse_cov4$Group <- "COV4"
colnames(df_ellipse_cov4) <- c("X1", "X2", "Group")

# Calculate LDA
lda_model <- lda(cluster ~ X1 + X2, data = df)
lda_coef <- lda_model$scaling[, 1]
origin <- colMeans(df[, 1:2])
axis_length <- 20
axis_start <- origin - axis_length * lda_coef
axis_end <- origin + axis_length * lda_coef

df_lda <- data.frame(
  X1 = c(axis_start[1], axis_end[1]),
  X2 = c(axis_start[2], axis_end[2]),
  Group = "LDA"
)

# Calculate PCA
X_scaled <- scale(df[, 1:2])
pca_model <- prcomp(X_scaled, center = TRUE, scale. = TRUE)
first_component <- pca_model$rotation[, 1]
variance <- pca_model$sdev[1]^2
axis_start_pca <- colMeans(df[, 1:2]) - (axis_length-5) * first_component
axis_end_pca <- colMeans(df[, 1:2]) + (axis_length-5) * first_component

df_pca <- data.frame(
  X1 = c(axis_start_pca[1], axis_end_pca[1]),
  X2 = c(axis_start_pca[2], axis_end_pca[2]),
  Group = "PC1"
)

# Calculate ICS
ics_model <- ics(S1 = cov, S2 = covMcd, data = df[, 1:2], method = "whiten")
X_new <- ics_model$X.transformed
ics_coef <- ics_model$components[, 1]
axis_start_ics <- colMeans(df[, 1:2]) - axis_length * ics_coef
axis_end_ics <- colMeans(df[, 1:2]) + axis_length * ics_coef

df_ics <- data.frame(
  X1 = c(axis_start_ics[1], axis_end_ics[1]),
  X2 = c(axis_start_ics[2], axis_end_ics[2]),
  Group = "ICS"
)

# Combine all elements for plotting
df_plot <- rbind(
  df,
  df_ellipse_cov,
  df_ellipse_cov4,
  df_lda,
  df_pca,
  df_ics
)

# Create the plot
p <- ggplot(df_plot, aes(x = X1, y = X2, color = Group)) +
  geom_point(data = df, aes(x = X1, y = X2), alpha = 0.6) +
  geom_path(data = df_ellipse_cov, aes(x = X1, y = X2), color = "darkgreen", linewidth = 1) +
  geom_path(data = df_ellipse_cov4, aes(x = X1, y = X2), color = "purple", linewidth = 1) +
  geom_segment(data = df_lda, aes(x = X1, y = X2, xend = X1[2], yend = X2[2]),
               color = "cyan4", linewidth = 1) +
  geom_segment(data = df_pca, aes(x = X1, y = X2, xend = X1[2], yend = X2[2]),
               color = "hotpink", linewidth = 1) +
  geom_segment(data = df_ics, aes(x = X1, y = X2, xend = X1[2], yend = X2[2]),
               color = "red3", linewidth = 1) +
  scale_color_manual(values = c("Group_1" = "dodgerblue3",
                                "Group_2" = "darkorange",
                                "COV" = "darkgreen",
                                "COV4" = "purple",
                                "LDA" = "cyan4",
                                "PC1" = "hotpink",
                                "ICS" = "red3")) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.margin = margin(0, 0, 0, 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  coord_fixed(ratio = 1)

# Display the plot
print(p)

# Save the plot (optional)
# ggsave("ellipses.png", plot = p, width = 6, height = 5.25, dpi = 300)
