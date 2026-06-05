# *********************************************
# Authors: Aurore Archimbaud 
#          TBS Business School
# 05/06/2026
# *********************************************


# 1st example ----------

# load necessary package
library(MASS)

# 1. simulate data
set.seed(1234)
n <- 250
mu <- c(0, 0)
sigma <- matrix(c(4, 1.5, 
                  1.5, 2), nrow = 2)
x <- mvrnorm(n, mu, sigma)

# 2. center the data
x_mean <- colMeans(x)
x_centered <- sweep(x, 2, x_mean)

# 3. whitening
s <- cov(x)
eig <- eigen(s)
s_inv_sqrt <- eig$vectors %*% diag(1 / sqrt(eig$values)) %*% t(eig$vectors)
x_whitened <- t(s_inv_sqrt %*% t(x_centered))

# 4. pca (corrected!)
cov_centered <- cov(x_centered)
eig_pca <- eigen(cov_centered)
rotation <- eig_pca$vectors
x_pca <- x_centered %*% rotation

# 5. axis drawing helper with fixed length
draw_fixed_axis <- function(length = 2) {
  arrows(-length, 0, length, 0, length = 0.1, col = "gray")
  arrows(0, -length, 0, length, length = 0.1, col = "gray")
}

# 6. plotting

pdf("figures/whitenpca.pdf", paper = "special", width = 12, height = 4)
par(mfrow = c(1, 3), mar = c(1, 1, 1, 1)/10)

# panel 1: centered original data
plot(NA, asp = 1, xlab = "", ylab = "", xaxt = 'n', yaxt = 'n', 
     frame.plot = FALSE,
     xlim = c(-8, 8), ylim = c(-8, 8))
draw_fixed_axis(length = 8)
points(x_centered, pch = 16, cex = 0.8)

# panel 2: whitened data
plot(NA, asp = 1, xlab = "", ylab = "", xaxt = 'n', yaxt = 'n', frame.plot = FALSE,
     xlim = c(-8, 8), ylim = c(-8, 8))
draw_fixed_axis(length = 8)
points(x_whitened, pch = 16, cex = 0.8)

# panel 3: pca-rotated data
plot(NA, asp = 1, xlab = "", ylab = "", xaxt = 'n', yaxt = 'n', frame.plot = FALSE,
     xlim = c(-8, 8), ylim = c(-8, 8))
draw_fixed_axis(length = 8)
points(x_pca, pch = 16, cex = 0.8)
dev.off()



# 2nd example ----------

# Set seed for reproducibility
set.seed(1234)

# Define sample size and dimension
n <- 400   # number of observations
p <- 4     # number of variables

# Generate standard normal data
X0 <- matrix(rnorm(n * p), n, p)

# Copy data
X1 <- X0

# Impose non-Gaussian structure on first two variables:
# project points onto a circle (radius scaled)
R2 <- sqrt(2 / rowSums(X1[, 1:2]^2))
X1[, 1:2] <- R2 * X1[, 1:2]

# Check covariance 
cov(X1)
sqrt(mean((cov(X1) - diag(p))^2))  # deviation from identity

# Pairwise scatterplots of structured data
pairs(X1)

# Generate random mixing matrix
A <- matrix(rnorm(p^2), p, p)

# Mix the data (linear transformation)
X <- X1 %*% t(A)
colnames(X) <- paste("X", 1:4, sep = ".")

# Visualize mixed data
pairs(X)

# Whitening (sphering) step
library(BSSprep)
White <- BSSprep(X)$Y
colnames(White) <- paste("W", 1:4, sep = ".")

# Visualize whitened data
pairs(White)

# PCA transformation
PCA <- princomp(X)$scores
colnames(PCA) <- paste("PC", 1:4, sep = ".")

# Visualize PCA components
pairs(PCA)

# Invariant Coordinate Selection (ICS)
library(ICS)
ICs <- components(ICS(X))

# Visualize ICS components
pairs(ICs)


# Visualisation
library(GGally)
library(ggplot2)

# Plot: original mixed data
plotX <- ggpairs(as.data.frame(X),
                 lower = list(continuous = wrap("points", size = 0.5))) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )
plotX

# Plot: whitened data
plotWhite <- ggpairs(as.data.frame(White),
                     lower = list(continuous = wrap("points", size = 0.5))) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )
plotWhite

# Plot: PCA scores
plotPCA <- ggpairs(as.data.frame(PCA),
                   lower = list(continuous = wrap("points", size = 0.5))) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )
plotPCA

# Plot: ICS components
plotICS <- ggpairs(as.data.frame(ICs),
                   lower = list(continuous = wrap("points", size = 0.5))) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )
plotICS

ggsave("figures/StructureX.pdf", plot = plotX, width = 6, height = 6)
ggsave("figures/StructureWhite.pdf", plot = plotWhite, width = 6, height = 6)
ggsave("figures/StructurePCA.pdf", plot = plotPCA, width = 6, height = 6)
ggsave("figures/StructureICS.pdf", plot = plotICS, width = 6, height = 6)
                   


