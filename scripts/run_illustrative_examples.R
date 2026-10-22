# *********************************************
# Authors: Aurore Archimbaud and Klaus Nordhausen
#          TBS Business School
# 20/10/2026
# *********************************************

# 1st example: Figure 1 ----------
# Comparison of whitening and PCA 
# Simulated centered bivariate Gaussian dataset (n = 250) with non-identity
# covariance matrix. 

# load necessary packages
library(MASS)
library(ggplot2)

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

# 4. pca
cov_centered <- cov(x_centered)
eig_pca <- eigen(cov_centered)
rotation <- eig_pca$vectors
x_pca <- x_centered %*% rotation

# 5. combine the three versions into one long data frame
panel_labels <- c("Observed data", "After whitening", "After PCA")

df <- rbind(
  data.frame(x = x_centered[, 1], y = x_centered[, 2], panel = panel_labels[1]),
  data.frame(x = x_whitened[, 1], y = x_whitened[, 2], panel = panel_labels[2]),
  data.frame(x = x_pca[, 1],      y = x_pca[, 2],      panel = panel_labels[3])
)
df$panel <- factor(df$panel, levels = panel_labels)  # keep panel order

# 6. plotting: one plot, three facets
axis_length <- 8
axis_arrow <- arrow(length = unit(0.1, "inches"), ends = "both")

p <- ggplot(df, aes(x, y)) +
  # fixed-length axes (drawn in every facet)
  annotate("segment", x = -axis_length, xend = axis_length, y = 0, yend = 0,
           colour = "gray", arrow = axis_arrow) +
  annotate("segment", x = 0, xend = 0, y = -axis_length, yend = axis_length,
           colour = "gray", arrow = axis_arrow) +
  geom_point(size = 1.2) +
  facet_wrap(~ panel, nrow = 1) +
  coord_fixed(xlim = c(-axis_length, axis_length),
              ylim = c(-axis_length, axis_length)) +
  theme_void() +
  theme(
    strip.text = element_text(size = 14,
                              margin = margin(b = 6)),
    panel.spacing = unit(1, "lines")
  )

ggsave("figures/whitenpca.pdf", p, width = 12, height = 4)

# 2nd example: Figure 3 ----------
# Simulated dataset (n=400, p=4) obtained by projecting the first two 
# coordinates of standard Gaussian observations onto a circle and then applying 
# a random linear transformation.

library(BSSprep)
library(ICS)
library(GGally)
library(ggplot2)
library(patchwork)

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

# Generate random mixing matrix
A <- matrix(rnorm(p^2), p, p)

# Mix the data (linear transformation)
X <- X1 %*% t(A)
colnames(X) <- paste("X", 1:4, sep = ".")

# Whitening (sphering) step
White <- BSSprep(X)$Y
colnames(White) <- paste("W", 1:4, sep = ".")

# PCA transformation
PCA <- princomp(X)$scores
colnames(PCA) <- paste("PC", 1:4, sep = ".")

# Invariant Coordinate Selection (ICS)
ICs <- components(ICS(X))
colnames(ICs) <- paste("IC", 1:4, sep = ".")


### Visualisation ----------

# Helper: one scatterplot matrix with its panel header inside the plot
make_pairs <- function(data, title) {
  p <- ggpairs(as.data.frame(data),
               title = title,
               lower = list(continuous = wrap("points", size = 0.5))) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(size = 13, hjust = 0.5),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  # convert the ggmatrix into a grob so that patchwork can arrange it
  wrap_elements(full = ggmatrix_gtable(p))
}

plotX     <- make_pairs(X,     "Observed data")
plotWhite <- make_pairs(White, "After whitening")
plotPCA   <- make_pairs(PCA,   "After PCA")
plotICS   <- make_pairs(ICs,   "After ICS")

# One figure, four panels (2 x 2)
fig_structure <- (plotX | plotWhite) / (plotPCA | plotICS)
fig_structure

ggsave("figures/Structure.pdf", plot = fig_structure, width = 12, height = 12)



# 3rd example: Figure 6 ------------------------------------------------------------

# Simulated example for Section 3.4: ICS as an ICA method
# Two sources (Laplace and t_6) have the same classical kurtosis, so that
# COV-COV4 (FOBI) cannot separate them, while COV-MLC can.

library(ICS)
library(JADE)
library(ggplot2)

# standardized independent sources with kurtoses 3, 3, 0 and -1.2
rsources <- function(n) {
  laplace <- rexp(n) * sample(c(-1, 1), n, replace = TRUE)
  scale(cbind(laplace, rt(n, df = 6), rnorm(n), runif(n)))
}
src_names <- c("Laplace", "t[6]", "Normal", "Uniform")

pairs <- list(
  "COV-COV[4]" = list(S1 = ICS_cov, S2 = ICS_cov4, S2_args = list()),
  "COV-MLC"    = list(S1 = ICS_cov, S2 = ICS_tM,   S2_args = list(df = 1))
)
fit_pair <- function(X, p) ICS(X, S1 = p$S1, S2 = p$S2, S2_args = p$S2_args)

### Visualization: absolute correlations between ICS components and sources -----
set.seed(2026)
n <- 2000
S <- rsources(n)
A <- matrix(rnorm(16), 4, 4)
X <- tcrossprod(S, A)

dat <- do.call(rbind, lapply(names(pairs), function(nm) {
  fit <- fit_pair(X, pairs[[nm]])
  print(round(gen_kurtosis(fit), 3))
  cat(nm, "MD:", round(MD(coef(fit), A), 3), "\n")
  C <- abs(cor(components(fit), S))
  data.frame(component = factor(rep(paste0("IC.", 1:4), 4), levels = paste0("IC.", 4:1)),
             source = factor(rep(src_names, each = 4), levels = src_names),
             value = as.vector(C),
             pair = nm)
}))
dat$pair <- factor(dat$pair, levels = names(pairs))

p <- ggplot(dat, aes(x = source, y = component, fill = value)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3.2) +
  scale_fill_gradient(low = "white", high = "grey40", limits = c(0, 1),
                      name = "|correlation|") +
  scale_x_discrete(labels = scales::label_parse()) +
  facet_wrap(~ pair, labeller = label_parsed) +
  labs(x = "Source", y = "ICS component") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank())
ggsave("figures/ica_example.pdf", p, width = 7, height = 3.2)

### Separation accuracy (minimum distance index) -----
set.seed(1)
nrep <- 100
res <- sapply(c(1000, 5000, 20000), function(n) {
  rowMeans(replicate(nrep, {
    S <- rsources(n)
    A <- matrix(rnorm(16), 4, 4)
    X <- tcrossprod(S, A)
    sapply(pairs, function(p) MD(coef(fit_pair(X, p)), A))
  }))
})
colnames(res) <- c("n=1000", "n=5000", "n=20000")
round(res, 3)