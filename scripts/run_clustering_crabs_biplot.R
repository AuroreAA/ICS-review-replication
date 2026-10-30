# Packages ----------------------------------------------------------------
# If needed, run once, then restart R:
# install.packages(c("ggplot2", "patchwork", "ggrepel"))

library(ICS)
library(ICSClust)
library(rrcov)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(ggrepel)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# Data --------------------------------------------------------------------
data("crabs", package = "MASS")

# Log transformation
data_crabs <- as.matrix(log(crabs[, 4:8]))

# Groups defined by colour and sex
clusters <- factor(paste(crabs$sp, crabs$sex, sep = "_"))

colors <- setNames(
  ggthemes::colorblind_pal()(4)[c(3, 1, 2, 4)],
  levels(clusters)
)

# ICS: TCOV2-COV -----------------------------------------------------------
ics_res <- ICS(
  data_crabs,
  S1 = ICS_tcov,
  S2 = ICS_cov
)

# PCA: standardized log measurements --------------------------------------
PCA <- PcaClassic(data_crabs, scale = TRUE)

# Common theme ------------------------------------------------------------
theme_paper <- theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# Screeplots --------------------------------------------------------------
make_screeplot <- function(values, prefix, title, ylab) {
  df <- data.frame(
    component = factor(
      paste0(prefix, seq_along(values)),
      levels = paste0(prefix, seq_along(values))
    ),
    value = as.numeric(values)
  )
  
  ggplot(df, aes(component, value, group = 1)) +
    geom_point() +
    geom_line() +
    scale_y_continuous(limits = c(0, 5)) +
    labs(title = title, x = NULL, y = ylab) +
    theme_paper
}

plot_eig_ics <- make_screeplot(
  ics_res$gen_kurtosis,
  prefix = "IC.",
  title = expression("ICS (TCOV"[2] * "-COV)"),
  ylab = "Generalized eigenvalues"
)

plot_eig_pca <- make_screeplot(
  PCA@eigenvalues,
  prefix = "PC.",
  title = "PCA",
  ylab = "Eigenvalues"
)

# Biplot coordinates ------------------------------------------------------
# A: observation coordinates; B: variable coordinates.
# With all axes retained, A %*% t(B) reconstructs the centred,
# standardized log measurements.

Xc <- scale(data_crabs, center = TRUE, scale = FALSE)
sd_log <- apply(data_crabs, 2, sd)
Y <- sweep(Xc, 2, sd_log, "/")

# ICS directions are stored in ROWS of W.
G <- as.matrix(ics_res$W)
A_ics <- Xc %*% t(G)
B_ics <- sweep(solve(G), 1, sd_log, "/")

# PCA coordinates
A_pca <- as.matrix(PCA@scores)
B_pca <- as.matrix(PCA@loadings)

# Check the reconstruction identities
relative_error <- function(x, y) {
  sqrt(sum((x - y)^2)) / max(1, sqrt(sum(y^2)))
}

stopifnot(
  relative_error(A_ics %*% t(B_ics), Y) < 1e-6,
  relative_error(A_pca %*% t(B_pca), Y) < 1e-6
)

# Common biplot function ---------------------------------------------------
make_biplot <- function(A, B, axes, prefix) {
  A <- A[, axes, drop = FALSE]
  B <- B[, axes, drop = FALSE]
  
  # Reciprocal scaling preserves the product A %*% t(B).
  balance <- sqrt(
    max(sqrt(rowSums(B^2))) /
      max(sqrt(rowSums(A^2)))
  )
  
  points <- data.frame(
    x = A[, 1] * balance,
    y = A[, 2] * balance,
    clusters = clusters
  )
  
  arrows <- data.frame(
    x = B[, 1] / balance,
    y = B[, 2] / balance,
    variable = colnames(data_crabs)
  )
  
  bound <- 1.22 * max(
    abs(c(points$x, points$y, arrows$x, arrows$y))
  )
  
  ggplot(points, aes(x, y)) +
    geom_hline(
      yintercept = 0, colour = "grey80", linewidth = 0.35
    ) +
    geom_vline(
      xintercept = 0, colour = "grey80", linewidth = 0.35
    ) +
    geom_point(
      aes(colour = clusters),
      size = 2,
      alpha = 0.8
    ) +
    geom_segment(
      data = arrows,
      aes(x = 0, y = 0, xend = x, yend = y),
      inherit.aes = FALSE,
      colour = "grey15",
      linewidth = 0.55,
      arrow = grid::arrow(length = grid::unit(0.18, "cm"))
    ) +
    geom_text_repel(
      data = arrows,
      aes(x, y, label = variable),
      inherit.aes = FALSE,
      seed = 1,
      min.segment.length = 0,
      max.overlaps = Inf,
      box.padding = 0.4,
      colour = "black",
      size = 5
    ) +
    scale_colour_manual(
      values = colors,
      breaks = levels(clusters),
      name = "Group"
    ) +
    guides(colour = guide_legend(nrow = 1)) +
    coord_equal(
      xlim = c(-bound, bound),
      ylim = c(-bound, bound)
    ) +
    labs(
      x = sprintf("%s%d (x %.3g)", prefix, axes[1], balance),
      y = sprintf("%s%d (x %.3g)", prefix, axes[2], balance)
    ) +
    theme_paper
}

# ICS: IC1-IC2; PCA: PC2-PC3 ----------------------------------------------
plot_ics <- make_biplot(A_ics, B_ics, axes = c(1, 2), prefix = "IC.")
plot_pca <- make_biplot(A_pca, B_pca, axes = c(2, 3), prefix = "PC.")

# Assemble: screeplots above, biplots below --------------------------------
fig_crabs <- (
  plot_eig_ics + plot_eig_pca +
    plot_ics + plot_pca +
    plot_layout(
      ncol = 2,
      guides = "collect",
      heights = c(1, 1.3)
    )
) & theme(legend.position = "bottom")

print(fig_crabs)

# Save --------------------------------------------------------------------
ggsave(
  "figures/crabs_ics_pca_biplots.pdf",
  plot = fig_crabs,
  width = 13,
  height = 11
)

ggsave(
  "figures/crabs_ics_pca_biplots.png",
  plot = fig_crabs,
  width = 13,
  height = 11,
  dpi = 300
)