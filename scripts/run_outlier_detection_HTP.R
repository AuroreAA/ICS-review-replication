# *********************************************
# Authors: Aurore Archimbaud
#          TBS Business School
# 20/10/2026
# *********************************************
# Reproduce Figure 5

# Libraries --------------------------------------------------------------

# For ICS
library(ICSOutlier)
library(ICSClust)


# For PCA
library(rrcov)

# For plotting
library(ggplot2)
library(ggrepel)
library(patchwork)


# Data  --------------------------------------------------------------------
# HTP data: the observations 581 and 619 are considered as outliers
data(HTP)
outliers <- c(581, 619)
n <- nrow(HTP)
type <- factor(ifelse(seq_len(n) %in% outliers, "outlier", "normal"),
               levels = c("normal", "outlier"))

colors <- c(normal = "grey50", outlier = "black")

# Common theme for all panels
theme_paper <- theme_minimal(base_size = 16) +
  theme(strip.text = element_text(size = 15))


# Methods ------------------------------------------------------------------

# ICS COV-COV4
ics_cov_cov4 <- ICS(HTP)

# ICS reversed: COV4-COV (only the last component is kept)
ics_cov4_cov <- ICS(HTP, S1 = ICS_cov4, S2 = ICS_cov)

# Robust scatter: reweighted MCD with alpha = 0.75
MCD <- ICS_mcd_rwt(HTP, location = TRUE, alpha = 0.75)

# ROBPCA
ROBPCA <- PcaHubert(HTP)


# Panels -------------------------------------------------------------------
# Each panel carries its header as a facet strip (parsed, so that indices
# like COV[4] or ICSD[k==1]^2 are rendered as math).

# Panel 1: ICS screeplot (COV-COV4), first 5 and last 6 components
sel <- c(1:5, 83:88)
df_eig <- data.frame(
  component = factor(paste0("IC.", sel), levels = paste0("IC.", sel)),
  value     = ics_cov_cov4$gen_kurtosis[sel],
  panel     = '"ICS (COV-COV"[4]*"): screeplot"'
)

plot_eig <- ggplot(df_eig, aes(x = component, y = value)) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 5.5, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 6.5)) +
  facet_wrap(~ panel, labeller = label_parsed) +
  labs(x = "ICS Components", y = "Generalized kurtosis") +
  theme_paper +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panels 2-5: index plots of squared distances
# (mahalanobis() returns the squared distance, like ICSD^2)
make_index_plot <- function(z, header, ylab = "Squared distance") {
  df <- data.frame(ID = seq_len(n), Type = type, Z = z, panel = header)
  ggplot(df, aes(x = ID, y = Z, color = Type, shape = Type)) +
    geom_point(data = subset(df, Type == "normal"), alpha = 0.7, size = 2.5) +
    geom_point(data = subset(df, Type == "outlier"), size = 3) +
    geom_text_repel(data = subset(df, Type == "outlier"),
                    aes(label = ID), point.size = 2, size = 5,
                    max.overlaps = Inf, show.legend = FALSE) +
    facet_wrap(~ panel, labeller = label_parsed) +
    scale_color_manual(values = colors) +
    labs(x = "Observation number", y = ylab,
         color = NULL, shape = NULL) +
    theme_paper
}

plot_ics <- make_index_plot(ics_cov_cov4$scores[, 1]^2,
                            header = '"ICS (COV-COV"[4]*")"',
                            ylab   = expression(ICSD[k == 1]^2))

plot_ics_rev <- make_index_plot(ics_cov4_cov$scores[, 88]^2,
                                header = '"ICS (COV"[4]*"-COV)"',
                                ylab   = expression(ICSD[k == 88]^2))

plot_MD <- make_index_plot(mahalanobis(HTP, colMeans(HTP), cov(HTP)),
                           '"Mahalanobis distance"')

plot_RMD <- make_index_plot(mahalanobis(HTP, MCD$location, MCD$scatter),
                            '"Robust Mahalanobis distance (MCD"[0.75]*")"')

# Panel 6: ROBPCA diagnostic plot
df_PCA <- data.frame(ID = seq_len(n), Type = type,
                     SD = ROBPCA@sd, OD = ROBPCA@od,
                     panel = '"ROBPCA"')

plot_PCA <- ggplot(df_PCA, aes(x = SD, y = OD, color = Type, shape = Type)) +
  geom_point(data = subset(df_PCA, Type == "normal"), alpha = 0.7, size = 2.5) +
  geom_point(data = subset(df_PCA, Type == "outlier"), size = 3) +
  geom_vline(xintercept = ROBPCA@cutoff.sd, linetype = "dashed") +
  geom_hline(yintercept = ROBPCA@cutoff.od, linetype = "dashed") +
  geom_text_repel(data = subset(df_PCA, Type == "outlier"),
                  aes(label = ID), point.size = 2, size = 5,
                  max.overlaps = Inf, show.legend = FALSE) +
  facet_wrap(~ panel, labeller = label_parsed) +
  scale_color_manual(values = colors) +
  labs(x = "Score distance", y = "Orthogonal distance",
       color = NULL, shape = NULL) +
  theme_paper


# One figure: 3 rows x 2 columns, single shared legend ----------------------
fig_HTP <- plot_eig + plot_ics +
  plot_ics_rev + plot_MD +
  plot_RMD + plot_PCA +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "top")
#fig_HTP

ggsave("figures/HTP_outliers.pdf", plot = fig_HTP, width = 12, height = 14)
