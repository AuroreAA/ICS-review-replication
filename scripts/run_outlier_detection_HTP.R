# *********************************************
# Authors: Aurore Archimbaud 
#          TBS Business School
# *********************************************


# Libraries --------------------------------------------------------------

# For ICS
library(ICSOutlier)
library(ICSClust)
library(ICSOutlier)

# For PCA
library(rrcov)

# For plotting
library(ggplot2)
library(ggrepel)

# For pre processing
library(magrittr)


# Data  --------------------------------------------------------------------
# HTP data: the observations 581 and 619 are considered as outliers
data(HTP)
outliers <- c(581, 619)
boxplot(HTP)
clusters <- rep("normal", nrow(HTP))
clusters[outliers] <- "outlier"

# ICS COV-COV4 ---------------------------------------------------------------------
ics_res <- ICS(HTP)

plot_ics_eigenvalues <- select_plot(ics_res, select = c(1:5, 83:88),
                                    type = "lines")
plot_ics_eigenvalues <- plot_ics_eigenvalues + 
  geom_vline(aes(xintercept = 5.5), linetype = "dashed") +
  ylim(0,6.5)+
  theme_minimal() +
  labs(x = "ICS-Components", fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position="top") 

pdf("figures/HTP_outliers_ics_eigenvalues.pdf",  width = 6.5, height = 5.75)
print(plot_ics_eigenvalues)
dev.off()


# screeplot: only one component to keep
df_scores <- data.frame(ID = 1:nrow(ics_res$scores),
                        Type = clusters,
                        Z = ics_res$scores[,1]^2)


# Create the y-axis label based on the combination of scatter matrices and
# the number of selected components

k_details <- 1
y_axis <- bquote("COV-COV"[4] : ICSD[ ~ 'k' == .(k_details) ] ^2 )

colors <-  c("grey50", "black")
plot_ics <- df_scores  %>%
  ggplot(aes(x = ID, y =  Z)) +
  geom_point(data = df_scores[df_scores$Type != "outlier",],
             aes(color = Type, shape = Type), alpha = 0.7, size = 2.5) +
  geom_point(data = df_scores[df_scores$Type == "outlier",],
             aes(color = Type, shape = Type), alpha = 1, size = 3) +
  scale_color_manual(values = colors) +
  ylim(0,700)+
  theme_minimal() +
  labs(x = "Observation Number", y = y_axis, fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position="top") +
  ggrepel::geom_text_repel(data = df_scores[df_scores$Type == "outlier",], 
                           aes(label = ID), point.size = 2, size = 5, 
                           max.overlaps = Inf) 
plot_ics


# save plot to pdf

pdf("figures/HTP_outliers_ics.pdf", width = 6.5, height = 5.75)
print(plot_ics)
dev.off()



# ICS reversed: COV4-COV ---------------------------------------------------------------------
# Run ICS with reversing order or the scatters: 
ics_res <- ICS(HTP, S1 = ICS_cov4, S2 = ICS_cov)

# screeplot: only one component to keep (the last one)
select_plot(ics_res)

df_scores <- data.frame(ID = 1:nrow(ics_res$scores),
                        Type = clusters,
                        Z = ics_res$scores[,88]^2)


# Create the y-axis label based on the combination of scatter matrices and
# the number of selected components

k_details <- 88
y_axis <- bquote(COV[4]*"- COV : ICSD"[ ~ 'k' == .(k_details) ] ^2 )


plot_ics <- df_scores  %>%
  ggplot(aes(x = ID, y =  Z)) +
  geom_point(data = df_scores[df_scores$Type != "outlier",],
             aes(color = Type, shape = Type), alpha = 0.7, size = 2.5) +
  geom_point(data = df_scores[df_scores$Type == "outlier",],
             aes(color = Type, shape = Type), alpha = 1, size = 3) +
  scale_color_manual(values = colors) +
  ylim(0,700)+
  theme_minimal() +
  labs(x = "Observation Number", y = y_axis, fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position="top") +
  ggrepel::geom_text_repel(data = df_scores[df_scores$Type == "outlier",], 
                           aes(label = ID), point.size = 2, size = 5,
                           max.overlaps = Inf) 
plot_ics


# save plot to pdf

pdf("figures/HTP_outliers_ics_rev.pdf", width = 6.5, height = 5.75)
print(plot_ics)
dev.off()





# Mahalanobs Distance ----------------------------------------------------------------------
df_MD <- data.frame(ID = 1:nrow(ics_res$scores),
                        Type = clusters,
                        Z = mahalanobis(HTP, colMeans(HTP), cov(HTP)))

plot_MD <- df_MD  %>%
  ggplot(aes(x = ID, y =  Z)) +
  geom_point(data =  df_MD[df_scores$Type != "outlier",],
             aes(color = Type, shape = Type), alpha = 0.7, size = 2.5) +
  geom_point(data =  df_MD[df_scores$Type == "outlier",],
             aes(color = Type, shape = Type), alpha = 1, size = 3) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Observation Number", y = "Mahalanobis distance", fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position="top") +
  ggrepel::geom_text_repel(data = df_MD[df_MD$Type == "outlier",], 
                           aes(label = ID), point.size = 2, size = 5,
                           max.overlaps = Inf) 
plot_MD


# save plot to pdf

pdf("figures/HTP_outliers_MD.pdf", width = 6.5, height = 5.75)
print(plot_MD)
dev.off()




# Robust Mahalanobis distance with MCD 0.75 ----------------------------------------------------------------------
MCD  <- ICS_mcd_rwt(HTP, location = TRUE, alpha = 0.75)
df_RMD <- data.frame(ID = 1:nrow(ics_res$scores),
                     Type = clusters,
                     Z = mahalanobis(HTP, MCD$location, MCD$scatter))

plot_RMD <- df_RMD  %>%
  ggplot(aes(x = ID, y =  Z)) +
  geom_point(data =  df_RMD[df_scores$Type != "outlier",],
             aes(color = Type, shape = Type), alpha = 0.7, size = 2.5) +
  geom_point(data =  df_RMD[df_scores$Type == "outlier",],
             aes(color = Type, shape = Type), alpha = 1, size = 3) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Observation Number", y = "Robust Mahalanobis distance", fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position="top") +
  ggrepel::geom_text_repel(data = df_RMD[df_RMD$Type == "outlier",], 
                           aes(label = ID), point.size = 2, size = 5,
                           max.overlaps = Inf) 
plot_RMD


# save plot to pdf

pdf("figures/HTP_outliers_RMD75.pdf", width = 6.5, height = 5.75)
print(plot_RMD)
dev.off()




# ROBPCA ---------------------------------------------------------------------
ROBPCA  <- PcaHubert(HTP)
plot(ROBPCA)
df_PCA <- data.frame(ID = 1:nrow(ics_res$scores),
                     Type = clusters,
                     SD = ROBPCA@sd,
                     OD = ROBPCA@od)


plot_PCA <- df_PCA  %>%
  ggplot(aes(x = SD, y =  OD)) +
  geom_point(data =   df_PCA[df_scores$Type != "outlier",],
             aes(color = Type, shape = Type), alpha = 0.7, size = 2.5) +
  geom_point(data =   df_PCA[df_scores$Type == "outlier",],
             aes(color = Type, shape = Type), alpha = 1, size = 3) +
  geom_vline(aes(xintercept = ROBPCA@cutoff.sd)) +
  geom_hline(aes(yintercept = ROBPCA@cutoff.od)) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "Score distance", y = "Orthogonal distance", fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        legend.position="top") +
  ggrepel::geom_text_repel(data = df_PCA[df_PCA$Type == "outlier",], 
                           aes(label = ID), point.size = 2, size = 5,
                           max.overlaps = Inf) 
plot_PCA


# save plot to pdf

pdf("figures/HTP_outliers_ROBPCA.pdf", width = 6.5, height = 5.75)
print(plot_PCA)
dev.off()

