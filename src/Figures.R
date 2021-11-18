library(dplyr)
library(ggplot2)
library(ggpubr)
library(plotly)
library(forestplot)


# Forest plot -------------------------------------------------------------

dat <- read.delim("results/output/ForestData.txt")

dat[, 2:6] <- lapply(dat[, 2:6], as.numeric)
text <- data.frame(Phenotype = unique(dat$Phenotype))
text$OR <- NA
for (i in seq(1, nrow(dat), 3)) {
  if (i == 1) {
    j = 1
  } else {
    j = (i + 2) / 3
  }
  text[j, 2] <- paste0(round(dat[i, 4], 2), " (", round(dat[i, 5], 2), ";", round(dat[i, 6], 2), ")\n",
                       round(dat[i+1, 4], 2), " (", round(dat[i+1, 5], 2), ";", round(dat[i+1, 6], 2), ")\n",
                       round(dat[i+2, 4], 2), " (", round(dat[i+2, 5], 2), ";", round(dat[i+2, 6], 2), ")")
  
}
text[1, 2] <- "OR (95% CI)"

pdf("results/figures/Figure3-Forest_plot.pdf", width = 10, height = 10)
dat %>% group_by(CM) %>%
  forestplot(labeltext = text, 
             mean = OR, lower = LCI, upper = UCI, xlog = TRUE, 
             fn.ci_norm = c(fpDrawNormalCI, fpDrawDiamondCI, fpDrawCircleCI),
             shapes_gp = fpShapesGp(box = c("#ffd167", "#168ab2", "#ef476f") %>% 
                                      lapply(function(x) gpar(fill = x, col = "#555555")),
                                    default = gpar(vertices = TRUE)),
             boxsize = .1, graph.pos = 2, align = "r",
             is.summary = c(TRUE, rep(FALSE, 18)), xlab = "Odds ratio",
             col = fpColors(box = "black", line = "darkgrey"))
dev.off()

# Venn Diagram phenotypes -------------------------------------------------

acm <- df %>% filter(CM == "ACM")
grid.newpage()
v1 <- draw.quad.venn(area1 = nrow(subset(acm, Heart_Failure_sum == "Yes")), 
                     area2 = nrow(subset(acm, Cardiomyopathy_sum == "Yes")), 
                     area3 = nrow(subset(acm, Chronic_ischaemic_heart_disease_sum == "Yes")), 
                     area4 = nrow(subset(acm, Ventricular_arrhythmias_sum == "Yes")),
                     n12 = nrow(subset(acm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes")), 
                     n13 = nrow(subset(acm, Heart_Failure_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes")), 
                     n14 = nrow(subset(acm, Heart_Failure_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")), 
                     n23 = nrow(subset(acm, Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes")), 
                     n24 = nrow(subset(acm, Cardiomyopathy_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")), 
                     n34 = nrow(subset(acm, Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")),
                     n123 = nrow(subset(acm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes")), 
                     n124 = nrow(subset(acm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")), 
                     n134 = nrow(subset(acm, Heart_Failure_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")), 
                     n234 = nrow(subset(acm, Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")), 
                     n1234 = nrow(subset(acm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")),
                     fill = c("#ffd167", "#4cbd97", "#168ab2", "#ef476f"), 
                     cat.cex = rep(3, 4), cex = rep(3, 15), alpha = rep(0.8, 4),
                     fontfamily = rep("Helvetica", 15), 
                     cat.fontfamily = rep("Helvetica", 4), 
                     category = c("HF", "CM", "Ischaemia", "VA"))
# ggsave("results/figures/ACM_venn.svg", v1)

dcm <- df %>% filter(CM == "DCM")
grid.newpage()
v2 <- draw.quad.venn(area1 = nrow(subset(dcm, Heart_Failure_sum == "Yes")), 
                     area2 = nrow(subset(dcm, Cardiomyopathy_sum == "Yes")), 
                     area3 = nrow(subset(dcm, Chronic_ischaemic_heart_disease_sum == "Yes")), 
                     area4 = nrow(subset(dcm, Ventricular_arrhythmias_sum == "Yes")),
                     n12 = nrow(subset(dcm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes")), 
                     n13 = nrow(subset(dcm, Heart_Failure_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes")), 
                     n14 = nrow(subset(dcm, Heart_Failure_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")), 
                     n23 = nrow(subset(dcm, Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes")), 
                     n24 = nrow(subset(dcm, Cardiomyopathy_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")), 
                     n34 = nrow(subset(dcm, Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")),
                     n123 = nrow(subset(dcm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes")), 
                     n124 = nrow(subset(dcm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")), 
                     n134 = nrow(subset(dcm, Heart_Failure_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")), 
                     n234 = nrow(subset(dcm, Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")), 
                     n1234 = nrow(subset(dcm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")),
                     fill = c("#ffd167", "#4cbd97", "#168ab2", "#ef476f"), 
                     cat.cex = rep(3, 4), cex = rep(3, 15), alpha = rep(0.8, 4),
                     fontfamily = rep("Helvetica", 15), 
                     cat.fontfamily = rep("Helvetica", 4), 
                     category = c("HF", "CM", "Ischaemia", "VA"))
# ggsave("results/figures/DCM_venn.svg", v2)

hcm <- df %>% filter(CM == "HCM")
grid.newpage()
v3 <- draw.quad.venn(area1 = nrow(subset(hcm, Heart_Failure_sum == "Yes")), 
                     area2 = nrow(subset(hcm, Cardiomyopathy_sum == "Yes")), 
                     area3 = nrow(subset(hcm, Chronic_ischaemic_heart_disease_sum == "Yes")), 
                     area4 = nrow(subset(hcm, Ventricular_arrhythmias_sum == "Yes")),
                     n12 = nrow(subset(hcm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes")), 
                     n13 = nrow(subset(hcm, Heart_Failure_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes")), 
                     n14 = nrow(subset(hcm, Heart_Failure_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")), 
                     n23 = nrow(subset(hcm, Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes")), 
                     n24 = nrow(subset(hcm, Cardiomyopathy_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")), 
                     n34 = nrow(subset(hcm, Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")),
                     n123 = nrow(subset(hcm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes")), 
                     n124 = nrow(subset(hcm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")), 
                     n134 = nrow(subset(hcm, Heart_Failure_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")), 
                     n234 = nrow(subset(hcm, Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")), 
                     n1234 = nrow(subset(hcm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")),
                     fill = c("#ffd167", "#4cbd97", "#168ab2", "#ef476f"), 
                     cat.cex = rep(3, 4), cex = rep(3, 15), alpha = rep(0.8, 4),
                     fontfamily = rep("Helvetica", 15), 
                     cat.fontfamily = rep("Helvetica", 4), 
                     category = c("HF", "CM", "Ischaemia", "VA"))
# ggsave("results/figures/HCM_venn.svg", v3)

ggarrange(v1, v2, v3, labels = c("A) ACM", "B) DCM", "C) HCM"), ncol = 3,
          font.label = list(size = 35, color = "black", face = "bold", family = "Helvetica"))
ggsave("results/figures/Figure4-Venn.pdf", width = 30, height = 8)


# Boxplots ----------------------------------------------------------------

# Create dataframe with desired variables
ex <- c("LV", "RV")
cmr <- df %>% select(starts_with("RV"), starts_with("LV"), contains("Wall_thickness"), 
                     contains("Ecc", ignore.case = FALSE), 
                     contains("Ell", ignore.case = FALSE)) %>% 
  select(-any_of(ex)) %>% names()
new <- df %>% select(f.eid, BSA, CM, Pheno, any_of(cmr)) %>% filter(Pheno == "Non-Diagnosed")

wt <- new %>% select(f.eid, contains("segment"))
wt$Max_WT <- NA
for (i in 1:nrow(wt)) {
  wt[i, "Max_WT"] <- max(wt[i, 2:(ncol(wt)-1)])
  wt$Max_WT <- as.numeric(wt$Max_WT)
}
wt <- wt %>% select(f.eid, Max_WT)
new <- merge(new, wt)
new$Max_WT[new$Max_WT == 0] <- NA

# Make plots
p1 <- ggplot(new, aes(x = factor(CM, levels = c("Controls", "ACM", "DCM", "HCM")), 
                      y = LVEF)) +
  geom_boxplot(outlier.shape = 18) +
  labs(x = "Cardiomyopathy", y = "LVEF (%)") +
  annotate("segment", x = c(1, 1, 3), xend = c(3, 1, 3), y = c(90, 88, 88), yend = c(90, 90, 90)) +
  annotate("text", x = 2, y = 92, label = "*") +
  my_theme()

p1.1 <- ggplot(new, aes(x = factor(CM, levels = c("Controls", "ACM", "DCM", "HCM")), 
                        y = RVEF)) +
  geom_boxplot(outlier.shape = 18) +
  labs(x = "Cardiomyopathy", y = "RVEF (%)") +
  annotate("segment", x = c(1, 1, 4), xend = c(4, 1, 4), y = c(92, 90, 90), yend = c(92, 92, 92)) +
  annotate("text", x = 2.5, y = 94, label = "**") +
  my_theme()

p2 <- ggplot(new, aes(x = factor(CM, levels = c("Controls", "ACM", "DCM", "HCM")), 
                      y = LVEDVi)) +
  geom_boxplot(outlier.shape = 18) +
  labs(x = "Cardiomyopathy", y = "LVEDVi (ml/m2)") +
  my_theme()

p2.1 <- ggplot(new, aes(x = factor(CM, levels = c("Controls", "ACM", "DCM", "HCM")), 
                        y = RVEDVi)) +
  geom_boxplot(outlier.shape = 18) +
  labs(x = "Cardiomyopathy", y = "RVEDVi (mm)") +
  my_theme()

p3 <- ggplot(new, aes(x = factor(CM, levels = c("Controls", "ACM", "DCM", "HCM")), 
                      y = Max_WT), fill = factor(CM, levels = c("Controls", "ACM", "DCM", "HCM"))) +
  geom_boxplot(outlier.shape = 18) +
  labs(x = "Cardiomyopathy", y = "Maximum wall thickness (mm)") +
  scale_fill_manual(breaks = c("Controls", "ACM", "DCM", "HCM"),
                    values = c("green", "yellow", "blue", "red")) +
  my_theme()

p4 <- ggplot(new, 
             aes(x = factor(CM, levels = c("Controls", "ACM", "DCM", "HCM")), 
                 y = peakEll4Ch)) +
  geom_boxplot(outlier.shape = 18) +
  labs(x = "Cardiomyopathy", y = "Peak longitudinal strain (%)") +
  annotate("segment", x = c(1, 1, 3), xend = c(3, 1, 3), y = c(-5, -7, -7), yend = c(-5, -5, -5)) +
  annotate("text", x = 2, y = -4, label = "**") +
  my_theme()

ggarrange(p1, p1.1, p2, p2.1, p3, p4, 
          common.legend = TRUE, legend = "bottom",
          labels = c("A", "B", "C", "D", "E", "F"), 
          ncol = 2, nrow = 3)
ggsave("Figure6-CMR_box.pdf", width = 20, height = 30)
