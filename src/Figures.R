library(dplyr)
library(ggplot2)
library(ggpubr)
library(plotly)
library(forestplot)
library(VennDiagram)
source("src/functions.R")
pix <- 0.393700787


# Forest plot -------------------------------------------------------------

pheno <- c("Outcomes", "Heart_Failure", "Cardiomyopathy", 
           "Ventricular_arrhythmias", "Atrial_arrhythmias", "Heart_Arrhythmia", 
           "Chronic_ischaemic_heart_disease", "Angina", 
           "Cardiovascular_Death", "All_cause_mortality")
dat <- read.delim("results/output/CVD_enrichment_fisher.tsv") %>%
  filter(CM %in% c("ACM", "DCM", "HCM")) %>% filter(Phenotype %in% pheno) %>% arrange(factor(Phenotype, levels = pheno)) 
dat$CM <- paste0(dat$CM, " G+")

dat[, 2:4] <- lapply(dat[, 2:4], as.numeric)
text <- data.frame()
for (i in pheno) {

  if (i == "Outcomes") {
    t <- c(i, "OR (95% CI)")
  } else {
    temp <- dat %>% filter(Phenotype == i)
    t <- c(i, paste0(format(round(temp[1, 2], 2), nsmall = 2), " (", format(round(temp[1, 3], 2), nsmall = 2), ";", format(round(temp[1, 4], 2), nsmall = 2), ")\n",
                     format(round(temp[2, 2], 2), nsmall = 2), " (", format(round(temp[2, 3], 2), nsmall = 2), ";", format(round(temp[2, 4], 2), nsmall = 2), ")\n",
                     format(round(temp[3, 2], 2), nsmall = 2), " (", format(round(temp[3, 3], 2), nsmall = 2), ";", format(round(temp[3, 4], 2), nsmall = 2), ")"))
  }
  text <- rbind(text, t)
  text[1:2] <- lapply(text[1:2], as.character)
}
text[, 1] <- gsub("_", " ", text[, 1])
text[,1][text[, 1] == "Heart Failure"] <- "Heart failure"
text[,1][text[, 1] == "Heart Arrhythmia"] <- "Heart arrhythmia"
text[,1][text[, 1] == "Cardiovascular Death"] <- "Cardiovascular death"

svglite("results/figures/Figure3-Forest_plot.svg", width = 10 * pix, height = 8 * pix)
dat %>% group_by(CM) %>% 
  forestplot(labeltext = text, mean = OR, lower = LCI, upper = UCI, zero = 1, 
             fn.ci_norm = c(fpDrawNormalCI, fpDrawDiamondCI, fpDrawCircleCI),
             shapes_gp = fpShapesGp(box = c("#ffd167", "#168ab2", "#ef476f") %>% 
                                      lapply(function(x) gpar(fill = x, col = x)),
                                    line = c("#ffd167", "#168ab2", "#ef476f") %>% 
                                      lapply(function(x) gpar(fill = x, col = x)),
                                    default = gpar(vertices = TRUE)),
             boxsize = .1, graph.pos = 2, align = "r", lwd.xaxis = .5,
             lwd.zero = .5, lwd.ci = .5, xlog = TRUE,
             txt_gp = fpTxtGp(label = gpar(cex = .3), xlab = gpar(cex = .35),
                              ticks = gpar(cex= .3), legend = gpar(cex = .3)),
             is.summary = c(TRUE, rep(FALSE, 18)), xlab = "Odds ratio (95% CI)")
dev.off()


# Venn Diagram phenotypes -------------------------------------------------

df <- read.delim("data/processed/MCM_clean_final.tsv")
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
                     cat.cex = rep(.3, 4), cex = rep(.3, 15), 
                     alpha = rep(0.8, 4), lwd = rep(.2, 4),
                     fontfamily = rep("Helvetica", 15), 
                     cat.fontfamily = rep("Helvetica", 4), 
                     category = c("HF", "CM", "Ischemia", "VA"))
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
                     cat.cex = rep(.3, 4), cex = rep(.3, 15), 
                     alpha = rep(0.8, 4), lwd = rep(.2, 4),
                     fontfamily = rep("Helvetica", 15), 
                     cat.fontfamily = rep("Helvetica", 4), 
                     category = c("HF", "CM", "Ischemia", "VA"))
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
                     cat.cex = rep(.3, 4), cex = rep(.3, 15), 
                     alpha = rep(0.8, 4), lwd = rep(.2, 4),
                     fontfamily = rep("Helvetica", 15), 
                     cat.fontfamily = rep("Helvetica", 4), 
                     category = c("HF", "CM", "Ischemia", "VA"))
# ggsave("results/figures/HCM_venn.svg", v3)

ggarrange(v1, v2, v3, labels = c("A) ACM G+", "B) DCM G+", "C) HCM G+"), ncol = 3,
          font.label = list(size = 4, color = "black", face = "bold", 
                            family = "Helvetica"))
ggsave("results/figures/Figure4-Venn.svg", width = 14 * pix, height = 4 * pix)
rm(v1, v2, v3, acm, dcm, hcm)

# Boxplots ----------------------------------------------------------------

# Create dataframe with desired variables
ex <- c("LV", "RV")
cmr <- df %>% select(starts_with("RV"), starts_with("LV"), contains("Wall_thickness"), 
                     contains("Ecc", ignore.case = FALSE), 
                     contains("Ell", ignore.case = FALSE)) %>% 
  select(-any_of(ex)) %>% names()
new <- df %>% select(f.eid, BSA, CM, Pheno, Sex, any_of(cmr)) %>% filter(Pheno == "Non-Diagnosed")

wt <- new %>% select(f.eid, contains("segment"))
wt$Max_WT <- NA
for (i in 1:nrow(wt)) {
  wt[i, "Max_WT"] <- max(wt[i, 2:(ncol(wt)-1)])
  wt$Max_WT <- as.numeric(wt$Max_WT)
}
wt <- wt %>% select(f.eid, Max_WT)
new <- merge(new, unique(wt))
new$Max_WT[new$Max_WT == 0] <- NA
rm(wt, cmr, ex, i)

# Make plots
# ggplot(new, aes(x = Max_WT, fill = Sex)) +
#   geom_histogram(position = "dodge") + 
#   labs(x = "Maximum LV wall thickness (mm)", y = "Count",
#        title = "Distribution of Maximum Left Ventricular wall thickness") +
#   scale_fill_manual(values = c("tomato4", "cornflowerblue")) +
#   my_theme() 
# ggsave("results/figures/Max_wall_thickness.svg", width = 6 * pix, height = 5  * pix)

 p1 <- ggplot(new, 
             aes(x = factor(CM, levels = c("Controls", "ACM", "DCM", "HCM"),
                            labels = c("Controls G-P-", "ACM G+P-", "DCM G+P-", 
                                       "HCM G+P-")), 
                 y = LVEF, fill = CM)) +
  geom_boxplot(outlier.shape = 18, outlier.size = .4, alpha = 0.9, 
               lwd = .2, fatten = .8) +
  labs(x = "Cardiomyopathy", y = "LVEF (%)") +
  annotate("segment", x = c(1, 1, 3), xend = c(3, 1, 3), 
           y = c(90, 88, 88), yend = c(90, 90, 90), lwd = .2) +
  annotate("text", x = 2, y = 92, label = "p=0.009", size = 1) +
  scale_fill_manual(breaks = c("Controls", "ACM", "DCM", "HCM"),
                    values = c("#4cbd97", "#ffd167", "#168ab2", "#ef476f")) +
  my_theme() +
  theme(legend.position="none", text = element_text(size = 3),
        line = element_line(size = .2))

p1.1 <- ggplot(new, 
               aes(x = factor(CM, levels = c("Controls", "ACM", "DCM", "HCM"),
                              labels = c("Controls G-P-", "ACM G+P-", "DCM G+P-", 
                                         "HCM G+P-")), 
                   y = RVEF, fill = CM)) +
  geom_boxplot(outlier.shape = 18, outlier.size = .4, alpha = 0.9, 
               lwd = .2, fatten = .8) +
  labs(x = "Cardiomyopathy", y = "RVEF (%)") +
  annotate("segment", x = c(1, 1, 4), xend = c(4, 1, 4), 
           y = c(92, 90, 90), yend = c(92, 92, 92), lwd = .2) +
  annotate("text", x = 2.5, y = 94, label = "p=0.034", size = 1) +
  scale_fill_manual(breaks = c("Controls", "ACM", "DCM", "HCM"),
                    values = c("#4cbd97", "#ffd167", "#168ab2", "#ef476f")) +
  my_theme() +
  theme(legend.position="none", text = element_text(size = 3),
        line = element_line(size = .2))

p2 <- ggplot(new, 
             aes(x = factor(CM, levels = c("Controls", "ACM", "DCM", "HCM"),
                            labels = c("Controls G-P-", "ACM G+P-", "DCM G+P-", 
                                       "HCM G+P-")), 
                 y = LVEDVi, fill = CM)) +
  geom_boxplot(outlier.shape = 18, outlier.size = .4, alpha = 0.9, 
               lwd = .2, fatten = .8) +
  labs(x = "Cardiomyopathy", y = "LVEDVi (ml/m2)") +
  scale_fill_manual(breaks = c("Controls", "ACM", "DCM", "HCM"),
                    values = c("#4cbd97", "#ffd167", "#168ab2", "#ef476f")) +
  my_theme() +
  theme(legend.position="none", text = element_text(size = 3),
        line = element_line(size = .2))

p2.1 <- ggplot(new, 
               aes(x = factor(CM, levels = c("Controls", "ACM", "DCM", "HCM"),
                              labels = c("Controls G-P-", "ACM G+P-", "DCM G+P-", 
                                         "HCM G+P-")), 
                   y = RVEDVi, fill = CM)) +
  geom_boxplot(outlier.shape = 18, outlier.size = .4, alpha = 0.9, 
               lwd = .2, fatten = .8) +
  labs(x = "Cardiomyopathy", y = "RVEDVi (ml/m2)")  +
  annotate("segment", x = c(1, 1, 3), xend = c(3, 1, 3), 
           y = c(166, 160, 160), yend = c(166, 166, 166), lwd = .2) +
  annotate("text", x = 2, y = 170, label = "p=0.048", size = 1) +
  scale_fill_manual(breaks = c("Controls", "ACM", "DCM", "HCM"),
                    values = c("#4cbd97", "#ffd167", "#168ab2", "#ef476f")) +
  my_theme() +
  theme(legend.position="none", text = element_text(size = 3),
        line = element_line(size = .2))

p3 <- ggplot(new, 
             aes(x = factor(CM, levels = c("Controls", "ACM", "DCM", "HCM"),
                            labels = c("Controls G-P-", "ACM G+P-", "DCM G+P-", 
                                       "HCM G+P-")), 
                 y = Max_WT, fill = CM)) +
  geom_boxplot(outlier.shape = 18, outlier.size = .4, alpha = 0.9, 
               lwd = .2, fatten = .8) +
  labs(x = "Cardiomyopathy", y = "Maximum wall thickness (mm)") +
  scale_fill_manual(breaks = c("Controls", "ACM", "DCM", "HCM"),
                    values = c("#4cbd97", "#ffd167", "#168ab2", "#ef476f")) +
  my_theme() +
  theme(legend.position="none", text = element_text(size = 3),
        line = element_line(size = .2))

p4 <- ggplot(new, 
             aes(x = factor(CM, levels = c("Controls", "ACM", "DCM", "HCM"),
                            labels = c("Controls G-P-", "ACM G+P-", "DCM G+P-", 
                                       "HCM G+P-")), 
                 y = peakEll4Ch, fill = CM)) +
  geom_boxplot(outlier.shape = 18, outlier.size = .4, alpha = 0.9, 
               lwd = .2, fatten = .8) +
  labs(x = "Cardiomyopathy", y = "Peak longitudinal strain (%)") +
  annotate("segment", x = c(1, 1, 3), xend = c(3, 1, 3), 
           y = c(-5, -7, -7), yend = c(-5, -5, -5), lwd = .2) +
  annotate("text", x = 2, y = -4, label = "p=0.009", size = 1) +
  scale_fill_manual(breaks = c("Controls", "ACM", "DCM", "HCM"),
                    values = c("#4cbd97", "#ffd167", "#168ab2", "#ef476f")) +
  my_theme() +
  theme(legend.position="none", text = element_text(size = 3),
        line = element_line(size = .2))

ggarrange(p1, p1.1, p2, p2.1, p3, p4, labels = "AUTO", 
          font.label = list(size = 3.5, color = "black", face = "bold", 
                            family = "Helvetica"),
          ncol = 2, nrow = 3)
ggsave("results/figures/Figure6-CMR_box.svg", width = 8 * pix, height = 12 * pix)
rm(p1, p1.1, p2, p2.1, p3, p4, new)
