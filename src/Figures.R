library(dplyr)
library(ggplot2)
library(ggpubr)
library(plotly)
library(forestplot)
library(VennDiagram)
library(svglite)
source("src/functions.R")
pix <- 0.393700787


# Forest plot -------------------------------------------------------------

pheno <- c("Outcomes", "Heart_Failure", "Cardiomyopathy", "HFCM", "Pheno",
           "Ventricular_arrhythmias", "Atrial_arrhythmias", "Heart_Arrhythmia", 
           "Chronic_ischaemic_heart_disease", "Angina", 
           "Cardiovascular_Death", "All_cause_mortality")
dat <- read.delim("~/surfdrive/Mendelian_CM_WES_UKB/Difference_testing_v4.tsv") %>%
  filter(CM %in% c("ACM", "DCM", "HCM", "strict HCM")) %>%
  filter(Phenotype %in% pheno) %>% arrange(factor(Phenotype, levels = pheno))
dat$CM[dat$CM == "ACM"] <- "ARVC"
dat$CM <- paste0(dat$CM, " G+")
dat$CM[dat$CM == "strict HCM G+"] <- "HCM* G+"

dat[, 2:4] <- lapply(dat[, 2:4], as.numeric)
dat <- dat[order(dat[, 1], dat[, 7]),] %>% arrange(factor(Phenotype, levels = pheno)) 

text <- data.frame()
for (i in pheno) {
  
  if (i == "Outcomes") {
    t <- c(i, "OR (95% CI)")
  } else {
    temp <- dat %>% filter(Phenotype == i)
    for (g in 1:nrow(temp)) {
      # Paste together the OR and 95% CI into text
      s <- paste0(format(round(temp[g, 2], 2), nsmall = 2), " (", format(round(temp[g, 3], 2), nsmall = 2), ";", format(round(temp[g, 4], 2), nsmall = 2), ")")
      if (g == 1) {
        # Create variable 
        p <- s
      } else {
        # Paste with newline to variable p
        p <- paste0(p, "\n", s)
      } # end if first gene
      t <- c(i, p)
    } # iterate over CMs  
  }
  text <- rbind(text, t)
  text[1:2] <- lapply(text[1:2], as.character)
}
text[, 1] <- gsub("_", " ", text[, 1])
text[,1][text[, 1] == "Heart Failure"] <- "Heart failure"
text[,1][text[, 1] == "HFCM"] <- "Heart failure + cardiomyopathy"
text[,1][text[, 1] == "Pheno"] <- "Phenotype positive"
text[,1][text[, 1] == "Angina"] <- "Angina pectoris"
text[,1][text[, 1] == "Chronic ischaemic heart disease"] <- "Chronic ischemic heart disease"
text[,1][text[, 1] == "Heart Arrhythmia"] <- "Self-reported heart arrhythmias"
text[,1][text[, 1] == "Cardiovascular Death"] <- "Cardiovascular death"
text[,1][text[, 1] == "All cause mortality"] <- "All-cause mortality"

pdf("~/surfdrive/Mendelian_CM_WES_UKB/FigureS1.pdf", paper = "a4",
    width = 10 * pix, height = 12 * pix)
dat %>% group_by(CM) %>% 
  forestplot(labeltext = text, mean = Estimate, lower = LCI, upper = UCI, zero = 1, 
             fn.ci_norm = c(fpDrawNormalCI, fpDrawDiamondCI, fpDrawCircleCI, fpDrawCircleCI),
             shapes_gp = fpShapesGp(box = c("#ffd167", "#168ab2", "#ef476f", "black") %>% 
                                      lapply(function(x) gpar(fill = x, col = x)),
                                    line = c("#ffd167", "#168ab2", "#ef476f", "black") %>% 
                                      lapply(function(x) gpar(fill = x, col = x)),
                                    default = gpar(vertices = TRUE)),
             boxsize = .1, graph.pos = 2, align = "r", lwd.xaxis = .5,
             lwd.zero = .5, lwd.ci = .5, xlog = TRUE,
             txt_gp = fpTxtGp(label = gpar(cex = .3), xlab = gpar(cex = .35),
                              ticks = gpar(cex= .3), legend = gpar(cex = .3)),
             is.summary = c(TRUE, rep(FALSE, 24)), xlab = "Odds ratio (95% CI)")
dev.off()

dat <- dat %>% filter(CM != "HCM* G+")
text <- data.frame()
for (i in pheno) {
  
  if (i == "Outcomes") {
    t <- c(i, "OR (95% CI)")
  } else {
    temp <- dat %>% filter(Phenotype == i)
    for (g in 1:nrow(temp)) {
      # Paste together the OR and 95% CI into text
      s <- paste0(format(round(temp[g, 2], 2), nsmall = 2), " (", format(round(temp[g, 3], 2), nsmall = 2), ";", format(round(temp[g, 4], 2), nsmall = 2), ")")
      if (g == 1) {
        # Create variable 
        p <- s
      } else {
        # Paste with newline to variable p
        p <- paste0(p, "\n", s)
      } # end if first gene
      t <- c(i, p)
    } # iterate over CMs  
  }
  text <- rbind(text, t)
  text[1:2] <- lapply(text[1:2], as.character)
}
text[, 1] <- gsub("_", " ", text[, 1])
text[,1][text[, 1] == "Heart Failure"] <- "Heart failure"
text[,1][text[, 1] == "HFCM"] <- "Heart failure + cardiomyopathy"
text[,1][text[, 1] == "Pheno"] <- "Phenotype positive"
text[,1][text[, 1] == "Angina"] <- "Angina pectoris"
text[,1][text[, 1] == "Chronic ischaemic heart disease"] <- "Chronic ischemic heart disease"
text[,1][text[, 1] == "Heart Arrhythmia"] <- "Self-reported heart arrhythmias"
text[,1][text[, 1] == "Cardiovascular Death"] <- "Cardiovascular death"
text[,1][text[, 1] == "All cause mortality"] <- "All-cause mortality"

pdf("~/surfdrive/Mendelian_CM_WES_UKB/Figure4.pdf", paper = "a4",
    width = 10 * pix, height = 9 * pix)
dat %>% group_by(CM) %>% 
  forestplot(labeltext = text, mean = Estimate, lower = LCI, upper = UCI, zero = 1, 
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
rm(dat, temp, text, i, pheno, t)


# Venn Diagram phenotypes -------------------------------------------------

df <- read.delim("data/processed/MCM_clean_final.tsv")
acm <- df %>% filter(CM == "ACM")
grid.newpage()
v1 <- draw.quad.venn(area1 = round(nrow(subset(acm, Heart_Failure_sum == "Yes")) / nrow(acm) * 100, 1),
                     area2 = round(nrow(subset(acm, Cardiomyopathy_sum == "Yes")) / nrow(acm) * 100, 1),
                     area3 = round(nrow(subset(acm, Chronic_ischaemic_heart_disease_sum == "Yes")) / nrow(acm) * 100,  1),
                     area4 = round(nrow(subset(acm, Ventricular_arrhythmias_sum == "Yes")) / nrow(acm) * 100, 1),
                     n12 = round(nrow(subset(acm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes")) / nrow(acm) * 100,  1),
                     n13 = round(nrow(subset(acm, Heart_Failure_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes")) / nrow(acm) * 100,  1),
                     n14 = round(nrow(subset(acm, Heart_Failure_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(acm) * 100,  1),
                     n23 = round(nrow(subset(acm, Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes")) / nrow(acm) * 100,  1),
                     n24 = round(nrow(subset(acm, Cardiomyopathy_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(acm) * 100,  1),
                     n34 = round(nrow(subset(acm, Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(acm) * 100, 1),
                     n123 = round(nrow(subset(acm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes")) / nrow(acm) * 100,  1),
                     n124 = round(nrow(subset(acm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(acm) * 100,  1),
                     n134 = round(nrow(subset(acm, Heart_Failure_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(acm) * 100,  1),
                     n234 = round(nrow(subset(acm, Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(acm) * 100,  1),
                     n1234 = round(nrow(subset(acm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(acm) * 100, 1),
                     fill = c("#ffd167", "#4cbd97", "#168ab2", "#ef476f"), 
                     cat.cex = rep(.3, 4), cex = rep(.3, 15), 
                     alpha = rep(0.8, 4), lwd = rep(.2, 4),
                     fontfamily = rep("Helvetica", 15), 
                     cat.fontfamily = rep("Helvetica", 4), 
                     category = c("HF", "CM", "Ischemia", "VA"))
# ggsave("results/figures/ACM_venn.svg", v1)

dcm <- df %>% filter(CM == "DCM")
grid.newpage()
v2 <- draw.quad.venn(area1 = round(nrow(subset(dcm, Heart_Failure_sum == "Yes")) / nrow(dcm) * 100, 1),
                     area2 = round(nrow(subset(dcm, Cardiomyopathy_sum == "Yes")) / nrow(dcm) * 100, 1), 
                     area3 = round(nrow(subset(dcm, Chronic_ischaemic_heart_disease_sum == "Yes")) / nrow(dcm) * 100, 1), 
                     area4 = round(nrow(subset(dcm, Ventricular_arrhythmias_sum == "Yes")) / nrow(dcm) * 100, 1),
                     n12 = round(nrow(subset(dcm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes")) / nrow(dcm) * 100, 1), 
                     n13 = round(nrow(subset(dcm, Heart_Failure_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes")) / nrow(dcm) * 100, 1), 
                     n14 = round(nrow(subset(dcm, Heart_Failure_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(dcm) * 100, 1), 
                     n23 = round(nrow(subset(dcm, Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes")) / nrow(dcm) * 100, 1), 
                     n24 = round(nrow(subset(dcm, Cardiomyopathy_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(dcm) * 100, 1), 
                     n34 = round(nrow(subset(dcm, Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(dcm) * 100, 1),
                     n123 = round(nrow(subset(dcm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes")) / nrow(dcm) * 100, 1), 
                     n124 = round(nrow(subset(dcm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(dcm) * 100, 1), 
                     n134 = round(nrow(subset(dcm, Heart_Failure_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(dcm) * 100, 1), 
                     n234 = round(nrow(subset(dcm, Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(dcm) * 100, 1), 
                     n1234 = round(nrow(subset(dcm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(dcm) * 100, 1),
                     fill = c("#ffd167", "#4cbd97", "#168ab2", "#ef476f"), 
                     cat.cex = rep(.3, 4), cex = rep(.3, 15), 
                     alpha = rep(0.8, 4), lwd = rep(.2, 4),
                     fontfamily = rep("Helvetica", 15), 
                     cat.fontfamily = rep("Helvetica", 4), 
                     category = c("HF", "CM", "Ischemia", "VA"))
# ggsave("results/figures/DCM_venn.svg", v2)

hcm <- df %>% filter(CM == "HCM")
grid.newpage()
v3 <- draw.quad.venn(area1 = round(nrow(subset(hcm, Heart_Failure_sum == "Yes")) / nrow(hcm) * 100, 1),
                     area2 = round(nrow(subset(hcm, Cardiomyopathy_sum == "Yes")) / nrow(hcm) * 100, 1), 
                     area3 = round(nrow(subset(hcm, Chronic_ischaemic_heart_disease_sum == "Yes")) / nrow(hcm) * 100, 1), 
                     area4 = round(nrow(subset(hcm, Ventricular_arrhythmias_sum == "Yes")) / nrow(hcm) * 100, 1),
                     n12 = round(nrow(subset(hcm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes")) / nrow(hcm) * 100, 1), 
                     n13 = round(nrow(subset(hcm, Heart_Failure_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes")) / nrow(hcm) * 100, 1), 
                     n14 = round(nrow(subset(hcm, Heart_Failure_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(hcm) * 100, 1), 
                     n23 = round(nrow(subset(hcm, Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes")) / nrow(hcm) * 100, 1), 
                     n24 = round(nrow(subset(hcm, Cardiomyopathy_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(hcm) * 100, 1), 
                     n34 = round(nrow(subset(hcm, Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(hcm) * 100, 1),
                     n123 = round(nrow(subset(hcm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes")) / nrow(hcm) * 100, 1), 
                     n124 = round(nrow(subset(hcm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(hcm) * 100, 1), 
                     n134 = round(nrow(subset(hcm, Heart_Failure_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(hcm) * 100, 1), 
                     n234 = round(nrow(subset(hcm, Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(hcm) * 100, 1), 
                     n1234 = round(nrow(subset(hcm, Heart_Failure_sum == "Yes" & Cardiomyopathy_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes" & Ventricular_arrhythmias_sum == "Yes")) / nrow(hcm) * 100, 1),
                     fill = c("#ffd167", "#4cbd97", "#168ab2", "#ef476f"), 
                     cat.cex = rep(.3, 4), cex = rep(.3, 15), 
                     alpha = rep(0.8, 4), lwd = rep(.2, 4),
                     fontfamily = rep("Helvetica", 15), 
                     cat.fontfamily = rep("Helvetica", 4), 
                     category = c("HF", "CM", "Ischemia", "VA"))
# ggsave("results/figures/HCM_venn.svg", v3)

ggarrange(v1, v2, v3, labels = c("A) ARVC G+", "B) DCM G+", "C) HCM G+"), ncol = 3,
          font.label = list(size = 4, color = "black", face = "bold", 
                            family = "Helvetica"))
ggsave("results/figures/Figure5.pdf", paper = "a4",
       width = 14 * pix, height = 4 * pix)
rm(v1, v2, v3, acm, dcm, hcm)

# Boxplots ----------------------------------------------------------------

# Create dataframe with desired variables
ex <- c("LV", "RV")
cmr <- df %>% select(starts_with("RV"), starts_with("LV"), contains("Wall_thickness"), 
                     contains("Ecc", ignore.case = FALSE), 
                     contains("Ell", ignore.case = FALSE)) %>% 
  select(-any_of(ex)) %>% names()
new <- df %>% select(f.eid, BSA, CM, Pheno, Sex, any_of(cmr)) %>% 
  filter(Pheno == "Non-Diagnosed") %>% filter(CM %in% c("Controls", "ACM", "DCM", "HCM"))

wt <- new %>% select(1, contains("segment"))
wt$Max_WT <- NA
for (i in 1:nrow(wt)) {
  wt[i, "Max_WT"] <- max(wt[i, 2:(ncol(wt)-1)])
  wt$Max_WT <- as.numeric(wt$Max_WT)
}
wt <- wt %>% select(1, Max_WT)
new <- merge(new, unique(wt))
new$Max_WT[new$Max_WT == 0] <- NA
rm(wt, cmr, ex, i)

# Make plots
ggplot(new, aes(x = Max_WT, fill = Sex)) +
  geom_histogram(position = "dodge") + 
  labs(x = "Maximum LV wall thickness (mm)", y = "Count",
       title = "Distribution of Maximum Left Ventricular wall thickness") +
  scale_fill_manual(values = c("tomato4", "cornflowerblue")) +
  my_theme() 
ggsave("/hpc/dhl_ec/data/uk_biobank/projects/LoF_CMR/results/figures/Wall_thickness_all.svg",
       width = 6 * pix, height = 5  * pix)

p1 <- ggplot(new, 
             aes(x = factor(CM, levels = c("Controls", "ACM", "DCM", "HCM"),
                            labels = c("Controls G-P-", "ARVC G+P-", "DCM G+P-", 
                                       "HCM G+P-")), 
                 y = LVEF, fill = CM)) +
  geom_boxplot(outlier.shape = 18, outlier.size = .4, alpha = 0.9, 
               lwd = .2, fatten = .8) +
  labs(x = "Cardiomyopathy", y = "LVEF (%)") +
  # annotate("segment", x = c(1, 1, 3), xend = c(3, 1, 3),
  #          y = c(90, 88, 88), yend = c(90, 90, 90), lwd = .2) +
  # annotate("text", x = 2, y = 92, label = "p=0.009", size = 1) +
  scale_fill_manual(breaks = c("Controls", "ACM", "DCM", "HCM"),
                    values = c("#4cbd97", "#ffd167", "#168ab2", "#ef476f")) +
  my_theme() +
  theme(legend.position="none", text = element_text(size = 3),
        line = element_line(size = .2))

p1.1 <- ggplot(new, 
               aes(x = factor(CM, levels = c("Controls", "ACM", "DCM", "HCM"),
                              labels = c("Controls G-P-", "ARVC G+P-", "DCM G+P-", 
                                         "HCM G+P-")), 
                   y = RVEF, fill = CM)) +
  geom_boxplot(outlier.shape = 18, outlier.size = .4, alpha = 0.9, 
               lwd = .2, fatten = .8) +
  labs(x = "Cardiomyopathy", y = "RVEF (%)") +
  # annotate("segment", x = c(1, 1, 4), xend = c(4, 1, 4),
  #          y = c(92, 90, 90), yend = c(92, 92, 92), lwd = .2) +
  # annotate("text", x = 2.5, y = 94, label = "p=0.034", size = 1) +
  scale_fill_manual(breaks = c("Controls", "ACM", "DCM", "HCM"),
                    values = c("#4cbd97", "#ffd167", "#168ab2", "#ef476f")) +
  my_theme() +
  theme(legend.position="none", text = element_text(size = 3),
        line = element_line(size = .2))

p2 <- ggplot(new, 
             aes(x = factor(CM, levels = c("Controls", "ACM", "DCM", "HCM"),
                            labels = c("Controls G-P-", "ARVC G+P-", "DCM G+P-", 
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
                              labels = c("Controls G-P-", "ARVC G+P-", "DCM G+P-", 
                                         "HCM G+P-")), 
                   y = RVEDVi, fill = CM)) +
  geom_boxplot(outlier.shape = 18, outlier.size = .4, alpha = 0.9, 
               lwd = .2, fatten = .8) +
  labs(x = "Cardiomyopathy", y = "RVEDVi (ml/m2)")  +
  # annotate("segment", x = c(1, 1, 3), xend = c(3, 1, 3),
  #          y = c(166, 160, 160), yend = c(166, 166, 166), lwd = .2) +
  # annotate("text", x = 2, y = 170, label = "p=0.048", size = 1) +
  scale_fill_manual(breaks = c("Controls", "ACM", "DCM", "HCM"),
                    values = c("#4cbd97", "#ffd167", "#168ab2", "#ef476f")) +
  my_theme() +
  theme(legend.position="none", text = element_text(size = 3),
        line = element_line(size = .2))

p3 <- ggplot(new, 
             aes(x = factor(CM, levels = c("Controls", "ACM", "DCM", "HCM"),
                            labels = c("Controls G-P-", "ARVC G+P-", "DCM G+P-", 
                                       "HCM G+P-")), 
                 y = Max_wall_thickness, fill = CM)) +
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
                            labels = c("Controls G-P-", "ARVC G+P-", "DCM G+P-", 
                                       "HCM G+P-")), 
                 y = peakEll4Ch, fill = CM)) +
  geom_boxplot(outlier.shape = 18, outlier.size = .4, alpha = 0.9, 
               lwd = .2, fatten = .8) +
  labs(x = "Cardiomyopathy", y = "Peak longitudinal strain (%)") +
  # annotate("segment", x = c(1, 1, 3), xend = c(3, 1, 3),
  #          y = c(-5, -7, -7), yend = c(-5, -5, -5), lwd = .2) +
  # annotate("text", x = 2, y = -4, label = "p=0.009", size = 1) +
  scale_fill_manual(breaks = c("Controls", "ACM", "DCM", "HCM"),
                    values = c("#4cbd97", "#ffd167", "#168ab2", "#ef476f")) +
  my_theme() +
  theme(legend.position="none", text = element_text(size = 3),
        line = element_line(size = .2))

ggarrange(p1, p1.1, p2, p2.1, p3, p4, labels = "AUTO", 
          font.label = list(size = 3.5, color = "black", face = "bold", 
                            family = "Helvetica"),
          ncol = 2, nrow = 3)
ggsave("results/figures/Figure6.pdf", paper = "a4",
       width = 8 * pix, height = 12 * pix)
rm(p1, p1.1, p2, p2.1, p3, p4, new)


# Incidence matrix all tests ----------------------------------------------

### Prepare data ###
# Make correct format
dat <- read.delim("~/surfdrive/Mendelian_CM_WES_UKB/Difference_testing_v4.tsv") %>%
  select(Phenotype, p, CM) %>% filter(!is.na(p)) %>%
  tidyr::pivot_wider(names_from = Phenotype, values_from = p)

# Change column names for beauty
dat$CM <- as.factor(dat$CM)
levels(dat$CM) <- c("ACM G+", "ACM/DCM G+ overlap", "Undiagnosed ACM G+", 
                    "DCM G+", "DCM/ACM G+ overlap", "DCM/HCM G+ overlap", 
                    "Undiagnosed DCM G+", "All diagnosed G+", "HCM G+", 
                    "HCM/DCM G+ overlap", "Undiagnosed HCM G+", "Strict HCM G+",
                    "Strict undiagnosed HCM G+")
names(dat)[2:19] <- c("BMI", "MET minutes per week for walking", 
                      "MET minutes per week for moderate activity", 
                      "MET minutes per week for vigorous activity", 
                      "Total MET minutes per week", "Total Cholesterol", "HDL", 
                      "LDL", "Mean systolic blood pressure", 
                      "Mean diastolic blood pressure", "ECG heart rate", 
                      "P duration", "P axis", "PQ interval", "QRS duration", 
                      "R axis", "QTC interval", "T axis")
names(dat)[43:44] <- c("LVEDV/RVEDV", "LVESV/RVESV")
names(dat)[c(45:63, 81:83, 85:89, 93)] <- gsub ("_", " ", 
                                                names(dat)[c(45:63, 81:83, 85:89, 93)])
names(dat)[c(72, 75:77, 84, 90, 92, 94, 95)] <- c("Ever smoked", 
                                                  "Family heart disease", 
                                                  "Cardiac problem", 
                                                  "Heart failure",
                                                  "Acute myocardial infarction",
                                                  "Heart arrhythmia", 
                                                  "Cardiovascular death",
                                                  "Heart failure + cardiomyopathy",
                                                  "Phenotype positive")

# Remove and relocate columns
out <- c("ECG heart rate", "Obesity")
move <- c("Hypertension", "Diabetes", "Ever smoked", "Hypercholesterolaemia", 
          "Family heart disease")
dat <- dat %>% select(-any_of(out)) %>% relocate(any_of(move), .before = BMI)

# Include categories
cat <- c(rep("RISK FACTORS", 15), rep("ECG", 7),
         rep("CMR MEASUREMENTS", 50), rep("CARDIAC OUTCOMES", 20))

pdf("~/surfdrive/Mendelian_CM_WES_UKB/Figures and Tables/Incidence_matrix_v4.pdf",
    height = 17 * pix, width = 7 * pix, paper = "a4")
inc_mat(dat, sig1 = 0.05/nrow(dat)/(ncol(dat) - 1), xas = "CM", legend = "right", cat = cat)
dev.off()
ggsave("~/surfdrive/Mendelian_CM_WES_UKB/Figures and Tables/Incidence_matrix_v4.svg",
       height = 16 * pix, width = 5 * pix)
