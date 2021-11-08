library(plotly)

m <- list(l = 50, r = 50, b = 300, t = 100, pad = 4)

pie <- plot_ly(tmp, 
               labels = paste0(tmp$Gene, "\n", tmp$N, " (", round(tmp$prop), "%)"),
               values = ~N, type = "pie", showlegend = FALSE,
               textposition = ifelse(tmp$prop < 4, "outside", "inside"), 
               textinfo = "label",
               insidetextfont = list(color = "black"),
               marker = list(colors = viridis(nrow(tmp), alpha = 0.8),
                             line = list(color = "#FFFFFF", width = 1)))
pie <- pie %>% layout(title = paste0("Mutated genes in ", cm),
                      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                      autosize = FALSE, margin = m)
pie


# Create dataframe with desired variables
cmr <- df %>% select(starts_with("RV"), starts_with("LV"), contains("Wall_thickness"), 
                     contains("Ecc", ignore.case = FALSE), 
                     contains("Ell", ignore.case = FALSE)) %>% names()
new <- df %>% select(f.eid, BSA, CM, Pheno, any_of(cmr)) %>% filter(Pheno == "Non-Diagnosed")

# Create right variables
new$LVEDVi <- new$LVEDV / new$BSA
new$RVEDVi <- new$RVEDV / new$BSA

wt <- new %>% select(f.eid, contains("segment"))
wt$Max_WT <- NA
for (i in 1:nrow(wt)) {
  wt[i, "Max_WT"] <- max(wt[i, 2:(ncol(wt)-1)])
}
wt <- wt %>% select(f.eid, Max_WT)
new <- merge(new, wt)
new$Max_WT[new$Max_WT == 0] <- NA

# Make plots
library(ggplot2)
library(ggpubr)

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
                      y = Max_WT)) +
  geom_boxplot(outlier.shape = 18) +
  labs(x = "Cardiomyopathy", y = "Maximum wall thickness (mm)") +
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
ggsave("CMR_box.svg", width = 20, height = 30)


# Make Venn Diagram
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
               lty = "blank", cat.cex = rep(3, 4), cex = rep(3, 15),
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
               lty = "blank", cat.cex = rep(3, 4), cex = rep(3, 15),
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
               lty = "blank", cat.cex = rep(3, 4), cex = rep(3, 15),
               category = c("HF", "CM", "Ischaemia", "VA"))
# ggsave("results/figures/HCM_venn.svg", v3)

ggarrange(v1, v2, v3, labels = c("A) ACM", "B) DCM", "C) HCM"), ncol = 3)
ggsave("results/figures/Fig4-Venn.svg", width = 30, height = 8)
