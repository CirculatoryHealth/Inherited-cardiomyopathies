#!/usr/bin/env Rscript
##
## Script information ----------------------------------------------------------
##
## Script name: Create_Figure5.R
##
## Purpose of script: Create boxplots of some CMR parameters (Figure 5)
##
## Author: M. van Vugt
##
## Date Created: 2022/04/15
##
## Copyright (c) M. van Vugt, 2022
## Email: m.vanvugt-2@umcutrecht.nl
##

## Options ---------------------------------------------------------------------

options(scipen = 6, digits = 4) # view outputs in non-scientific notation
pix <- 0.393700787


## Loading packages ------------------------------------------------------------

suppressMessages(library(dplyr))
source("src/functions.R")
suppressMessages(library(ggpubr))
# Set ggplot theme for all plots created in this script
theme_set(my_theme())


## Loading arguments -----------------------------------------------------------

# Arguments expected:
#     #1 -- Path to and name of the full phenotype file,
#           for example /hpc/dhl_ec/mvanvugt/UKBB/data/processed/MCM_final_pheno.tsv
#     #2 -- Path to and name of the output file,
#           for example /hpc/dhl_ec/mvanvugt/UKBB/results/figures/Figure5.pdf

args = commandArgs(trailingOnly = TRUE)
input    = args[1]
output   = args[2]


## Loading files ---------------------------------------------------------------

message(paste0("Loading data from ", input))
df <- read.delim(input)


## Prepare data ----------------------------------------------------------------

message("Prepare data")
# Create dataframe with desired variables
ex <- c("LV", "RV")
cmr <- df %>% select(starts_with("RV"), starts_with("LV"), contains("Wall_thickness"), 
                     contains("Ecc", ignore.case = FALSE), 
                     contains("Ell", ignore.case = FALSE)) %>% 
select(-any_of(ex)) %>% names()
# Create df with only P- participants
new <- df %>% select(f.eid, BSA, CM, Pheno, Sex, any_of(cmr)) %>% 
    filter(Pheno == "Non-Diagnosed") %>% 
    filter(CM %in% c("Controls", "ACM", "DCM", "HCM"))


## Create boxplots -------------------------------------------------------------

# LVEF boxplot
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
theme(legend.position="none", text = element_text(size = 3),
      line = element_line(size = .2))

# RVEF boxplot
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
theme(legend.position="none", text = element_text(size = 3),
      line = element_line(size = .2))

# LVEDVi boxplot
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
theme(legend.position="none", text = element_text(size = 3),
      line = element_line(size = .2))

# RVEDVi boxplot
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
theme(legend.position="none", text = element_text(size = 3),
      line = element_line(size = .2))

# Max Wall thickness boxplot
p3 <- ggplot(new, 
             aes(x = factor(CM, levels = c("Controls", "ACM", "DCM", "HCM"),
                            labels = c("Controls G-P-", "ARVC G+P-", "DCM G+P-", 
                                       "HCM G+P-")), 
                 y = Maximum_wall_thickness, fill = CM)) +
geom_boxplot(outlier.shape = 18, outlier.size = .4, alpha = 0.9, 
             lwd = .2, fatten = .8) +
labs(x = "Cardiomyopathy", y = "Maximum wall thickness (mm)") +
scale_fill_manual(breaks = c("Controls", "ACM", "DCM", "HCM"),
                  values = c("#4cbd97", "#ffd167", "#168ab2", "#ef476f")) +
theme(legend.position="none", text = element_text(size = 3),
      line = element_line(size = .2))

# Peak longitudinal strain boxplot
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
theme(legend.position="none", text = element_text(size = 3),
      line = element_line(size = .2))


## Combine and save figure -----------------------------------------------------

message(paste0("Combine and save figure in ", output))
ggarrange(p1, p1.1, p2, p2.1, p3, p4, labels = "AUTO", 
          font.label = list(size = 3.5, color = "black", face = "bold", 
                            family = "Helvetica"),
          ncol = 2, nrow = 3)
ggsave(output, paper = "a4", width = 8 * pix, height = 12 * pix)

