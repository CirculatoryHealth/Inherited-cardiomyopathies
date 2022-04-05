#!/usr/bin/env Rscript
##
## Script information ----------------------------------------------------------
##
## Script name: Create_FigureS1.R
##
## Purpose of script: Create figure S1 (Venn Diagrams with overlapping 
##                      phenotypes)
##
## Author: M. van Vugt
##
## Date Created: 2022/04/14
##
## Copyright (c) M. van Vugt, 2022
## Email: m.vanvugt-2@umcutrecht.nl
##

## Options ---------------------------------------------------------------------

options(scipen = 6, digits = 4) # view outputs in non-scientific notation
pix <- 0.393700787


## Loading packages ------------------------------------------------------------

suppressMessages(library(ggpubr))
suppressMessages(library(VennDiagram))
suppressMessages(library(svglite))
source("src/functions.R")


## Loading arguments -----------------------------------------------------------

# Arguments expected:
#     #1 -- Path to and name of the full phenotype file,
#           for example /hpc/dhl_ec/mvanvugt/UKBB/data/processed/MCM_final_pheno.tsv
#     #2 -- Path to and name of the output file,
#           for example /hpc/dhl_ec/mvanvugt/UKBB/results/figures/FigureS1.pdf

args   = commandArgs(trailingOnly = TRUE)
input  = args[1]
output = args[2]


## Loading data ----------------------------------------------------------------

message(paste0("Loading data from ", input))
df <- read.delim(input)


## Make ARVC VennDiagram -------------------------------------------------------

message("Creating figures")
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


## Make DCM VennDiagram --------------------------------------------------------

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


## Make HCM VennDiagram --------------------------------------------------------

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


## Combine and save figure -----------------------------------------------------

message(paste0("Combining and saving figure in ", output))
ggarrange(v1, v2, v3, labels = c("A) ARVC G+", "B) DCM G+", "C) HCM G+"), ncol = 3,
          font.label = list(size = 4, color = "black", face = "bold", 
                            family = "Helvetica"))
ggsave(output, paper = "a4", width = 14 * pix, height = 4 * pix)

