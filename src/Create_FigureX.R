#!/usr/bin/env Rscript
##
## Script information ----------------------------------------------------------
##
## Script name: Create_FigureX.R
##
## Purpose of script: Create incidence matrix showing all difference testing
##                      results in one overview
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

suppressMessages(library(ggpubr))
source("src/functions.R")


## Loading arguments -----------------------------------------------------------

# Arguments expected:
#     #1 -- Path to and name of the file with all difference testing results,
#           for example /hpc/dhl_ec/mvanvugt/results/output/TableS6_S7.tsv
#     #2 -- Path to and name of the output file,
#           for example /hpc/dhl_ec/mvanvugt/results/figures/FigureX.pdf

args = commandArgs(trailingOnly = TRUE)
input    = args[1]
output   = args[2]


## Loading files ---------------------------------------------------------------

message(paste0("Loading data from ", input))
dat <- read.delim(input) %>%
    select(Phenotype, Estimate, p.value, CM) %>% filter(!is.na(p.value))
dat$p.value[dat$Estimate <= 1] <- -dat$p.value[dat$Estimate <= 1]
dat <- dat %>% select(Phenotype, p.value, CM) %>% 
    tidyr::pivot_wider(names_from = Phenotype, values_from = p.value)


## Prepare data ----------------------------------------------------------------

message("Preparing data")
# Change column names for beauty
dat$CM <- as.factor(dat$CM)
names(dat)[2:19] <- c("BMI", "MET minutes per week for walking", 
                      "MET minutes per week for moderate activity", 
                      "MET minutes per week for vigorous activity", 
                      "Total MET minutes per week", "Total Cholesterol", "HDL", 
                      "LDL", "Mean systolic blood pressure", 
                      "Mean diastolic blood pressure", "ECG heart rate", 
                      "P duration", "P axis", "PQ interval", "QRS duration", 
                      "R axis", "QTC interval", "T axis")
names(dat)[44:45] <- c("LVEDV/RVEDV", "LVESV/RVESV")
names(dat)[c(46:64, 82, 83, 86:90)] <- gsub ("_", " ", 
                                             names(dat)[c(46:64, 82, 83, 86:90)])
names(dat)[c(73, 76:78, 84, 85, 91:96)] <- c("Ever smoked", 
                                             "Family heart disease", 
                                             "Cardiac problem", 
                                             "Heart failure",
                                             "Chronic ischemic heart disease",
                                             "Acute myocardial infarction",
                                             "Heart arrhythmia", 
                                             "Angina pectoris",
                                             "Cardiovascular death",
                                             "All-cause mortality",
                                             "Heart failure + cardiomyopathy",
                                             "Phenotype positive")

# Remove and relocate columns
out <- c("ECG heart rate", "Obesity")
move <- c("Hypertension", "Diabetes", "Ever smoked", "Hypercholesterolaemia", 
          "Family heart disease")
dat <- dat %>% select(-any_of(out)) %>% relocate(any_of(move), .before = BMI)

# Include categories
cat <- c(rep("RISK FACTORS", 15), rep("ECG", 7),
         rep("CMR MEASUREMENTS", 51), rep("CARDIAC OUTCOMES", 20))


## Create and save figure ------------------------------------------------------

message(paste0("Saving figure in ", output))
pdf(output, height = 17 * pix, width = 7 * pix, paper = "a4")
inc_mat(dat, sig1 = 0.05/nrow(dat)/(ncol(dat) - 1), pdif = "shape", 
        xas = "CM", legend = "right", odif = "color", 
        oname = "Effect direction", cat = cat)
dev.off()

