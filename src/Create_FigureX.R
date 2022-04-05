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
#           for example /hpc/dhl_ec/mvanvugt/results/output/
#     #2 -- Path to and name of the output file,
#           for example /hpc/dhl_ec/mvanvugt/results/figures/FigureX.pdf

args = commandArgs(trailingOnly = TRUE)
input    = args[1]
output   = args[2]


## Loading files ---------------------------------------------------------------

message(paste0("Loading data from ", input))
dat <- read.delim(input) %>%
    select(Phenotype, p, CM) %>% filter(!is.na(p)) %>%
    tidyr::pivot_wider(names_from = Phenotype, values_from = p)


## Prepare data ----------------------------------------------------------------

message("Preparing data")
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


## Create and save figure ------------------------------------------------------

message(paste0("Saving figure in ", output))
pdf(output, height = 17 * pix, width = 7 * pix, paper = "a4")
inc_mat(dat, sig1 = 0.05/nrow(dat)/(ncol(dat) - 1), xas = "CM", 
        legend = "right", cat = cat)
dev.off()

# Also save as svg
output <- gsub("pdf", "svg", output)
message(paste0("Also saving figure in ", output))
ggsave(output, height = 16 * pix, width = 5 * pix)



