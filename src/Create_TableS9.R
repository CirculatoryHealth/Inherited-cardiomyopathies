#!/usr/bin/env Rscript
##
## Script information ----------------------------------------------------------
##
## Script name: Create_TableS9.R
##
## Purpose of script: Test differences in cardiovascular outcomes between 
##                      variant carriers and controls stratified by gene
##
## Author: M. van Vugt
##
## Date Created: 2022/03/22
##
## Copyright (c) M. van Vugt, 2022
## Email: m.vanvugt-2@umcutrecht.nl
##

## Options ---------------------------------------------------------------------

options(scipen = 6, digits = 4) # view outputs in non-scientific notation


## Loading packages ------------------------------------------------------------

suppressMessages(library(dplyr))
source("src/functions.R")
suppressMessages(library(data.table))


## Loading arguments -----------------------------------------------------------

# Arguments expected:
#     #1 -- Path to and name of the input file,
#           for example /hpc/dhl_ec/mvanvugt/UKBB/data/processed/MCM_final_pheno.tsv
#     #2 -- Path to and name of the output file,
#           for example /hpc/dhl_ec/mvanvugt/UKBB/results/output/TableS9.tsv

args = commandArgs(trailingOnly = TRUE)
input    = args[1]
output   = args[2]


## Loading files ---------------------------------------------------------------

message(paste0("Loading data from ", input))
df <- data.table(read.delim(input, header = TRUE, sep = "\t", 
                            stringsAsFactors = FALSE))
df$CM <- as.factor(df$CM)


## Testing differences ---------------------------------------------------------

out <- df %>% select(Heart_Failure_sum, Cardiomyopathy_sum, HFCM, Pheno, 
                     Ventricular_arrhythmias_sum, Atrial_arrhythmias_sum, 
                     Heart_Arrhythmia_sum, Chronic_ischaemic_heart_disease_sum, 
                     Angina_sum, Cardiovascular_Death_sum, All_cause_mortality_sum) %>% names
genes <- unique(c(levels(as.factor(df$Gene_1)), levels(as.factor(df$Gene_2))))
pval <- data.frame()
for (cm in c("ACM", "DCM", "HCM")) {

    for (gen in genes) {

        for (dis in out) {

            # Select data 
            tmp <- df %>% select(CM, Gene_1, Gene_2, !!dis) %>%
                filter(CM %in% c(cm, "Controls")) %>%
                filter(CM == "Controls" | Gene_1 == gen | Gene_2 == gen) %>% 
                select(CM, !!dis)
            tmp <- data.frame(droplevels(tmp))

            # Change codings
            if (dis == "Pheno") {
                tmp$Pheno <- ifelse(tmp$Pheno == "Diagnosed", 1, 0)
            } else {
                tmp[, dis] <- ifelse(tmp[, dis] == "Yes", 1, 0)
            } # change coding to 1/0
            tmp$CM <- ifelse(tmp$CM == "Controls", 0, 1)

            if (nrow(table(tmp)) > 1) {
                # Perform test
                test <- fisher.test(table(tmp), workspace = 1e9)
                # Save OR and confidence interval
                or <- test$estimate
                ci <- test$conf.int
                # Save results in dataframe
                ps <- c(gsub("_sum", "", dis), or, ci, test$p.value, cm, gen)
                pval <- rbind(pval, ps)
            }

        } # iterate over outcomes
    } # iterate over genes

} # iterate over the cardiomyopathies

# Clean up the df
names(pval) <- c("Phenotype", "OR", "95% LCI", "95% UCI", "pvalue", "CM", "Gene") 
pval$CM[pval$CM == "ACM"] <- "ARVC G+ vs G-"
pval$CM[pval$CM == "DCM"] <- "DCM G+ vs G-"
pval$CM[pval$CM == "HCM"] <- "HCM G+ vs G-"

pval[names(pval)[2:4]] <- lapply(pval[names(pval)[2:4]], as.numeric)
pval[names(pval)[c(1, 6, 7)]] <- lapply(pval[names(pval)[c(1, 6, 7)]], as.factor)
pval <- pval %>% mutate(across(where(is.numeric), round, 3))
pval$pvalue <- as.numeric(pval$pvalue)
pval$pvalue <- ifelse(pval$pvalue < 0.001, signif(pval$pvalue, 2), round(pval$pvalue, 3)) 


## Save data -------------------------------------------------------------------

message(paste0("Saving results in ", output))
write.table(pval, output, sep = "\t", row.names = FALSE, quote = FALSE, 
            col.names = TRUE)

