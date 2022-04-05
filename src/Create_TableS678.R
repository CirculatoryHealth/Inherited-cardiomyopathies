#!/usr/bin/env Rscript
##
## Script information ----------------------------------------------------------
##
## Script name: Create_TableS678.R
##
## Purpose of script: Test for differences in outcomes between the G+ and G- 
##                      groups (Table S6 and S7)
##
## Author: M. van Vugt
##
## Date Created: 2021/06/16
##
## Copyright (c) M. van Vugt, 2022
## Email: m.vanvugt-2@umcutrecht.nl
##

## Options ---------------------------------------------------------------------

options(scipen = 6, digits = 4) # view outputs in non-scientific notation


## Loading packages ------------------------------------------------------------

suppressMessages(library(plyr))
source("src/functions.R")
suppressMessages(library(data.table))
suppressMessages(library(rstatix))


## Loading arguments -----------------------------------------------------------

# Arguments expected:
#     #1 -- Path to and name of the phenotype file,
#           for example /hpc/dhl_ec/mvanvugt/UKB/data/processed/MCM_final_pheno.tsv
#     #2 -- Path to and name of the first output file,
#           for example /hpc/dhl_ec/mvanvugt/UKB/results/output/TableS6_S7.tsv
#     #3 -- Path to and name of output file without overlapping genes,
#           for example /hpc/dhl_ec/mvanvugt/UKB/results/output/TableS8.tsv

args = commandArgs(trailingOnly = TRUE)
input    = args[1]
output   = args[2]
out2     = args[3]


## Loading files ---------------------------------------------------------------

message(paste0("Loading data from ", input))
df <- data.table(read.delim(input, header = TRUE, sep = "\t", stringsAsFactors = FALSE))


## Make strict HCM group -------------------------------------------------------

message("Preparing data")
# For subanalysis without the two most frequent HCM SNPs, we remove these from 
# the dataframe. 
hcm <- df %>% filter(CM == "HCM") %>% 
    filter(!SNP %in% c("1:201359245:G:A", "11:47332274:D:25"))
hcm$CM <- "strict HCM"
df <- rbind(df, hcm)

# Change ACM to ARVC
df$CM[df$CM == "ACM"] <- "ARVC"
df$CM <- as.factor(df$CM)


## Testing differences ---------------------------------------------------------

# Select all variables for which to perform tests
cmr <- df %>% select(starts_with("RV"), starts_with("LV"), contains("Wall_thickness"), 
                     contains("Ecc", ignore.case = FALSE), 
                     contains("Ell", ignore.case = FALSE)) %>% 
select(-c("RV", "LV")) %>% names()
ecg <- c("ECG_heart_rate.0_mean", "P_duration", "P_axis.2.0",
         "PQ_interval.2.0", "QRS_duration", "R_axis.2.0",
         "QTC_interval.2.0", "T_axis.2.0")
met <- df %>% select(contains("MET", ignore.case = FALSE)) %>% names()
bp <- df %>% select(Total_Cholesterol, HDL, LDL,
                    contains("blood_pressure_mean")) %>% names()
cvd <- df %>% select(ends_with("sum"), HFCM, Pheno) %>% names()
dia <- c("Heart_Failure_sum", "Cardiomyopathy_sum", "HCM_sum", "DCM_sum",
         "Chronic_ischaemic_heart_disease_sum", "HFCM", "Pheno")
cols <- c("BMI", met, bp, ecg, cmr, cvd)

# Make list of continuous variables
cont <- cols[!cols %in% cvd]

# Prepare results dataframe
pval <- data.frame(Phenotype = rep("Outcomes", 8), Estimate = rep(NA, 8), 
                   LCI = rep(NA, 8), UCI = rep(NA, 8), p = rep(NA, 8),
                   test = rep(NA, 8), 
                   CM = c("ARVC G+ vs G-", "ARVC G+P- vs G-P-", 
                          "DCM G+ vs G-", "DCM G+P- vs G-P-", "HCM G+ vs G-", 
                          "HCM G+P- vs G-P-", "strict HCM G+ vs G-", 
                          "strict HCM G+P- vs G-P-"))
prim <- c("ARVC", "DCM", "HCM", "strict HCM")

message("Testing differences")
for (cm in prim) {

    # Make df of only this CM and controls
    f <- df %>% select(any_of(cols), CM) %>% 
        filter(CM %in% c(cm,  "Controls"))
    f$CM <- droplevels(f$CM)

    # Make df of only non-diagnosed subjects
    fh <- f %>% filter(Pheno == "Non-Diagnosed")
    fh$CM <- droplevels(fh$CM)

    # Iterate over outcomes
    for (x in cols) {
        if (x %in% cvd) {

            # Select data and change coding
            tmp <- f %>% dplyr::select(any_of(x), CM) %>% as.data.frame
            if (x == "Pheno") {
                tmp[, x] <- ifelse(tmp[, x] == "Diagnosed", 1, 0)
            } else {
                tmp[, x] <- ifelse(tmp[, x] == "Yes", 1, 0)
            }
            tmp$CM <- ifelse(tmp$CM == "Controls", 0, 1)

            tried <- try(fisher.test(table(tmp), workspace = 1e9))

            if (inherits(tried, "try-error")) {
                # If the test fails, replace OR, CI and p-value with NA
                ps <- c(gsub("_sum", "", x), NA, NA, NA, "Fisher", 
                        paste0(cm, " G+ vs G-"))
            } else {
                # Perform test
                test <- fisher.test(table(tmp), workspace = 1e9)
                # Save OR and confidence interval
                or <- test$estimate
                ci <- test$conf.int

                # Save results in dataframe
                ps <- c(gsub("_sum", "", x), or, ci, test$p.value, "Fisher", 
                        paste0(cm, " G+ vs G-"))
            }
            pval <- rbind(pval, ps)

            if (!x %in% dia) {

                # Select data and change coding
                tmp <- fh %>% dplyr::select(any_of(x), CM) %>% as.data.frame
                tmp[, x] <- ifelse(tmp[, x] == "Yes", 1, 0)
                tmp$CM <- ifelse(tmp$CM == "Controls", 0, 1)

                tried <- try(fisher.test(table(tmp), workspace = 1e9))

                if (inherits(tried, "try-error")) {

                    # If the test fails, replace OR, CI and p-value with NA
                    ps <- c(gsub("_sum", "", x), NA, NA, NA, "Fisher", 
                            paste0(cm, " G+P- vs G-P-"))
                } else {
                    # Perform Fisher's exact test on outcome x for undiagnosed subjects
                    test <- fisher.test(table(tmp), workspace = 1e9)
                    # Save OR and confidence interval
                    or <- test$estimate
                    ci <- test$conf.int
                    # Save results in dataframe
                    ps <- c(gsub("_sum", "", x), or, ci, test$p.value, "Fisher", 
                            paste0(cm, " G+P- vs G-P-"))
                }
                pval <- rbind(pval, ps)

            } # Only variables not used in diagnosis test on diagnosed groups

        } else if (x %in% cont) {

            if (!x %in% c(cmr, ecg)) {
                # See if test fails or not
                tried <- try(wilcox.test(get(x) ~ CM, data = f), silent = TRUE)

                if (inherits(tried, "try-error")) {
                    # If test fails, replace Estimate, CI and p-value with NA
                    vals <- c(x, NA, NA, NA, NA, "MWU", cm)
                } else {
                    # Otherwise, perform test
                    test <- wilcox.test(get(x) ~ CM, data = f, conf.int = TRUE)
                    # Calculate effect size 
                    eff <- wilcox_effsize(f, as.formula(paste0(x, " ~ CM")))
                    # Save results 
                    vals <- c(x, eff$effsize, test$conf.int, test$p.value, 
                              "MWU", paste0(cm, " G+ vs G-"))
                }
                pval <- rbind(pval, vals)
            } # Don't test ECG and CMR differences for full groups

            # See if test fails or not
            tried <- try(wilcox.test(get(x) ~ CM, data = fh), silent = TRUE)

            if (inherits(tried, "try-error")) {
                # If test fails, replace Estimate, CI and p-value with NA
                vals <- c(x, NA, NA, NA, NA, "MWU", 
                          paste0(cm, " G+P- vs G-P-"))
            } else {
                # Otherwise, perform test
                test <- wilcox.test(get(x) ~ CM, data = fh, conf.int = TRUE)
                # Calculate effect size 
                eff <- wilcox_effsize(fh, as.formula(paste0(x, " ~ CM")))
                # Save results 
                vals <- c(x, eff$effsize, test$conf.int, test$p.value, "MWU", 
                          paste0(cm, " G+P- vs G-P-"))
            }
            pval <- rbind(pval, vals)

        } # End choose type of test (Fisher or MWU)

    } # End iteration outcomes

}

# Clean up df
pval[, c("Phenotype", "test", "CM")] <- lapply(pval[, c("Phenotype", "test", "CM")], as.factor)
pval[, c("Estimate", "LCI", "UCI", "p")] <- lapply(pval[, c("Estimate", "LCI", "UCI", "p")], as.numeric)
pval <- pval %>% mutate(across(where(is.numeric), round, 3))
pval$p <- ifelse(pval$p < 0.001, signif(pval$p, 2), round(pval$p, 3)) 
names(pval) <- c("Phenotype", "Estimate", "95% LCI", "95% UCI", "p-value", 
                 "test", "CM")


## Save data -------------------------------------------------------------------

message(paste0("Saving results in ", output))
write.table(pval, output, sep = "\t", col.names = TRUE, quote = FALSE, 
            row.names = FALSE)


## Prepare data without overlapping genes --------------------------------------

message("Removing overlapping genes")
# To perform testing on cardiomyopathies excluding the overlapping genes,
# we make subgroups for these 
# vectors for overlapping genes between ACM and DCM and DCM and HCM
dag <- c("DES", "DSP", "PLN")
hag <- c("ACTC1", "JPH2", "MYH7", "TNNC1", "TNNI3", "TNNT2", "TPM1")

# make list containing the subset dfs
dfs <- list()
dfs[["ARVC G+ vs G-"]] <- df %>% filter(CM %in% c("ARVC", "Controls")) %>% filter(!Gene_1 %in% dag)
# DCM has overlapping genes with both ACM and HCM, so we create two datasets
dfs[["DCM G+ vs G- (without ARVC genes)"]] <- df %>% 
    filter(CM %in% c("DCM", "Controls")) %>% filter(!Gene_1 %in% dag)
dfs[["DCM G+ vs G- (without HCM genes)"]] <- df %>% 
    filter(CM %in% c("DCM", "Controls")) %>% filter(!Gene_1 %in% hag)
dfs[["HCM G+ vs G-"]] <- df %>% filter(CM %in% c("HCM", "Controls")) %>% 
    filter(!Gene_1 %in% hag)

# Define outcomes to test for
cvd <- df %>% select(ends_with("sum")) %>% names()


## Test differences ------------------------------------------------------------

message("Testing differences without overlapping genes")
pval <- data.frame()
for (cm in names(dfs)) {

    f <- dfs[[cm]]

    for (x in cvd) {

        # Select correct data
        tmp <- f %>% select(CM, !!x) %>% as.data.frame
        tmp <- droplevels(tmp)
        # Change coding
        tmp[, x] <- ifelse(tmp[, x] == "Yes", 1, 0)
        tmp$CM <- ifelse(tmp$CM == "Controls", 0, 1)

        # Perform Fisher's exact test
        test <- fisher.test(table(tmp), workspace = 1e9)
        # Save OR and confidence interval
        or <- test$estimate
        ci <- test$conf.int
        # Save results in dataframe
        ps <- c(gsub("_sum", "", x), or, ci, test$p.value, cm)
        pval <- rbind(pval, ps)

    } # end iteration outcomes

} # end iteration dataframes

# Clean up df
names(pval) <- c("Phenotype", "OR", "95% LCI", "95% UCI", "pvalue", "CM") 
pval[names(pval)[2:4]] <- lapply(pval[names(pval)[2:4]], as.numeric)
pval[names(pval)[c(1, 6)]] <- lapply(pval[names(pval)[c(1, 6)]], as.factor)
pval <- pval %>% mutate(across(where(is.numeric), round, 3))
pval$pvalue <- as.numeric(pval$pvalue)
pval$pvalue <- ifelse(pval$pvalue < 0.001, signif(pval$pvalue, 2), round(pval$pvalue, 3)) 
names(pval)[5] <- "p-value"


## Save data -------------------------------------------------------------------

message(paste0("Saving results in ", out2))
write.table(pval, out2, sep = "\t", row.names = FALSE, quote = FALSE, 
            col.names = TRUE)


