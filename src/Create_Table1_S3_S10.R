#!/usr/bin/env Rscript
##
## Script information ----------------------------------------------------------
##
## Script name: MCM_baselinetable.R
##
## Purpose of script: Make a baseline table of the specified outcomes, 
##                      stratified by CM and diagnosis (Table 1 and S3) and 
##                      for only CMR P- participants (Table S10). Also the
##                      differences for Table S10 are calculated, to make sure
##                      these groups are similar
##
## Author: M. van Vugt
##
## Date Created: 2021/07/15
##
## Copyright (c) M. van Vugt, 2022
## Email: m.vanvugt-2@umcutrecht.nl
##

## Options ---------------------------------------------------------------------

options(scipen = 6, digits = 4) # view outputs in non-scientific notation


## Loading packages ------------------------------------------------------------

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(tableone))
suppressMessages(library(ggpubr))
suppressMessages(library(rstatix))


## Loading arguments -----------------------------------------------------------

# Arguments expected:
#     #1 -- Path to and name of the phenotype file,
#           for example /hpc/dhl_ec/mvanvugt/UKBB/data/processed/MCM_final_pheno.tsv
#     #2 -- Path to output directory, for example results/output
#     #3 -- Prefix of the output, for example MCM

args = commandArgs(trailingOnly = TRUE)
input    = args[1]
output   = args[2]
prefix   = args[3]


## Loading files ---------------------------------------------------------------

message(paste0("Loading data from ", input))
df <- data.table(read.table(input, sep = "\t", header = TRUE, quote = ""))


## Preparing data --------------------------------------------------------------

message("Preparing data")
# Calculating Total MET
df$Total_MET_minutes_per_week <- df$MET_minutes_per_week_for_moderate_activity.0.0 + df$MET_minutes_per_week_for_vigorous_activity.0.0 + df$MET_minutes_per_week_for_walking.0.0

# Categorising Sex in the correct way
df$Sex <- factor(df$Sex, levels = c("Male", "Female"))

# Defining phenotype positive people
# Individuals with a diagnosis for CM / DCM / HCM / HF and not ischemia are
# marked as phenotype positive (Diagnosed)
df$HFCM[df$Heart_Failure_sum == "Yes" | df$Cardiomyopathy_sum == "Yes" | df$HCM_sum == "Yes" | df$DCM_sum == "Yes"] <- "Yes"
df$HFCM[df$Heart_Failure_sum == "No" & df$Cardiomyopathy_sum == "No" & df$HCM_sum == "No" & df$HCM_sum == "No"] <- "No"

df$Pheno[df$Chronic_ischaemic_heart_disease_sum == "No" & df$HFCM == "Yes"] <- "Diagnosed"
df$Pheno[df$Chronic_ischaemic_heart_disease_sum == "Yes" & df$HFCM == "Yes"] <- "Non-Diagnosed"
df$Pheno[df$HFCM == "No"] <- "Non-Diagnosed"
df$Pheno <- as.factor(df$Pheno)

# Calculate some CMR parameters
# Max wall thickness
wt <- df %>% select(1, contains("segment"))
wt$Maximum_wall_thickness <- NA
for (i in 1:nrow(wt)) {
  wt[i, "Maximum_wall_thickness"] <- max(wt[i, 2:(ncol(wt)-1)])
  wt$Maximum_wall_thickness <- as.numeric(wt$Maximum_wall_thickness)
}
wt <- wt %>% select(1, Maximum_wall_thickness)
df <- merge(df, unique(wt))
df$Maximum_wall_thickness[df$Maximum_wall_thickness == 0] <- NA
rm(wt, i)

# Septal wall thickness
df$Septal_wall_thickness <- (df$Wall_thickness_segment_2 + df$Wall_thickness_segment_3 + df$Wall_thickness_segment_8 + df$Wall_thickness_segment_9 + df$Wall_thickness_segment_14) / 5

# BSA-indexed parameters
df$RVEDVi <- df$RVEDV / df$BSA
df$RVESVi <- df$RVESV / df$BSA
df$RVSVi  <- df$RVSV / df$BSA
df$LVEDVi <- df$LVEDV / df$BSA
df$LVESVi <- df$LVESV / df$BSA
df$LVSVi  <- df$LVSV / df$BSA
df$LVEDMi <- df$LVEDM / df$BSA
# CMR ratios
df$LVMVR       <- df$LVEDM / df$LVEDV
df$LVEDV.RVEDV <- df$LVEDV / df$RVEDV
df$LVESV.RVESV <- df$LVESV / df$RVESV

# Preparing vectors for inclusion
# The vectors are created with column names which we want to include in the
# baseline-table. If a column name is not represented in the vector cols, it
# will not be included in the baseline-table, so this should be changed here.
cmr <- df %>% select(starts_with("RV"), starts_with("LV"), 
                     contains("Wall_thickness"), 
                     contains("Ecc", ignore.case = FALSE), 
                     contains("Ell", ignore.case = FALSE)) %>% 
  select(-c("RV", "LV")) %>% names()
ecg <- c("ECG_heart_rate.0_mean", "P_duration", "P_axis.2.0",
         "PQ_interval.2.0", "QRS_duration", "R_axis.2.0",
         "QTC_interval.2.0", "T_axis.2.0")
met <- df %>% select(starts_with("MET")) %>% names()
bp <- df %>% select(Total_Cholesterol, HDL, LDL,
                    contains("blood_pressure_mean")) %>% names()
diag <- df %>% select(ends_with("sum")) %>% names()
cols <- c("Sex", "Age_when_attended_assessment_centre.0.0", "Ethnicity",
          "BMI", met, bp, diag, ecg, cmr, "CM", "Total_MET_minutes_per_week",
          "Pheno", "HFCM", "ECG", "CMR")

# Create column with ECG availability
tmp <- df %>% select(f.eid, any_of(ecg)) %>% select(-ECG_heart_rate.0_mean)
tmp$ECG <- rowMeans(tmp[, 2:ncol(tmp)], na.rm = TRUE)
tmp$ECG <- ifelse(is.na(tmp$ECG), "No", "Yes")
tmp <- tmp %>% select(f.eid, ECG)
df <- merge(df, unique(tmp))

# Save data
write.table(df, input, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")


## Make strict HCM group -------------------------------------------------------

# For subanalysis without the two most frequent HCM SNPs, we remove these from 
# the dataframe. 
hcm <- df %>% filter(CM == "HCM") %>% 
    filter(!SNP %in% c("1:201359245:G:A", "11:47332274:D:25"))
hcm$CM <- "strict HCM G+"
df <- rbind(df, hcm)


## Continue data preparations --------------------------------------------------

# A new df (df) is created, with all columns listed in cols
message("Selecting variables for baseline table")
df <- df %>% dplyr::select(f.eid, any_of(cols), ends_with("FirstDate"))

# Now we'll create vectors with column names for the baseline-table.
# fac.col contains all column names that should be regarded as factors
# nn.col contains all columns that should be regarded as nonnormal
fac.col <- c("Sex", "Ethnicity", diag, "CM", "CMR", "ECG")
df <- as.data.frame(df)
df[fac.col] <- lapply(df[fac.col], as.factor)
nn.col <- cols[!cols %in% fac.col]

# Change ACM to ARVC and add carrier status (G+/G-)
df$CM <- as.character(df$CM)
df$CM[df$CM == "ACM"]      <- "ARVC G+"
df$CM[df$CM == "DCM"]      <- "DCM G+"
df$CM[df$CM == "HCM"]      <- "HCM G+"
df$CM[df$CM == "Controls"] <- "Controls G-"
df$CM <- as.factor(df$CM)


## Create CM-stratified baseline table -----------------------------------------

message("Creating the baseline tables")
tab1 <- CreateTableOne(vars = cols, data = df, factorVars = fac.col,
                       strata = "CM", addOverall = FALSE)
ex1 <- print(tab1, showAllLevels = FALSE, missing = TRUE, nonnormal = nn.col,
             formatOptions = list(big.mark = ","), quote = FALSE,
             noSpaces = TRUE, printToggle = FALSE)


## Create pheno-stratified tables per CM ---------------------------------------

# Subset the df
acm <- subset(df, CM == "ARVC G+")
dcm <- subset(df, CM == "DCM G+")
hcm <- subset(df, CM == "HCM G+")
scm <- subset(df, CM == "strict HCM G+")
con <- subset(df, CM == "Controls G-")

# Create the table
taba <- CreateTableOne(vars = cols, data = acm, factorVars = fac.col,
                       strata = "Pheno", addOverall = TRUE)
exa <- print(taba, showAllLevels = FALSE, missing = TRUE, nonnormal = nn.col,
             formatOptions = list(big.mark = ","), quote = FALSE,
             noSpaces = TRUE, printToggle = FALSE)
tabd <- CreateTableOne(vars = cols, data = dcm, factorVars = fac.col,
                       strata = "Pheno", addOverall = TRUE)
exd <- print(tabd, showAllLevels = FALSE, missing = TRUE, nonnormal = nn.col,
             formatOptions = list(big.mark = ","), quote = FALSE,
             noSpaces = TRUE, printToggle = FALSE)
tabh <- CreateTableOne(vars = cols, data = hcm, factorVars = fac.col,
                       strata = "Pheno", addOverall = TRUE)
exh <- print(tabh, showAllLevels = FALSE, missing = TRUE, nonnormal = nn.col,
             formatOptions = list(big.mark = ","), quote = FALSE,
             noSpaces = TRUE, printToggle = FALSE)
tabs <- CreateTableOne(vars = cols, data = scm, factorVars = fac.col,
                       strata = "Pheno", addOverall = TRUE)
exs <- print(tabs, showAllLevels = FALSE, missing = TRUE, nonnormal = nn.col,
             formatOptions = list(big.mark = ","), quote = FALSE,
             noSpaces = TRUE, printToggle = FALSE)
tabc <- CreateTableOne(vars = cols, data = con, factorVars = fac.col,
                       strata = "Pheno", addOverall = TRUE)
exc <- print(tabc, showAllLevels = FALSE, missing = TRUE, nonnormal = nn.col,
             formatOptions = list(big.mark = ","), quote = FALSE,
             noSpaces = TRUE, printToggle = FALSE)


## Saving tables and df with data ----------------------------------------------

message(paste0("Saving data in ", output, "/*_Table1.csv"))
write.csv(ex1, paste0(file.path(output, prefix), "_Table1.csv"))
write.csv(exa, file.path(output, "ACM_Table1.csv"))
write.csv(exd, file.path(output, "DCM_Table1.csv"))
write.csv(exh, file.path(output, "HCM_Table1.csv"))
write.csv(exs, file.path(output, "strictHCM_Table1.csv"))
write.csv(exc, file.path(output, "Controls_Table1.csv"))

# And last but not least, the final dataframe used to create the baseline-
# table is exported as well.
write.table(df, paste0(file.path(output, prefix), "_table1_data.tsv"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
saveRDS(df, paste0(file.path(output, prefix), "_table1_data.rds"))


## Baseline table for only CMR P- participants ---------------------------------

# A new df (df) is created, with all columns listed in cols
message("Selecting variables for baseline table for CMR P- participants")
f <- df %>% filter(Pheno == "Non-Diagnosed") %>% 
    dplyr::select(f.eid, any_of(cols), ends_with("FirstDate")) %>%
    filter(CMR == "Yes")

# Now we'll create vectors with column names for the baseline-table.
# fac.col contains all column names that should be regarded as factors
# nn.col contains all columns that should be regarded as nonnormal
fac.col <- c("Sex", "Ethnicity", diag, "CM", "ECG")
f <- as.data.frame(f)
f[fac.col] <- lapply(f[fac.col], as.factor)
nn.col <- cols[!cols %in% fac.col]
f$Sex <- factor(f$Sex, levels = c("Male", "Female"))

message("Creating the baseline table for only CMR individuals")
tab1 <- CreateTableOne(vars = cols, data = f, factorVars = fac.col,
                       strata = "CM", addOverall = FALSE)
ex1 <- print(tab1, showAllLevels = FALSE, missing = TRUE, nonnormal = nn.col,
             formatOptions = list(big.mark = ","), quote = FALSE,
             noSpaces = TRUE, printToggle = FALSE)


## Saving results --------------------------------------------------------------

message(paste0("Saving results in ", file.path(output, prefix), "_Table1_CMR.csv"))
write.csv(ex1, paste0(file.path(output, prefix), "_Table1_CMR.csv"))


## Testing differences ---------------------------------------------------------

message("Testing for differences among CMR P- participants")
# Select the outcomes to test for
con <- c("BMI", "Total_MET_minutes_per_week")
diag <- f %>% select(ends_with("sum")) %>% 
    select(!c("Obesity_sum", "Cardiac_arrest_sum", 
              "Conduction_disorders_sum", "Valvular_disease_sum", 
              "Congenital_heart_disease_sum", "Pulmonary_obstructive_disease_sum", 
              "Heart_Arrhythmia_sum")) %>% names

# Test the differences
pval <- data.frame()
for (cm in c("ARVC G+", "DCM G+", "HCM G+")) {

    for (dis in c(con, diag)) {

        # Select data
        tmp <- f %>% filter(CM %in% c(cm, "Controls G-")) %>%
            select(CM, !!dis)
        tmp <- droplevels(tmp)

        if (dis %in% con) {

            # Perform MWU for continuous outcomes
            test <- wilcox.test(get(dis) ~ CM, data = tmp, conf.int = TRUE)
            eff <- wilcox_effsize(tmp, as.formula(paste0(dis, " ~ CM")))
            # Save results
            p <- c(dis, eff$effsize, test$conf.int, test$p.value, cm)

        } else if (dis %in% diag) {

            # Change coding
            tmp[, dis] <- ifelse(tmp[, dis] == "Yes", 1, 0)
            tmp$CM <- ifelse(tmp$CM == "Controls G-", 0, 1)

            # Only perform test if there's observations
            if (nrow(table(tmp)) > 1 & ncol(table(tmp)) > 1) {
                # Perform Fisher's exact test
                test <- fisher.test(table(tmp), workspace = 1e9)
                # Save results
                p <- c(dis, test$estimate, test$conf.int, test$p.value, cm)
            }
        } # decide on test to perform

        # Add results to df
        pval <- rbind(pval, p)
        pval[names(pval)] <- lapply(pval[names(pval)], as.character)

    } # iterate over outcomes

} # iterate over CMs

# Clean up results df
names(pval) <- c("Phenotype", "Estimate", "LCI", "UCI", "pvalue", "CM")
pval[names(pval)[2:4]] <- lapply(pval[names(pval)[2:4]], as.numeric)
pval[names(pval)[c(1, 6)]] <- lapply(pval[names(pval)[c(1, 6)]], as.factor)
pval <- pval %>% mutate(across(where(is.numeric), round, 3))
pval$pvalue <- as.numeric(pval$pvalue)
pval$pvalue <- ifelse(pval$pvalue < 0.001, signif(pval$pvalue, 2), 
                      round(pval$pvalue, 3)) 


## Saving results --------------------------------------------------------------

message(paste0("Saving data in ", file.path(output, prefix), "_Differences_CMR.tsv"))
write.table(pval, paste0(file.path(output, prefix), "_Differences_CMR.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)


