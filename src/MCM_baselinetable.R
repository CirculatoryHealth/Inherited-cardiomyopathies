## Script information ---------------------------
##
## Script name: MCM_baselinetable.R
##
## Purpose of script: Make a baseline table of the UKB-phenotypes
##
## Author: M. van Vugt
##
## Date Created: 2021-07-15
##
## Copyright (c) M. van Vugt, 2021
## Email: m.vanvugt-2@umcutrecht.nl
##

## Options ---------------------------

options(scipen = 6, digits = 4) # view outputs in non-scientific notation


## Loading arguments --------------------------

# Arguments expected:
#     #1 -- Path to and name of the phenotype file,
#           for example /hpc/dhl_ec/mvanvugt/UKBB/Project1_ukb_phenotypes.tab
#     #3 -- Output directory, for example /hpc/dhl_ec/mvanvugt/UKBB
#     #4 -- Prefix of the output files

args = commandArgs(trailingOnly = TRUE)
input = args[1] # "data/processed/MCM_final_pheno.txt"
output = args[2] # "results/output/"
prefix = args[3] # "MCM"


## Loading packages ---------------------------

library(dplyr)
library(data.table)
library(tableone)
library(ggpubr)


# Loading data ------------------------------------------------------------

message("Loading data")
df <- data.table(read.table(input, sep = "\t", header = TRUE, quote = ""))


# Making Baselinetable ----------------------------------------------------

# The vectors are created with column names which we want to include in the
# baseline-table. If a column name is not represented in the vector cols, it
# will not be included in the baseline-table, so this should be changed here.
cmr <- df %>% select(starts_with("RV"), starts_with("LV")) %>% names()
ecg <- c("ECG_heart_rate.0_mean", "P_duration", "P_axis.2.0",
         "PQ_interval.2.0", "QRS_duration", "R_axis.2.0",
         "QTC_interval.2.0", "T_axis.2.0")
met <- df %>% select(starts_with("MET")) %>% names()
bp <- df %>% select(Total_Cholesterol, HDL, LDL,
                    contains("blood_pressure_mean")) %>% names()
diag <- df %>% select(ends_with("sum")) %>% names()
cols <- c("Sex", "Age_when_attended_assessment_centre.0.0", "Ethnicity",
          "BMI", met, bp, diag, ecg, cmr, "CM")

# A new df (df) is created, with all columns listed in cols
message("Selecting variables for baseline table")
df <- df %>% dplyr::select(f.eid, any_of(cols), ends_with("FirstDate"))


# Now we'll create vectors with column names for the baseline-table.
# fac.col contains all column names that should be regarded as factors
# nn.col contains all columns that should be regarded as nonnormal (for now
# these are all columns except the factors, but can be changed of course).
fac.col <- c("Sex", "Ethnicity", diag, "CM")
df <- as.data.frame(df)
df[fac.col] <- lapply(df[fac.col], as.factor)
nn.col <- cols[!cols %in% fac.col]

# ocol <- vector()
# norm <- list()
# for (n in nn.col) {
#   vec <- na.omit(df_named[n])
#   norm[[n]] <- list()
#   norm[[n]][["dens.plot"]] <- ggdensity(vec[,1])
#   norm[[n]][["qq.plot"]] <- ggqqplot(vec[,1])
#   norm[[n]][["test"]] <- ks.test(vec[,1], "pnorm")
#
#   if (norm[[n]][["test"]]$p.value < 0.05) {
#     ocol <- c(ocol, n)
#   } # Finish if-loop adding nonnormal variables
# } # Finish for-loop iterating over non-factor columns

# After manual inspection, some variables are removed anyways
# rem <- c("Total_Cholesterol", "LDL", "Systolic_blood_pressure_mean",
#          "Diastolic_blood_pressure_mean", "ECG_heart_rate.0_mean", "P_duration",
#          "RVEF", "RVPFR", "LVPER", "LVPFR")
# ocol <- ocol[!ocol %in% rem]

message("Creating the baseline table")
tab1 <- CreateTableOne(vars = cols, data = df, factorVars = fac.col,
                       strata = "CM", addOverall = FALSE)
ex1 <- print(tab1, showAllLevels = FALSE, missing = TRUE, nonnormal = nn.col,
             formatOptions = list(big.mark = ","), quote = FALSE,
             noSpaces = TRUE, printToggle = FALSE)

message("Saving data")
write.csv(ex1, paste0(output, prefix, "_Table1.csv"))

# And last but not least, the final dataframe used to create the baseline-
# table is exported as well.
write.table(df, paste0(output, prefix, "_table1_data.tsv"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
saveRDS(df, paste0(output, prefix, "_table1_data.rds"))
