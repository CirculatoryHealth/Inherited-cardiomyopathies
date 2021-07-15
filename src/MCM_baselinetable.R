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
#     #2 -- Directory to the helper files, 
#           for example /hpc/dhl_ec/mvanvugt/Software/UKB-pipeline-Utrecht
#     #3 -- Output directory, for example /hpc/dhl_ec/mvanvugt/UKBB
#     #4 -- Prefix of the output files

args = commandArgs(trailingOnly = TRUE)
input = args[1] # "Data/processed/MCM_final_pheno.txt"
output = args[3] # "results/output/"
prefix = args[4] # "MCM"


## Loading packages ---------------------------

library(dplyr)
library(data.table)
library(tableone)
library(ggpubr)


# Loading data ------------------------------------------------------------

message("Loading data")
df <- data.table(read.table(input, sep = " ", header = TRUE, quote = FALSE))


# Making Baselinetable ----------------------------------------------------

# In this section, the data for the baseline-table is prepared. 
# In this particular loop, the controls-df gets an extra column 
# (Total_variants), because all dfs need the same columns. 


# Sine vectors are created with column names which we want to include in 
# the baseline-table. If a column name is not represented in the vector cols,
# it will not be included in the baseline-table, so this can of course be
# changed here.
cmr <- c("RVEDVi", "RVESVi", "RVSV", "RVEF", "RVPER", "RVPFR", "RVPAFR",
         "LVEDVi", "LVESVi", "LVSV", "LVEF", "LVPER", "LVPFR", "LVPAFR",
         "LVEDMi", "LVEDM/LVEDV", "LV_RV_EDV")
ecg <- c("ECG_heart_rate.0_mean", "P_duration", "P_axis", "PQ_interval", 
         "QRS_duration", "R_axis", "QTC_interval", "T_axis")
met <- df_named %>% select(starts_with("MET")) %>% names()
bp <- df_named %>% select(Total_Cholesterol, HDL, LDL, 
                          contains("blood_pressure_mean")) %>% names()
diag <- df_named %>% select(ends_with("SR")) %>% names()
death <- df_named %>% select(All_cause_mortality,
                             ends_with("Death", ignore.case = FALSE)) %>% 
  names()
cols <- c("Sex", "Age_when_attended_assessment_centre.0.0", "Ethnicity", 
          "BMI", "Ever_Smoked", met, bp, diag, death, ecg, cmr, "CM")

# A new df (df_named) is created, with all columns listed in cols
df_named <- df_named %>% dplyr::select(f.eid, any_of(cols))
df <- df %>% dplyr::select(f.eid, all_of(cmr))
df_named <- merge(df_named, df, by = "f.eid", all = TRUE)

# To be able to distinguish the different phenotypes, an extra column is 
# added to the final df which contains the name of the original df.
df_named$CM <- toupper(names(dfs)[q])
all[[names(dfs)[q]]] <- df_named

# In this loop, we create one df with all phenotypes combined. The last added
# column will help us to still distinguish the phenotypes.
df_named <- data.frame()
for (i in 1:length(all)) {
  df_named <- rbind(df_named, all[[i]])
}
df_named$CM <- as.factor(df_named$CM)
rm(dfs)

# Now we'll create vectors with column names for the baseline-table.
# fac.col contains all column names that should be regarded as factors
# nn.col contains all columns that should be regarded as nonnormal (for now 
# these are all columns except the factors, but can be changed of course).
diag <- df_named %>% 
  select(All_cause_mortality, contains("Death"), ends_with("SR")) %>%
  names()
fac.col <- c("Sex", "Ethnicity", "Ever_Smoked", diag, "CM")
df_named[fac.col] <- lapply(df_named[fac.col], as.factor)
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
rem <- c("Total_Cholesterol", "LDL", "Systolic_blood_pressure_mean",
         "Diastolic_blood_pressure_mean", "ECG_heart_rate.0_mean", "P_duration",
         "RVEF", "RVPFR", "LVPER", "LVPFR")
ocol <- ocol[!ocol %in% rem]

# Then we start creating the baseline table. I noticed that it cannot handle 
# too many columns, so this is separated into two parts. First the first 70
# columns are summarized into a table and then the rest. After exporting them,
# the tables can simply be combined in excel.
tab1 <- CreateTableOne(vars = cols, data = df_named, factorVars = fac.col, 
                       strata = "CM", addOverall = TRUE)
#tab2 <- CreateTableOne(vars = cols[71:length(cols)], data = df_named, 
#                       factorVars = fac.col, strata = "CM", addOverall = TRUE)
ex1 <- print(tab1, showAllLevels = FALSE, missing = TRUE, nonnormal = nn.col,
             formatOptions = list(big.mark = ","), quote = FALSE, 
             noSpaces = TRUE, printToggle = FALSE)
#ex2 <- print(tab2, showAllLevels = FALSE, missing = TRUE, nonnormal = nn.col,
#             formatOptions = list(big.mark = ","), quote = FALSE, 
#             noSpaces = TRUE, printToggle = FALSE)

write.csv(ex1, paste0(out_path, "MCM_Table1_norm.csv"))
#write.csv(ex2, paste0(out_path, "CM_Table1b.csv"))

# And last but not least, the final dataframe used to create the baseline-
# table is exported as well.
write.table(df_named, paste0(out_path, "CM_baseline_all.tsv"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
saveRDS(all, paste0(out_path, "All_list_baseline.rds"))
rm(desc, df, df_named, ex, gene, nas, rep, summ, tab, tes, cols, fac.col, i, 
   mean, cmr, q, met, ecg, bp, diag, death, rem, ocol)

