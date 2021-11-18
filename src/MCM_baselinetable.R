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

# Uncomment if needed: Selecting only certain genes for HCM
# gen <- c("MYH7", "MYBPC3", "MYL2", "MYL3", "ACTC1", "TNNI3", "TNNT2", "TPM1",
#          "CSRP3", "TNNC1", "JPH2"  )
# hcm <- df %>% filter(CM == "HCM")
# df <- df %>% filter(CM != "HCM")
# hcm <- df %>% filter(!Gene_1 %in% gen)
# df <- rbind(df, hcm)


# CMR parameters ----------------------------------------------------------

# n <- grep("LVEDM.LVEDV", names(df))
# names(df)[n] <- "LVMVR"
n <- grep("LVEDV.RVEDV", names(df))
names(df)[n] <- "LVEDV/RVEDV"
rm(n)
# df$Total_MET_minutes_per_week <- df$MET_minutes_per_week_for_moderate_activity.0.0 + df$MET_minutes_per_week_for_vigorous_activity.0.0 + df$MET_minutes_per_week_for_walking.0.0


# Making Baselinetable ----------------------------------------------------

# The vectors are created with column names which we want to include in the
# baseline-table. If a column name is not represented in the vector cols, it
# will not be included in the baseline-table, so this should be changed here.
cmr <- df %>% select(starts_with("RV"), starts_with("LV"), contains("Wall_thickness"), 
                     contains("Ecc", ignore.case = FALSE), 
                     contains("Ell", ignore.case = FALSE)) %>% names()
ecg <- c("ECG_heart_rate.0_mean", "P_duration", "P_axis.2.0",
         "PQ_interval.2.0", "QRS_duration", "R_axis.2.0",
         "QTC_interval.2.0", "T_axis.2.0")
met <- df %>% select(starts_with("MET")) %>% names()
bp <- df %>% select(Total_Cholesterol, HDL, LDL,
                    contains("blood_pressure_mean")) %>% names()
diag <- df %>% select(ends_with("sum")) %>% names()
cols <- c("Sex", "Age_when_attended_assessment_centre.0.0", "Ethnicity",
          "BMI", met, bp, diag, ecg, cmr, "CM", "Total_MET_minutes_per_week",
          "ECG", "CMR")
rm(cmr, ecg, met, bp, diag)

# A new df (df) is created, with all columns listed in cols
message("Selecting variables for baseline table")
df <- df %>% dplyr::select(f.eid, any_of(cols), ends_with("FirstDate"))


# Now we'll create vectors with column names for the baseline-table.
# fac.col contains all column names that should be regarded as factors
# nn.col contains all columns that should be regarded as nonnormal (for now
# these are all columns except the factors, but can be changed of course).
fac.col <- c("Sex", "Ethnicity", diag, "CM", "CMR", "ECG")
df <- as.data.frame(df)
df[fac.col] <- lapply(df[fac.col], as.factor)
nn.col <- cols[!cols %in% fac.col]
df$Sex <- factor(df$Sex, levels = c("Male", "Female"))
df$HFCM[df$Heart_Failure_sum == "Yes" | df$Cardiomyopathy_sum == "Yes" | df$HCM_sum == "Yes" | df$DCM_sum == "Yes"] <- "Yes"
df$HFCM[df$Heart_Failure_sum == "No" & df$Cardiomyopathy_sum == "No" & df$HCM_sum == "No" & df$HCM_sum == "No"] <- "No"

df$Pheno[df$Chronic_ischaemic_heart_disease_sum == "No" & df$HFCM == "Yes"] <- "Diagnosed"
df$Pheno[df$Chronic_ischaemic_heart_disease_sum == "Yes" & df$HFCM == "Yes"] <- "Non-Diagnosed"
df$Pheno[df$HFCM == "No"] <- "Non-Diagnosed"
df$Pheno <- as.factor(df$Pheno)

# ocol <- vector()
# norm <- list()
# for (n in nn.col) {
#   vec <- na.omit(df[n])
#   norm[[n]] <- list()
#   norm[[n]][["dens.plot"]] <- ggdensity(vec[,1])
#   norm[[n]][["qq.plot"]] <- ggqqplot(vec[,1])
#   norm[[n]][["test"]] <- ks.test(vec[,1], "pnorm")
# 
#   if (norm[[n]][["test"]]$p.value < 0.05) {
#     ocol <- c(ocol, n)
#   } # Finish if-loop adding nonnormal variables
# } # Finish for-loop iterating over non-factor columns
# 
# # After manual inspection, some variables are removed anyways
# rem <- c("Total_Cholesterol", "LDL", "Systolic_blood_pressure_mean",
#          "Diastolic_blood_pressure_mean", "ECG_heart_rate.0_mean", "P_duration",
#          "RVEF", "RVPFR", "LVPER", "LVPFR")
# ocol <- ocol[!ocol %in% rem]
con <- df %>% select(any_of(nn.col))
p <- con %>% select_if(is.numeric) %>% tidyr::gather(cols, value) %>% 
  ggplot(aes(x = value)) + 
  geom_histogram(aes(y = ..density..)) +
  facet_wrap(.~cols, scales = "free") +
  labs(x = "Value", y = "Density", 
       title = "Distribution continuous data") +
  my_theme() + theme(text = element_text(size = 10))
ggsave("results/figures/Distribution_continuous_data.svg")
rm(con, p)

message("Creating the baseline tables")
tab1 <- CreateTableOne(vars = cols, data = df, factorVars = fac.col,
                       strata = "CM", addOverall = FALSE)
ex1 <- print(tab1, showAllLevels = FALSE, missing = TRUE, nonnormal = nn.col,
             formatOptions = list(big.mark = ","), quote = FALSE,
             noSpaces = TRUE, printToggle = FALSE)

acm <- subset(df, CM == "ACM")
dcm <- subset(df, CM == "DCM")
hcm <- subset(df, CM == "HCM")
con <- subset(df, CM == "Controls")

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
tabc <- CreateTableOne(vars = cols, data = con, factorVars = fac.col,
                       strata = "Pheno", addOverall = TRUE)
exc <- print(tabc, showAllLevels = FALSE, missing = TRUE, nonnormal = nn.col,
             formatOptions = list(big.mark = ","), quote = FALSE,
             noSpaces = TRUE, printToggle = FALSE)

message("Saving data")
write.csv(ex1, paste0(output, prefix, "_Table1.csv"))
write.csv(exa, paste0(output, "ACM_Table1.csv"))
write.csv(exd, paste0(output, "DCM_Table1.csv"))
write.csv(exh, paste0(output, "HCM_Table1.csv"))
write.csv(exc, paste0(output, "Controls_Table1.csv"))

# And last but not least, the final dataframe used to create the baseline-
# table is exported as well.
write.table(df, paste0(output, prefix, "_table1_data.tsv"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
saveRDS(df, paste0(output, prefix, "_table1_data.rds"))
