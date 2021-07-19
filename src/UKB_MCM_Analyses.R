## Script information ---------------------------
##
## Script name: UKB_MCM_Analyses.R
##
## Purpose of script:
##
## Author: M. van Vugt
##
## Date Created: 2021-06-16
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
# output = args[3] # "results/output/"
# prefix = args[4] # "MCM"


## Loading functions ---------------------------

library(plyr)
source("src/functions.R") 


## Loading packages ---------------------------

library(data.table)
library(dplyr)
library(ggpubr)
library(viridis)


# Loading Data ------------------------------------------------------------

message("Loading data")
df <- data.table(readRDS(input))
lof <- readRDS("data/processed/MCM_gene_summary.rds")


# Prevalence --------------------------------------------------------------

message("Calculating prevalence")
prev <- data.table()
for (cm in levels(df$CM)) {
  if (cm != "Controls") {
    
    f <- subset(df, CM == cm)
    tmp <- perc_var(f, c("Gene_1", "Gene_2")) %>% 
      na.omit() %>%
      mutate(perc = count / 200643)
    if (cm == "HCM") {
      dup <- tmp %>% filter(value == "MYH7")
      tmp <- tmp %>% filter(value != "MYH7")
      new <- c("MYH7", sum(dup$count), sum(dup$perc), "Gene")
      tmp <- rbind(tmp, new)
    } else if (cm == "DCM") {
      dup <- tmp %>% filter(value == "TNNI3")
      tmp <- tmp %>% filter(value != "TNNI3")
      new <- c("TNNI3", sum(dup$count), sum(dup$perc), "Gene")
      tmp <- rbind(tmp, new)
    } # End ifelse loop adding genes that occur twice
    tmp$name <- cm
    names(tmp) <- c("Gene", "N", "Prevalence_WES_UKB", "CM")
    tmp$N <- as.numeric(tmp$N)
    tmp$Prevalence_WES_UKB <- as.numeric(tmp$Prevalence_WES_UKB)
    prev <- rbind(prev, tmp)
    tmp <- tmp %>%
      arrange(desc(Gene)) %>%
      mutate(prop = N / sum(tmp$N) * 100) %>%
      mutate(ypos = cumsum(prop) - 0.5 * prop)
    
    source("src/Piechart_genes.R")
    
  }
} 
write.table(prev, "results/output/Prevalence.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE)
rm(cm, f, tmp, new, dup, m, pie, prev)


# Enrichment CVDs ---------------------------------------------------------

message("Test for enrichment in certain phenotypes")
cvd <- df %>% select(Sex, Ethnicity, ends_with("sum")) %>% names()
test <- list()
pval <- data.frame(Phenotype = cvd)
for (cm in levels(df$CM)) {
  if (cm != "Controls") {
    
    f <- df %>% select(any_of(cvd), CM) %>% filter(CM %in% c(cm,  "Controls"))
    f$CM <- droplevels(f$CM)
    ps <- vector()
    for (x in cvd) {
      tmp <- f %>% dplyr::select(any_of(x), CM)
      test[[x]] <- fisher.test(table(tmp), workspace = 1e9)
      ps <- c(ps, test[[x]]$p.value)
      write.table(table(tmp), "results/output/CVD_enrichment_XTab.tsv", 
                  sep = "\t", quote = FALSE, append = TRUE,
                  col.names = c(x, paste("Controls", cm, sep = "_")))
    }
    pval <- cbind(pval, ps)
    names(pval)[ncol(pval)] <- cm
    
  }
}
write.table(pval, "results/output/CVD_enrichment_fisher.tsv", sep = "\t", 
            quote = FALSE, row.names = FALSE)
rm(test, pval, f, ps, cm, tmp, x)


# Other statistics --------------------------------------------------------

cmr <- df %>% select(starts_with("RV"), starts_with("LV")) %>% 
  select(!c(LV, RV)) %>% names()
ecg <- c("ECG_heart_rate.0_mean", "P_duration", "P_axis.2.0", 
         "PQ_interval.2.0", "QRS_duration", "R_axis.2.0", 
         "QTC_interval.2.0", "T_axis.2.0")
met <- df %>% select(starts_with("MET")) %>% names()
bp <- df %>% select(Total_Cholesterol, HDL, LDL,
                    contains("blood_pressure_mean")) %>% names()
cols <- c("Age_when_attended_assessment_centre.0.0", "BMI", met, bp, ecg, cmr)

test <- list()
pval <- data.frame(Phenotype = cols)
for (cm in levels(df$CM)) {
  if (cm != "Controls") {
    
    f <- df %>% select(CM, any_of(cols)) %>% 
      filter(CM %in% c(cm,  "Controls")) %>% as.data.frame()
    f$CM <- droplevels(f$CM)
    ps <- vector()
    sta <- vector()
    for (x in names(f)) {
      if (x != "CM") {
        tmp <- subset(f, CM == cm)[, x]
        ref <- subset(f, CM == "Controls")[, x]
        test[[x]] <- wilcox.test(tmp, ref)
        ps <- c(ps, test[[x]]$p.value)
        sta <- c(sta, test[[x]]$statistic)
      } # End check for continuous variable
    } # End iteration variables
    pval <- cbind(pval, ps)
    names(pval)[ncol(pval)] <- paste0(cm, "_P")
    pval <- cbind(pval, sta)
    names(pval)[ncol(pval)] <- paste0(cm, "_statistic")
    
  } # End check CM-group
} # End iteration CMs
write.table(pval, "results/output/Differences_Continuous_Wilcox.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)
rm(cols, cvd, test, pval, f, ps, cm, x, tmp, ref, cmr, ecg, met, bp, sta)


# Filtering CM/HF diagnosed -----------------------------------------------

cms <- c("Heart_Failure_sum", "Cardiomyopathy_sum", "DCM_sum", "HCM_sum", 
         "Ventricular_arrhythmias_sum", "Cardiac_arrest_sum")
ecg <- c("ECG_heart_rate.0_mean", "P_duration", "P_axis.2.0", 
         "PQ_interval.2.0", "QRS_duration", "R_axis.2.0", 
         "QTC_interval.2.0", "T_axis.2.0")

d <- df %>% select(f.eid, all_of(cms))
d[, "CVD_diagnosis"] <- as.integer(apply(d[, 2:ncol(d)], 1, function(r) any(r %in% c(grep("Yes", r, value = T)))))
d <- d %>% select(f.eid, CVD_diagnosis)
d$CVD_diagnosis <- ifelse(d$CVD_diagnosis >= 1, "Yes", "No")
df <- merge(df, unique(d), by = "f.eid", all.x = TRUE)

d <- df %>% select(f.eid, starts_with("LV"), starts_with("RV")) %>%
  select(!c(LV, RV))
d$MRI <- rowMeans(d[, 2:ncol(d)], na.rm = TRUE)
d$MRI[!is.na(d$MRI)] <- "Yes"
d$MRI[d$MRI != "Yes"] <- "No"
d$MRI <- as.factor(d$MRI)
d <- d %>% select(f.eid, MRI)
df <- merge(df, unique(d), by = "f.eid", all.x = TRUE)

d <- df %>% dplyr::select(f.eid, any_of(ecg))
d$ECG <- rowMeans(d[, 2:ncol(d)], na.rm = TRUE)
d$ECG[!is.na(d$ECG)] <- "Yes"
d$ECG[d$ECG != "Yes"] <- "No"
d$ECG <- as.factor(d$ECG)
d <- d %>% select(f.eid, ECG)
df <- merge(df, unique(d), by = "f.eid", all.x = TRUE)

saveRDS(df, file = "data/processed/MCM_filter_ready_full.rds")


# Loss-of-Function --------------------------------------------------------

lof <- data.table(read.delim("data/temp/MCM_LoF_genes.txt", sep = "\t", header = T, stringsAsFactors=FALSE))
gen <- readRDS("data/processed/MCM_gene_summary.rds")[, 1:2]

lof <- merge(gen, lof, by = "f.eid", all.x = TRUE)
new <- lof %>% select(f.eid, CM) %>% data.table

for (i in 3:ncol(lof)) {
  g <- names(lof)[i]
  sub <- data.frame(f.eid = subset(lof, lof[, i] == 1)[, 1])
  if (nrow(sub) > 0) {
    sub[, paste0(g, "_1")] <- "Yes"
    new <- merge(new, unique(sub), by = "f.eid", all.x = TRUE)
  }
  sub <- data.frame(f.eid = subset(lof, lof[, i] == 2)[, 1])
  if (nrow(sub) > 0) {
    sub[, paste0(g, "_2")] <- "Yes"
    new <- merge(new, unique(sub), by = "f.eid", all.x = TRUE)
  }
}
new <- data.frame(new)
new[is.na(new)] <- "No"
new[2:ncol(new)] <- lapply(new[2:ncol(new)], as.factor)
cols <- names(new)[3:ncol(new)]

tab1 <- CreateTableOne(vars = cols, data = new, factorVars = cols, 
                       strata = "CM", addOverall = FALSE)

ex1 <- print(tab1, showAllLevels = FALSE, formatOptions = list(big.mark = ","), 
             quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

write.csv(ex1, "results/output/LoF_Table.csv")

write.table(new, "data/temp/LoF_Table_data.tsv", row.names = FALSE, 
            col.names = TRUE, quote = FALSE, sep = "\t")
saveRDS(new, "data/temp/LoF_Table_data.rds")
