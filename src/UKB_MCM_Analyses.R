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
input = args[1] # "data/processed/MCM_final_pheno.txt"
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
library(tableone)


# Loading Data ------------------------------------------------------------

message("Loading data")
df <- data.table(read.delim(input, header = TRUE, sep = "\t", stringsAsFactors = FALSE))
# lof <- readRDS("data/processed/MCM_gene_summary.rds")
df$CM <- as.factor(df$CM)


# Prevalence --------------------------------------------------------------

message("Calculating prevalence")
prev <- data.table()
for (cm in levels(df$CM)) {
  if (cm != "Controls") {
    
    f <- subset(df, CM == cm)
    tmp <- perc_var(f, c("Gene")) %>%
      na.omit() %>%
      mutate(perc = count / 200643)
    # if (cm %in% c("HCM", "DCM")) {
    #   dup <- tmp %>% filter(value == "TNNT2")
    #   tmp <- tmp %>% filter(value != "TNNT2")
    #   new <- c("TNNT2", sum(dup$count), sum(dup$perc), "Gene")
    #   tmp <- rbind(tmp, new)
    # } # End ifelse loop adding genes that occur twice
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
cvd <- df %>% select(ends_with("sum")) %>% names()
test <- list()
pval <- data.frame(Phenotype = cvd)
for (cm in levels(df$CM)) {
  if (cm != "Controls") {
    
    f <- df %>% select(any_of(cvd), CM) %>% filter(CM %in% c(cm,  "Controls"))
    f$CM <- droplevels(f$CM)
    
    fd <- df %>%
      select(any_of(cvd), Pheno, CM) %>%
      filter(CM %in% c(cm,  "Controls")) %>%
      filter(Pheno == "Diagnosed")
    fd$CM <- droplevels(fd$CM)
    
    fh <- df %>%
      select(any_of(cvd), Pheno, CM) %>%
      filter(CM %in% c(cm,  "Controls")) %>%
      filter(Pheno == "Non-Diagnosed")
    fh$CM <- droplevels(fh$CM)
    
    pa <- data.frame()
    for (x in cvd) {
      ps <- vector()
      
      tmp <- f %>% dplyr::select(any_of(x), CM)
      test[[x]] <- fisher.test(table(tmp), workspace = 1e9)
      or <- test[[x]]$estimate
      ci <- test[[x]]$conf.int
      if (or < 1) {
        or <- 1 / or
        ci <- 1 / test[[x]]$conf.int
      } 
      ps <- c(ps, test[[x]]$p.value, or, ci)
      
      tryCatch( {
        
        tmp <- fd %>% dplyr::select(any_of(x), CM)
        test[[x]] <- fisher.test(table(tmp), workspace = 1e9)
        or <- test[[x]]$estimate
        ci <- test[[x]]$conf.int
        if (or < 1) {
          or <- 1 / or
          ci <- 1 / test[[x]]$conf.int
        } 
        ps <- c(ps, test[[x]]$p.value, or, ci)
        
        tmp <- fh %>% dplyr::select(any_of(x), CM)
        test[[x]] <- fisher.test(table(tmp), workspace = 1e9)
        or <- test[[x]]$estimate
        ci <- test[[x]]$conf.int
        if (or < 1) {
          or <- 1 / or
          ci <- 1 / test[[x]]$conf.int
        } 
        ps <- c(ps, test[[x]]$p.value, or, ci)
        
      }, error = function(e) {
        
        ps <- c(ps, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
        
      })
      
      pa <- rbind(pa, ps)
      write.table(table(tmp), "results/output/CVD_enrichment_XTab.tsv",
                  sep = "\t", quote = FALSE, append = TRUE,
                  col.names = c(x, paste("Controls", cm, sep = "_")))
    }
    pval <- cbind(pval, pa)
    names(pval)[(ncol(pval)-11):ncol(pval)] <- c(cm, paste0(cm, "_OR"), paste0(cm, "_LCI"), paste0(cm, "_UCI"),
                                                paste0(cm, "_diagnosed"), paste0(cm, "_diag_OR"), paste0(cm, "_diag_LCI"), paste0(cm, "_diag_UCI"),
                                                paste0(cm, "_nondiagnosed"), paste0(cm, "_nondiag_OR"), paste0(cm, "_nondiag_LCI"), paste0(cm, "_nondiag_UCI"))
    
  }
}
write.table(pval, "results/output/CVD_enrichment_fisher.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)
rm(test, pval, f, ps, cm, tmp, x, pa, fh, fd, cvd)


# Other statistics --------------------------------------------------------

# n <- grep("LVEDM.LVEDV", names(df))
# names(df)[n] <- "LVMVR"
n <- grep("LVEDV.RVEDV", names(df))
names(df)[n] <- "LVEDV/RVEDV"

n <- grep("wall_thickness", names(df), ignore.case = TRUE)
df <- as.data.frame(df)
for (col in n) {
  t1 <- quantile(df[, col], 0.25, na.rm = TRUE) - (IQR(df[, col], na.rm = TRUE) * 3)
  t2 <- quantile(df[, col], 0.75, na.rm = TRUE) + (IQR(df[, col], na.rm = TRUE) * 3)
  
  df[, col][df[, col] > t2 | df[, col] < t1 ] <- NA
}
rm(n, t1, t2, col)

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
cols <- c("Age_when_attended_assessment_centre.0.0", "BMI", met, bp, ecg, cmr)

nn.col <- c("Age_when_attended_assessment_centre.0.0", "BMI", ecg, met,  
            "peakEll4Ch", "TPKEcc", "peakEcc", "peakEll2Ch", "RVESV", "RVSV", 
            "Septal_wall_thickness", "Global_wall_thickness", "LVEDM", "LVEDMi", 
            "LVEDV", "LVEDV/RVEDV", "LVEDVi", "LVESV", "LVESVi", "LVMVR", 
            "LVPPAFR", "LVSV", "LVSVi", "LV_RV_ESV", "LVEF")
test <- list()
pval <- data.frame(Phenotype = cols)
for (cm in levels(df$CM)) {
  if (cm != "Controls") {
    
    f <- df %>% select(CM, Pheno, any_of(cols)) %>%
      filter(CM %in% c(cm,  "Controls")) %>% as.data.frame()
    f$CM <- droplevels(f$CM)
    
    fd <- f %>% filter(Pheno == "Diagnosed")
    fd$CM <- droplevels(fd$CM)
    
    fh <- f %>% filter(Pheno == "Non-Diagnosed")
    fh$CM <- droplevels(fh$CM)
    
    row <- data.frame()
    
    for (x in names(f)) {
      vals <- vector()
      if (!x %in% c("CM", "Pheno", nn.col)) {
        tmp <- subset(f, CM == cm)[, x]
        ref <- subset(f, CM == "Controls")[, x]
        tried <- try(t.test(tmp, ref), silent = TRUE)
        
        if(inherits(tried, "try-error")) {
          vals <- c(vals, NA, NA, "T")
        } else {
          test[[x]] <- t.test(tmp, ref)
          vals <- c(vals, test[[x]]$p.value)
          vals <- c(vals, test[[x]]$statistic, "T")
        }
        
        tmp <- subset(fd, CM == cm)[, x]
        ref <- subset(fd, CM == "Controls")[, x]
        tried <- try(t.test(tmp, ref), silent = TRUE)
        
        if(inherits(tried, "try-error")) {
          vals <- c(vals, NA, NA, "T")
        } else {
          test[[x]] <- t.test(tmp, ref)
          vals <- c(vals, test[[x]]$p.value)
          vals <- c(vals, test[[x]]$statistic, "T")
        }
        
        tmp <- subset(fh, CM == cm)[, x]
        ref <- subset(fh, CM == "Controls")[, x]
        tried <- try(t.test(tmp, ref), silent = TRUE)
        
        if(inherits(tried, "try-error")) {
          vals <- c(vals, NA, NA, "T")
        } else {
          test[[x]] <- t.test(tmp, ref)
          vals <- c(vals, test[[x]]$p.value)
          vals <- c(vals, test[[x]]$statistic, "T")
        }
        row <- rbind(row, vals)
      } else if (!x %in% c("CM", "Pheno")) {
        tmp <- subset(f, CM == cm)[, x]
        ref <- subset(f, CM == "Controls")[, x]
        tried <- try(wilcox.test(tmp, ref), silent = TRUE)
        
        if(inherits(tried, "try-error")) {
          vals <- c(vals, NA, NA, "MWU")
        } else {
          test[[x]] <- wilcox.test(tmp, ref)
          vals <- c(vals, test[[x]]$p.value)
          vals <- c(vals, test[[x]]$statistic, "MWU")
        }
        
        tmp <- subset(fd, CM == cm)[, x]
        ref <- subset(fd, CM == "Controls")[, x]
        tried <- try(wilcox.test(tmp, ref), silent = TRUE)
        
        if(inherits(tried, "try-error")) {
          vals <- c(vals, NA, NA, "MWU")
        } else {
          test[[x]] <- wilcox.test(tmp, ref)
          vals <- c(vals, test[[x]]$p.value)
          vals <- c(vals, test[[x]]$statistic, "MWU")
        }
        
        tmp <- subset(fh, CM == cm)[, x]
        ref <- subset(fh, CM == "Controls")[, x]
        tried <- try(wilcox.test(tmp, ref), silent = TRUE)
        
        if(inherits(tried, "try-error")) {
          vals <- c(vals, NA, NA, "MWU")
        } else {
          test[[x]] <- wilcox.test(tmp, ref)
          vals <- c(vals, test[[x]]$p.value)
          vals <- c(vals, test[[x]]$statistic, "MWU")
        }
        row <- rbind(row, vals)
        row[names(row)] <- lapply(row[names(row)], as.character)
      } # End check for continuous normal variable
    } # End iteration variables
    names(row) <- c(paste0(cm, "_P"), paste0(cm, "_statistic"), 
                    paste0(cm, "_test"), paste0(cm, "_Diag_P"), 
                    paste0(cm, "_Diag_statistic"), paste0(cm, "_Diag_test"),
                    paste0(cm, "_NonDiag_P"), paste0(cm, "_NonDiag_statistic"), 
                    paste0(cm, "_NonDiag_test"))
    pval <- cbind(pval, row)
    
  } # End check CM-group
} # End iteration CMs
write.table(pval, "results/output/Differences_Continuous.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
rm(cols, cvd, test, pval, f, ps, cm, x, tmp, ref, cmr, ecg, met, bp, sta, vals, 
   nn.col, fd, fh, row, tried)


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
