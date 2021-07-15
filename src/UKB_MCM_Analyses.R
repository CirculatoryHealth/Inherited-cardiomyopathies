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


## Loading functions ---------------------------

library(plyr)
source("~/functions.R") 


## Loading packages ---------------------------

library(data.table)
library(ggpubr)
library(viridis)


# Loading Data ------------------------------------------------------------



# Prevalence --------------------------------------------------------------

prev <- data.table()
for (i in 1:3) {
  df <- all[[i]]
  tmp <- perc_var(df, c("Gene_1", "Gene_2")) %>% 
    na.omit() %>%
    mutate(perc = count / 200643)
  if (names(all)[i] == "hcm") {
    dup <- tmp %>% filter(value == "MYH7")
    tmp <- tmp %>% filter(value != "MYH7")
    new <- c("MYH7", sum(dup$count), sum(dup$perc), "Gene")
    tmp <- rbind(tmp, new)
  }
  tmp$name <- toupper(names(all)[i])
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
write.table(prev, "Results/output/Prevalence.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE)
rm(i, df, tmp)


# Enrichment CVDs ---------------------------------------------------------

cvd <- all$hcm %>% 
  select(Sex, Ethnicity, Cardiovascular_Death, All_cause_mortality, Ever_Smoked, 
         ends_with("SR")) %>% 
  names()
test <- list()
pval <- data.frame(Phenotype = cvd)
for (i in 1:3) {
  df <- rbind(all[[i]], all[[4]]) %>% dplyr::select(any_of(cvd), CM)
  df$CM <- droplevels(df$CM)
  ps <- vector()
  for (x in 1:length(cvd)) {
    tmp <- df[, c(x, ncol(df))]
    test[[cvd[x]]] <- fisher.test(table(tmp), workspace = 1e9)
    ps <- c(ps, test[[cvd[x]]]$p.value)
    write.table(table(tmp), "Results/output/CVD_enrichment_XTab.tsv", 
                sep = "\t", quote = FALSE, append = TRUE,
                col.names = c(cvd[x], paste("Controls", toupper(names(all)[i]), 
                                            sep = "_")))
  }
  pval <- cbind(pval, ps)
  names(pval)[ncol(pval)] <- toupper(names(all)[i])
}
write.table(pval, "Results/output/CVD_enrichment_fisher.tsv", sep = "\t", 
            quote = FALSE, row.names = FALSE)
rm(test, pval, df, ps, i, tmp, x)


# Other statistics --------------------------------------------------------

con <- all$hcm %>% 
  select(!c(f.eid, Ethnicity, Death, CM, any_of(cvd))) %>% 
  names()
test <- list()
pval <- data.frame(Phenotype = con)
for (i in 1:3) {
  df <- rbind(all[[i]], all[[4]]) %>% dplyr::select(any_of(con), CM)
  df$CM <- droplevels(df$CM)
  ps <- vector()
  sta <- vector()
  for (x in 1:(length(df)-1)) {
    tmp <- df[, x][df$CM == toupper(names(all)[i])]
    ref <- df[, x][df$CM == "CONTROLS"]
    test[[names(df)[x]]] <- wilcox.test(tmp, ref)
    ps <- c(ps, test[[names(df)[x]]]$p.value)
    sta <- c(sta, test[[names(df)[x]]]$statistic)
  }
  pval <- cbind(pval, ps)
  names(pval)[ncol(pval)] <- paste0(toupper(names(all)[i]), "_P")
  pval <- cbind(pval, sta)
  names(pval)[ncol(pval)] <- paste0(toupper(names(all)[i]), "_statistic")
}
write.table(pval, "Results/output/Differences_Continuous_Wilcox.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)
rm(con, cvd, test, pval, df, ps, i, x, tmp, ref)


# Filtering CM/HF diagnosed -----------------------------------------------

cms <- c("Heart_Failure_SR", "Cardiomyopathy_SR", "DCM_SR", "HCM_SR", 
         "Ventricular_arrhythmias_SR", "Cardiac_arrest_SR")
cmr <- c("RVEDVi", "RVESVi", "RVSV", "RVEF", "RVPER", "RVPFR", "RVPAFR",
         "LVEDVi", "LVESVi", "LVSV", "LVEF", "LVPER", "LVPFR", "LVPAFR",
         "LVEDMi", "LVEDM/LVEDV", "LV_RV_EDV")
ecg <- c("ECG_heart_rate.0_mean", "P_duration", "P_axis", "PQ_interval", 
         "QRS_duration", "R_axis", "QTC_interval", "T_axis")

lsi <- list()
for (i in names(all)) {
  df <- all[[i]] %>% dplyr::dplyr::select(f.eid, all_of(cms))
  sick <- data.frame()
  for (x in 2:ncol(df)) {
    tmp <- data.frame(df[, 1][df[, x] == "Yes"])
    names(tmp) <- "f.eid" 
    if (nrow(tmp) > 0) {
      tmp$Diagnose <- names(df)[x]
    }
    sick <- rbind(tmp, sick)
  }
  sick <- merge(all[[i]], sick, by = "f.eid", all.x = TRUE)
  sick$Diagnose[is.na(sick$Diagnose)] <- "No"
  sick$Diagnose <- as.factor(sick$Diagnose)
  lsi[[i]] <- unique(sick)
  
  df <- all[[i]] %>% dplyr::dplyr::select(f.eid, any_of(cmr))
  df$MRI <- rowMeans(df[, 2:ncol(df)], na.rm = TRUE)
  df$MRI[!is.na(df$MRI)] <- "Yes"
  df$MRI[df$MRI != "Yes"] <- "No"
  df$MRI <- as.factor(df$MRI)
  lsi[[i]] <- merge(lsi[[i]], df[, c(1, ncol(df))])
  
  df <- all[[i]] %>% dplyr::dplyr::select(f.eid, any_of(ecg))
  df$ECG <- rowMeans(df[, 2:ncol(df)], na.rm = TRUE)
  df$ECG[!is.na(df$ECG)] <- "Yes"
  df$ECG[df$ECG != "Yes"] <- "No"
  df$ECG <- as.factor(df$ECG)
  lsi[[i]] <- merge(lsi[[i]], df[, c(1, ncol(df))])
}

filt <- data.frame()
for (i in 1:length(lsi)) {
  filt <- rbind(filt, lsi[[i]])
}
saveRDS(filt, file = "Data/processed/All_Diagnose_ECG_MRI.rds")
