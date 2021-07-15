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

d <- df %>% dplyr::dplyr::select(f.eid, any_of(ecg))
d$ECG <- rowMeans(d[, 2:ncol(d)], na.rm = TRUE)
d$ECG[!is.na(d$ECG)] <- "Yes"
d$ECG[d$ECG != "Yes"] <- "No"
d$ECG <- as.factor(d$ECG)
d <- d %>% select(f.eid, ECG)
df <- merge(df, unique(d), by = "f.eid", all.x = TRUE)

saveRDS(df, file = "data/processed/MCM_filter_ready_full.rds")
