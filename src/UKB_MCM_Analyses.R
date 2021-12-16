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
library(rstatix)
library(ggplot2)
library(RColorBrewer)


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


# Making extra groups -----------------------------------------------------

# For subanalyses we make also groups of the overlapping individuals between CMs
# and two groups of diagnosed and non-diagnosed people (with and without controls)
# Comment this section out if only CMs should be analysed
acm <- subset(df, CM == "ACM")
ddcm <- df %>% filter(CM == "DCM")
hcm <- subset(df, CM == "HCM")
dcm <- df %>% filter(CM == "DCM") %>% select(f.eid)

ad <- merge(acm, dcm, all = FALSE)
ac <- anti_join(acm, ad)
d <- ad %>% select(f.eid)
d <- merge(d, ddcm, all.x = TRUE)
da <- anti_join(ddcm, d)
ad$CM <- "AD"
ac$CM <- "ACM_AD"
da$CM <- "DCM_AD"

hd <- merge(hcm, dcm, all = FALSE)
hc <- anti_join(hcm, hd)
d <- hd %>% select(f.eid)
d <- merge(d, ddcm, all.x = TRUE)
dh <- anti_join(ddcm, d)
hd$CM <- "HD"
hc$CM <- "HCM_HD"
dh$CM <- "DCM_HD"

new <- rbind(df, ad, ac, da, hd, hc, dh)

dia <- df %>% filter(Pheno == "Diagnosed") 
non <- df %>% filter(Pheno == "Non-Diagnosed") 
dia$CM <- "Diagnosed_all"
non$CM <- "Non-Diagnosed_all"
new <- rbind(new, dia, non)
dia <- df %>% filter(Pheno == "Diagnosed") %>% filter(CM != "Controls")
non <- df %>% filter(Pheno == "Non-Diagnosed") %>% filter(CM != "Controls")
dia$CM <- "Diagnosed_CM"
non$CM <- "Non-Diagnosed_CM"
new <- rbind(new, dia, non)
df <- new
df$CM <- as.factor(df$CM)
rm(acm, ddcm, hcm, dcm, ad, ac, d, da, hd, hc, dh, dia, non, new)


# Testing differences -----------------------------------------------------

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
cvd <- df %>% select(ends_with("sum")) %>% names()
dia <- c("Heart_Failure_sum", "Cardiomyopathy_sum", "HCM_sum", "DCM_sum",
         "Chronic_ischaemic_heart_disease_sum")
cols <- c("BMI", met, bp, ecg, cmr, cvd)

# Make list of continuous variables
cont <- cols[!cols %in% cvd]

# Prepare results dataframe
pval <- data.frame(Phenotype = rep("Outcomes", 6), Estimate = rep(NA, 6), 
                   LCI = rep(NA, 6), UCI = rep(NA, 6), p = rep(NA, 6),
                   test = rep(NA, 6), CM = c("ACM", "ACM_NonDiag", "DCM", 
                                             "DCM_NonDiag", "HCM", "HCM_NonDiag"))
prim <- c("ACM", "DCM", "HCM")
sec <- c("ACM_AD", "DCM_AD", "DCM_HD", "HCM_HD", "Diagnosed_CM")

message("Test for differences")
for (cm in c(prim, sec)) {
  if (cm %in% prim) {
    # Make df of only this CM and controls
    f <- df %>% select(any_of(cols), CM, Pheno) %>% 
      filter(CM %in% c(cm,  "Controls"))
    f$CM <- droplevels(f$CM)
    
    # Make df of only non-diagnosed subjects
    fh <- f %>% filter(Pheno == "Non-Diagnosed")
    fh$CM <- droplevels(fh$CM)
    
  } else if (cm %in% c("ACM_AD", "DCM_AD", "DCM_HD", "HCM_HD")) {
    # Make df of overlapping and non-overlapping subjects
    f <- df %>% select(any_of(cols), CM) %>% 
      filter(CM %in% c(cm,  unlist(strsplit(cm, "_"))[2]))
    f$CM <- droplevels(f$CM)
  } else {
    # Make df of diagnosed and non-diagnosed subjects
    f <- df %>% select(any_of(cols), CM) %>% 
      filter(CM %in% c(cm,  paste0("Non-", cm)))
    f$CM <- droplevels(f$CM)
  }# End preparing dfs based on CM
  
  for (x in cols) {
    if (x %in% cvd && cm %in% prim) {
      tryCatch( {
        # Perform Fisher's exact test on outcome x
        tmp <- f %>% dplyr::select(any_of(x), CM)
        test <- fisher.test(table(tmp), workspace = 1e9)
        # Save OR and confidence interval
        or <- test$estimate
        ci <- test$conf.int
        # Make OR for ACM in same direction as for DCM and HCM
        if (cm == "ACM") {
          or <- 1 / or
          ci <- 1 / ci
          ci <- ci[c(2,1)]
        } 
        # Save results in dataframe
        ps <- c(gsub("_sum", "", x), or, ci, test$p.value, "Fisher", cm)
        pval <- rbind(pval, ps)
        
      }, error = function(e) {
        
        # If the test fails, replace OR, CI and p-value with NA
        ps <- c(gsub("_sum", "", x), NA, NA, NA, "Fisher", cm)
        pval <- rbind(pval, ps)
        
      }) # End try to perform test
      
      if (!x %in% dia) {
        tryCatch( {
          # Perform Fisher's exact test on outcome x for undiagnosed subjects
          tmp <- fh %>% dplyr::select(any_of(x), CM)
          test <- fisher.test(table(tmp), workspace = 1e9)
          # Save OR and confidence interval
          or <- test$estimate
          ci <- test$conf.int
          # Make OR for ACM in same direction as for DCM and HCM
          if (cm == "ACM") {
            or <- 1 / or
            ci <- 1 / ci
            ci <- ci[c(2,1)]
          } 
          # Save results in dataframe
          ps <- c(gsub("_sum", "", x), or, ci, test$p.value, "Fisher", paste0(cm, "_NonDiag"))
          pval <- rbind(pval, ps)
          
        }, error = function(e) {
          
          # If the test fails, replace OR, CI and p-value with NA
          ps <- c(gsub("_sum", "", x), NA, NA, NA, "Fisher", paste0(cm, "_NonDiag"))
          pval <- rbind(pval, ps)
          
        }) # End try to perform test
        
      } # Only variables not used in diagnosis test on diagnosed groups
      
    } else if (x %in% cvd && cm %in% sec) {
      if (x %in% dia && cm %in% c("ACM_AD", "DCM_AD", "DCM_HD", "HCM_HD")) {
        tryCatch( {
          # Perform Fisher's exact test on outcome x
          tmp <- f %>% dplyr::select(any_of(x), CM)
          test <- fisher.test(table(tmp), workspace = 1e9)
          # Save OR and confidence interval
          or <- 1/ test$estimate
          ci <- 1/ test$conf.int
          ci <- ci[c(2,1)]
          # Save results in dataframe
          ps <- c(gsub("_sum", "", x), or, ci, test$p.value, "Fisher", cm)
          pval <- rbind(pval, ps)
          
        }, error = function(e) {
          
          # If the test fails, replace OR, CI and p-value with NA
          ps <- c(gsub("_sum", "", x), NA, NA, NA, "Fisher", cm)
          pval <- rbind(pval, ps)
          
        }) # End try to perform test
      } else if(!x %in% dia && cm %in% sec) {
        tryCatch( {
          # Perform Fisher's exact test on outcome x
          tmp <- f %>% dplyr::select(any_of(x), CM)
          test <- fisher.test(table(tmp), workspace = 1e9)
          # Save OR and confidence interval
          or <- 1/ test$estimate
          ci <- 1/ test$conf.int
          ci <- ci[c(2,1)]
          # Save results in dataframe
          ps <- c(gsub("_sum", "", x), or, ci, test$p.value, "Fisher", cm)
          pval <- rbind(pval, ps)
          
        }, error = function(e) {
          
          # If the test fails, replace OR, CI and p-value with NA
          ps <- c(gsub("_sum", "", x), NA, NA, NA, "Fisher", cm)
          pval <- rbind(pval, ps)
          
        }) # End try to perform test
      } # diagnosis testing only for overlapping
      
    } else if (x %in% cont) {
      
      if (!x %in% c(cmr, ecg) || cm %in% sec) {
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
          # Save results (WHAT CI???)
          vals <- c(x, eff$effsize, test$conf.int, test$p.value, "MWU", cm)
        }
        pval <- rbind(pval, vals)
      } # Don't test ECG and CMR differences for full groups
      
      if (cm %in% prim) {
        # See if test fails or not
        tried <- try(wilcox.test(get(x) ~ CM, data = fh), silent = TRUE)
        
        if (inherits(tried, "try-error")) {
          # If test fails, replace Estimate, CI and p-value with NA
          vals <- c(x, NA, NA, NA, NA, "MWU", paste0(cm, "_NonDiag"))
        } else {
          # Otherwise, perform test
          test <- wilcox.test(get(x) ~ CM, data = fh, conf.int = TRUE)
          # Calculate effect size 
          eff <- wilcox_effsize(fh, as.formula(paste0(x, " ~ CM")))
          # Save results (WHAT CI???)
          vals <- c(x, eff$effsize, test$conf.int, test$p.value, "MWU", paste0(cm, "_NonDiag"))
        }
        pval <- rbind(pval, vals)
      } # Only test for CMs
      
    } # End choose type of test (Fisher or MWU)
    
  } # End iteration outcomes
  
} # End iteration CMs

pval[, c("Phenotype", "test", "CM")] <- lapply(pval[, c("Phenotype", "test", "CM")], as.factor)
pval[, c("Estimate", "LCI", "UCI", "p")] <- lapply(pval[, c("Estimate", "LCI", "UCI", "p")], as.numeric)

write.table(pval, "results/output/Difference_testing.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)
rm(cols, test, f, cm, x, cmr, ecg, met, bp, vals, eff, tmp, fh, 
   tried, cont, dia, ci, cvd, or, ps)


# Visualizing p-values ----------------------------------------------------

# Create df to visualise p-value distributions
pval$Group <- "All (n=818)"
prim <- pval %>% filter(CM %in% prim)
prim$Group <- "Cardiomyopathies (n=102)"
cm <- pval %>% filter(CM %in% c("ACM_NonDiag", "DCM_NonDiag", "HCM_NonDiag"))
cm$Group <- "Undiagnosed cardiomyopathies (n=261)"
ov <- pval %>% filter(CM %in% c("ACM_AD", "DCM_AD", "DCM_HD", "HCM_HD"))
ov$Group <- "Overlap (n=368)"
dia <- pval %>% filter(CM %in% "Diagnosed_CM")
dia$Group <- "Diagnosed vs. undiagnosed (n=87)"
all <- rbind(pval, prim, cm, ov, dia)
all$Group <- as.factor(all$Group)

# Prepare p-value dataframe for annotation of histograms
ks <- data.frame(Group = levels(all$Group), 
                 pvalue = rep(NA, length(levels(all$Group))), 
                 p = rep(.9, length(levels(all$Group))),
                 label = rep(NA, length(levels(all$Group))))
ks$pvalue[ks$Group == "All (n=818)"] <- ks.test(pval$p, "punif", exact = TRUE)$p.value
ks$pvalue[ks$Group == "Cardiomyopathies (n=102)"] <- ks.test(prim$p, "punif", exact = TRUE)$p.value
ks$pvalue[ks$Group == "Undiagnosed cardiomyopathies (n=261)"] <- ks.test(cm$p, "punif", exact = TRUE)$p.value
ks$pvalue[ks$Group == "Overlap (n=368)"] <- ks.test(ov$p, "punif", exact = TRUE)$p.value
ks$pvalue[ks$Group == "Diagnosed vs. undiagnosed (n=87)"] <- ks.test(dia$p, "punif", exact = TRUE)$p.value
ks$label[ks$pvalue > 0.05] <- "p > 0.05"
ks$label[ks$pvalue >= .01 & ks$pvalue <= 0.05] <- paste0("p = ", round(ks$pvalue[ks$pvalue >= .01 & ks$pvalue <= 0.05], 2))
ks$label[ks$pvalue < .01] <- paste0("p = ", formatC(ks$pvalue[ks$pvalue < .01], format = "e", digits = 1))
# ks$label[ks$pvalue >= .01] <- round(ks$pvalue[ks$pvalue >= .01], 2)
# ks$label[ks$pvalue < .01] <- pretty10exp(ks$pvalue[ks$pvalue < .01], digits = 2)

# Visualization
ggplot(all, aes(x = p)) +
  geom_histogram(aes(y = ..density..), color = "black", fill = "grey", size = .2) +
  geom_hline(yintercept = 1, size = .2) +
  facet_wrap(~Group) +
  geom_text(data = ks, y = 4.9, label = ks$label, size = 1) +
  labs(x = "p-value", y = "Density") +
  my_theme() 
ggsave("results/figures/Multiple_testing_distribution.svg", width = 8 * pix, height =  5  * pix)

 # See 
cum <- NULL
for (x in levels(pval$Phenotype)) {
  new <- data.frame(Phenotype = x, Comparison = "CM-Controls", Ratio = NA)
  new$Ratio <- pval %>% filter(CM %in% c("ACM", "DCM", "HCM")) %>% 
    filter(Phenotype == x) %>% filter(p < .05) %>% nrow / 3
  cum <- rbind(cum, new)
}
cum$Phenotype <- gsub("_", " ", cum$Phenotype)
cum %>% filter(Ratio > 0) %>%
  ggplot(aes(x = reorder(Phenotype, -Ratio), y = Ratio, fill = Phenotype)) +
  geom_bar(stat = "identity", color = "black", size = .2) +
  labs(x = "Outcome", y = "Ratio p-values < 0.05",
       title = "Ratio significant CMs per phenotype") +
  scale_fill_manual(values = c("#FFD167", "#cff27e", "#4cbd97", "#e8e1ef", 
                               "#30362f", "#168ab2", "#1c1f33", "#161925", 
                               "2176ae", "#ef476f", "#881600", "#666370", 
                               "#4b2142", "#575d90", "#c16200", "red"),
                    labels = gsub("_", " ", levels(cum$Phenotype))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 1/3), 
                     labels = c("0", "1/3", "2/3", "3/3")) +
  pc_theme() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none")
ggsave("results/figures/Ratio_significance_CMs.svg", width = 6 * pix, height = 4 * pix)


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
