## Script information ---------------------------
##
## Script name: UKB_pheno_clean.R
##
## Purpose of script: Clean up a phenotype-summ of the UKB.
##
## Author: M. van Vugt
## Originally written by Mimount, but adapted for speed and readability.
##
## Date Created: 2021-05-04
##
## Copyright (c) M. van Vugt, 2021
## Email: m.vanvugt-2@umcutrecht.nl
##

## Options ---------------------------

options(scipen = 6, digits = 4) # view outputs in non-scientific notation
rm(list=ls())


## Loading arguments --------------------------

# Arguments expected:
#     #1 -- Path to and name of the phenotype file, 
#           for example /hpc/dhl_ec/mvanvugt/UKBB/Project1_ukb_phenotypes.tab
#     #2 -- Directory to the helper files, 
#           for example /hpc/dhl_ec/mvanvugt/Software/UKB-pipeline-Utrecht
#     #3 -- Output directory, for example /hpc/dhl_ec/mvanvugt/UKBB
#     #4 -- Prefix of the output files

args = commandArgs(trailingOnly = TRUE)
pheno = args[1] # "Data/temp/HCM_full.txt"
files = args[2] # "Data/raw"
output = args[3] # "Data/processed"
prefix = args[4] # "MCM"


## Loading packages ---------------------------

library(dplyr)
library(data.table)


# Loading data ------------------------------------------------------------

message("Loading data")
df <- data.table(readRDS(pheno)) 
df$CM <- as.factor(df$CM)
cols <- colnames(df)

setDT(df)
backup <- df

lof <- read.table(paste0(files, "/MCM_LoF_genes.txt"), 
                  header = TRUE, sep = "\t") 
lof <- lof %>% slice_tail(n = nrow(lof) - 2)
# lof[names(lof)[2:ncol(lof)]] <- lapply(lof[names(lof)[2:ncol(lof)]], as.factor)
gen <- read.table(paste0(files, "/CM_genes.txt"), col.names = c("Gene", "CM"))


# Start renaming ----------------------------------------------------------

dfs <- list()
# Summarizing the LoF-genes per CM
for (cm in levels(df$CM)) {
  
  ids <- df %>% filter(CM == cm) %>% select(f.eid)
  ids$f.eid <- as.character(ids$f.eid)
  gnc <- merge(ids, lof, all.x = TRUE)
  
  if (cm != "Controls") {
    
    cgn <- gen %>% filter(CM == cm)
    cgn <- cgn$Gene
    
    tmp <- gnc %>% select(f.eid, any_of(cgn)) 
    tmp[tmp == 0] <- NA
    tmp <- tmp %>% select_if(colSums(!is.na(.)) > 0)
    tmp[is.na(tmp)] <- 0
    
    vars <- df %>% filter(CM == cm) %>% 
      select(f.eid, ends_with("variant")) %>%
      as.data.frame()
    vars[names(vars)[1:ncol(vars)]] <- lapply(vars[names(vars)[1:ncol(vars)]], as.character)
    vars <- as.data.table(vars)
    vdf <- data.frame()
    for (line in 1:nrow(vars)) {
      var <- vars[line, f.eid]
      for (c in 2:ncol(vars)) {
        if (!is.na(vars[line, get(names(vars)[c])])) {
          var <- c(var, gsub("_variant", "", names(vars)[c]))
        } # End check value variant
      } # End iteration columns
    var <- c(var, rep(NA, (ncol(vars) - length(var))))
    vdf <- rbind(vdf, var)
    vdf[names(vdf)[1:ncol(vdf)]] <- lapply(vdf[names(vdf)[1:ncol(vdf)]], as.character)
    }
    vdf[names(vdf)[2:ncol(vdf)]] <- lapply(vdf[names(vdf)[2:ncol(vdf)]], as.factor)
    vdf <- vdf %>% select_if(colSums(!is.na(.)) > 0)
    names(vdf)[1] <- "f.eid"
    names(vdf)[2:ncol(vdf)] <- paste0("Gene_", 2:ncol(vdf)-1)
    
    dfs[[cm]] <- vdf
    dfs[[cm]]$CM <- cm
    
  } else {
    
    dfs[[cm]] <- ids
    dfs[[cm]]$CM <- cm
    tmp <- gnc %>% select_if(colSums(!is.na(.)) > 0)
    
  } # End checking for controls
  
  sdf <- data.frame()
  ddf <- data.frame()
  for (line in 1:nrow(tmp)) {
    single <- tmp[line, f.eid]
    double <- tmp[line, f.eid]
    for (c in 2:ncol(tmp)) {
      if (!is.na(tmp[line, get(names(tmp)[c])]) && tmp[line, get(names(tmp)[c])] == 1) {
        single <- c(single, names(tmp)[c])
      } else if (!is.na(tmp[line, get(names(tmp)[c])]) && tmp[line, get(names(tmp)[c])] == 2) {
        double <- c(double, names(tmp)[c])
      } # End check value LoF-gene
    } # End iteration columns
    single <- c(single, rep(NA, (ncol(tmp) - length(single))))
    sdf <- rbind(sdf, single)
    sdf[names(sdf)[1:ncol(sdf)]] <- lapply(sdf[names(sdf)[1:ncol(sdf)]], as.character)
    double <- c(double, rep(NA, (ncol(tmp) - length(double))))
    ddf <- rbind(ddf, double)
    ddf[names(ddf)[1:ncol(ddf)]] <- lapply(ddf[names(ddf)[1:ncol(ddf)]], as.character)
  }
  
  sdf[names(sdf)[2:ncol(sdf)]] <- lapply(sdf[names(sdf)[2:ncol(sdf)]], as.factor)
  sdf <- sdf %>% select_if(colSums(!is.na(.)) > 0)
  names(sdf)[1] <- "f.eid"
  if (ncol(sdf) > 1) {
    names(sdf)[2:ncol(sdf)] <- paste0("One_copy_LoF_", 2:ncol(sdf)-1)
  }
  
  ddf[names(ddf)[2:ncol(ddf)]] <- lapply(ddf[names(ddf)[2:ncol(ddf)]], as.factor)
  ddf <- ddf %>% select_if(colSums(!is.na(.)) > 0)
  names(ddf)[1] <- "f.eid"
  if (ncol(ddf) > 1) {
    names(ddf)[2:ncol(ddf)] <- paste0("Two_copy_LoF_", 2:ncol(ddf)-1)
  }
  
  new <- merge(sdf, ddf, by = "f.eid")
  dfs[[cm]] <- merge(dfs[[cm]], new, by = "f.eid")

} # End summarizing per CM


# To make one final dataframe, all columns should be the same and are therefore 
# created where necessary and reordered alphabetically
for (l in 1:(length(dfs) - 1)) {
  
  dfs[[l]] <- as.data.frame(dfs[[l]])
  col <- names(dfs[[l]])
  for (n in 2:length(dfs)) {
    dfs[[n]] <- as.data.frame(dfs[[n]])
    con <- names(dfs[[n]])
    if (length(con[!con %in% col]) > 0) {
      new <- con[!con %in% col]
      for (c in new) {
        dfs[[l]][c] <- NA
      } # End iteration new columns
    } # End check columns1
    if (length(col[!col %in% con]) > 0) {
      new <- col[!col %in% con]
      for (c in new) {
        dfs[[n]][c] <- NA
      } # End iteration new columns
    } # End check columns2
    dfs[[n]] <- dfs[[n]] %>% select(order(colnames(.)) )%>% 
      select(f.eid, everything())
    dfs[[n]][, 1] <- as.character(dfs[[n]][, 1])
    dfs[[n]][, 2:ncol(dfs[[n]])] <- lapply(dfs[[n]][, 2:ncol(dfs[[n]])], as.factor)
  } # End iteration next dfs
  dfs[[l]] <- dfs[[l]] %>% select(order(colnames(.))) %>% 
    select(f.eid, everything())
  dfs[[l]][, 1] <- as.character(dfs[[l]][, 1])
  dfs[[l]][, 2:ncol(dfs[[l]])] <- lapply(dfs[[l]][, 2:ncol(dfs[[l]])], as.factor)
} # End iteration first dfs
df_gen <- do.call(bind_rows, dfs)


# Saving data -------------------------------------------------------------

message("Saving the data")
saveRDS(df_gen, paste0("data/temp/", prefix, "_gene_summary.rds"))
write.table(df_gen, file = paste0("data/temp/", prefix, "_gene_summary.tsv"), 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

message("Enjoy the results!")

