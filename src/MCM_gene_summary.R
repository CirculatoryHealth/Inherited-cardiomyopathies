#!/usr/bin/env Rscript
##
## Script information ----------------------------------------------------------
##
## Script name: MCM_gene_summary.R
##
## Purpose of script: Add gene-columns containing the gene in which 
##                    individual carries variant
##
## Author: M. van Vugt
##
## Date Created: 2021/05/04
##
## Copyright (c) M. van Vugt, 2022
## Email: m.vanvugt-2@umcutrecht.nl
##

## Options ---------------------------------------------------------------------

options(scipen = 6, digits = 4) # view outputs in non-scientific notation
rm(list=ls())


## Loading packages ------------------------------------------------------------

suppressMessages(library(dplyr))
suppressMessages(library(data.table))


## Loading arguments -----------------------------------------------------------

# Arguments expected:
#     #1 -- Path to and name of the phenotype file,
#           for example /hpc/dhl_ec/mvanvugt/UKBB/data/raw/MCM_raw_full.rds
#     #2 -- Output directory, 
#           for example /hpc/dhl_ec/mvanvugt/UKBB/data/processed
#     #3 -- Prefix of the output files, for example MCM

args   = commandArgs(trailingOnly = TRUE)
pheno  = args[1] 
output = args[2] 
prefix = args[3]


## Loading files ---------------------------------------------------------------

message(paste0("Loading data from ", pheno))
df <- data.table(readRDS(pheno))
df$CM <- as.factor(df$CM)
cols <- colnames(df)

gen <- read.table(file.path(output, "CM_genes.txt"), col.names = c("Gene", "CM"))


## Start renaming --------------------------------------------------------------

message("Start summarizing")
d <- df %>% select(f.eid, CM, ends_with("variant")) %>% data.frame()
vars <- data.frame()

for (line in 1:nrow(d)) {
    row <- d[line, 1]

    for (col in 3:ncol(d)) {
        # If individual carries varian, add column with name of gene
        if (!is.na(d[line, col]) && d[line, col] == 1) {
            row <- c(row, gsub("_variant", "", names(d)[col]))
        }

    } # end iteration variants
    # Fill up with NA where necessary
    row <- c(row, rep(NA, (ncol(d) - length(row) - 2)))
    # Add to df
    vars <- rbind(vars, row)
    vars[names(vars)[1:ncol(vars)]] <- lapply(vars[names(vars)[1:ncol(vars)]], as.character)

} # End iteration rows

# Select only columns that contain information
vars <- vars %>% select_if(colSums(!is.na(.)) > 0)
# Rename the columns
names(vars)[1] <- "f.eid"
names(vars)[2:ncol(vars)] <- paste0("Gene_", 2:ncol(vars)-1)


## Saving data -----------------------------------------------------------------

message(paste0("Saving data in ", file.path(output, prefix), "_gene_summary.tsv"))
write.table(vars, paste0(file.path(output, prefix), "_gene_summary.tsv"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

