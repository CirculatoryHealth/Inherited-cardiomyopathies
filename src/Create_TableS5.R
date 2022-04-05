#!/usr/bin/env Rscript
##
## Script information ----------------------------------------------------------
##
## Script name: Create_TableS5.R
##
## Purpose of script: Calculate gene prevalences (Table S5)
##
## Author: M. van Vugt
##
## Date Created: 2021/06/16
##
## Copyright (c) M. van Vugt, 2022
## Email: m.vanvugt-2@umcutrecht.nl
##

## Options ---------------------------------------------------------------------

options(scipen = 6, digits = 4) # view outputs in non-scientific notation


## Loading packages ------------------------------------------------------------

suppressMessages(library(plyr))
source("src/functions.R")
suppressMessages(library(data.table))


## Loading arguments -----------------------------------------------------------

# Arguments expected:
#     #1 -- Path to and name of the phenotype file,
#           for example /hpc/dhl_ec/mvanvugt/UKB/data/processed/MCM_final_pheno.tsv
#     #2 -- Path to and name of the output file,
#           for example /hpc/dhl_ec/mvanvugt/UKB/results/output/TableS5.tsv

args = commandArgs(trailingOnly = TRUE)
input    = args[1]
output   = args[2]


## Loading files ---------------------------------------------------------------

message(paste0("Loading data from ", input))
df <- data.table(read.delim(input, header = TRUE, sep = "\t", 
                            stringsAsFactors = FALSE))
df$CM[df$CM == "ACM"] <- "ARVC"
df$CM <- as.factor(df$CM)


## Calculate gene prevalence ---------------------------------------------------

message("Calculating prevalence")
prev <- data.table()
# Iterate over the CMs, but not controls, because they are not carriers
for (cm in levels(df$CM)) {
    if (cm != "Controls") {

        # Subset and summarise the data per CM and gene
        f <- subset(df, CM == cm)
        tmp <- perc_var(f, c("Gene")) %>%
            na.omit() %>%
            mutate(perc = count / 200643)
        tmp$name <- cm
        # Clean up df
        names(tmp) <- c("Gene", "N", "Prevalence", "Cardiomyopathy")
        tmp$N <- as.numeric(tmp$N)
        tmp$Prevalence <- as.numeric(tmp$Prevalence)
        prev <- rbind(prev, tmp)
        tmp <- tmp %>%
            arrange(desc(Gene)) %>%
            mutate(Proportion = N / sum(tmp$N) * 100)
    }
}


## Save data -------------------------------------------------------------------

message(paste0("Saving data in ", output))
write.table(prev, output, sep = "\t", quote = FALSE, col.names = TRUE, 
            row.names = FALSE)
