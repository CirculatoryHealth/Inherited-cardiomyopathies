#!/usr/bin/env Rscript
##
## Script information ----------------------------------------------------------
##
## Script name: Extract_IID_WES_UKB.R
##
## Purpose of script: Clean up and summarize the variant-based information
##                    to gene-based information
##
## Author: M. van Vugt
##
## Date Created: 2022/04/14
##
## Copyright (c) M. van Vugt, 2022
## Email: m.vanvugt-2@umcutrecht.nl
##

## Options ---------------------------------------------------------------------

options(scipen = 6, digits = 4) # view outputs in non-scientific notation

## Loading packages ------------------------------------------------------------

suppressMessages(library(dplyr))
suppressMessages(library(data.table))


## Loading arguments -----------------------------------------------------------

# Arguments expected:
#     #1 -- Path to and name of the file with IDs, genes and CMR data,
#           for example /hpc/dhl_ec/mvanvugt/UKBB/data/temp/ACM_temp/ACM_SNPs_MRI.txt
#     #2 -- Path to and name of the file with all genes of this CM,
#           for example /hpc/dhl_ec/mvanvugt/UKBB/data/temp/ACM_ExtractIID_genes.txt
#     #3 -- Path to and name of the file with included genes,
#           for example /hpc/dhl_ec/mvanvugt/UKBB/data/temp/ACM_temp/ACM_genes.txt
#     #4 -- Path to and name of the output file,
#           for example /hpc/dhl_ec/mvanvugt/UKBB/data/temp/ACM_temp/ACM_IIDs_genes_variants.txt

args   = commandArgs(trailingOnly = TRUE)
input  = args[1]
genes  = args[2]
vector = args[3]
output = args[4]


## Loading files ---------------------------------------------------------------

cat("Loading the data\n")
df <- data.table(read.table(input, header = TRUE, stringsAsFactors = FALSE))
gen <- read.table(genes, header = FALSE, stringsAsFactors = FALSE)
gen <- na.omit(gen)
vec <- read.table(vector, header = FALSE, stringsAsFactors=FALSE)[, 1]

# Isolate ID and information on CMR availability
ind <- data.frame(ID = df$ID, LV = df$LV, RV = df$RV)


## Summarizing on gene-level ---------------------------------------------------

cat("Starting the loop of genes\n")
for (i in 1:nrow(gen)) {
    x <- gen[i,]
    if (x %in% vec) {

        # if individual carries variant in gene, gets 1 in column for gene x
        cat(paste0("\nWorking on ", x, ", gene ", i, " / ", nrow(gen), "\n"))
        lok <- names(df)[grep(x, df)]
        ind[, paste0(x, "_variant")] <- as.integer(apply(df[, ..lok], 1, function(r) any(r %in% c(grep("1", r, value = T)))))

    } else {
        cat(paste0("\n Gene ", x, " not found among the predefined genes, going to the next\n"))
        next
    }
}

# Add column with total number of variants
ind$Total_variants <- rowSums(ind[, 4:ncol(ind)], na.rm = TRUE)
ind <- slice_head(ind, n = (nrow(ind)-1))
# Remove individuals without any variants
ind <- ind %>% filter(Total_variants > 0)


## Saving data -----------------------------------------------------------------

cat("\nWriting output")
write.table(ind, output, col.names = TRUE, row.names = FALSE, quote = FALSE)


