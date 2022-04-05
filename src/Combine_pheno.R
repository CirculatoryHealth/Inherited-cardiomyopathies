#!/usr/bin/env Rscript
##
## Script information ----------------------------------------------------------
##
## Script name: Combine_pheno.R
##
## Purpose of script: Combine the phenotype files of controls and all CMs
##
## Author: M. van Vugt
##
## Date Created: 2021/07/14
##
## Copyright (c) M. van Vugt, 2022
## Email: m.vanvugt-2@umcutrecht.nl
##

## Options ---------------------------------------------------------------------

options(scipen = 6, digits = 4) # view outputs in non-scientific notation


## Loading packages ------------------------------------------------------------

suppressMessages(library(dplyr))


## Loading arguments -----------------------------------------------------------

# Arguments expected:
#     #1 -- Path to the phenotype files,
#           for example /hpc/dhl_ec/mvanvugt/UKBB
#     #2 -- Suffix of the phenotype files, for example _full.tsv
#     #3 -- Output file name, for example MCM_raw_full

args   = commandArgs(trailingOnly = TRUE)
path   = args[1] 
suffix = args[2]
output = args[3]


## Loading files ---------------------------------------------------------------

message(paste0("Loading data from ", path, "/*", suffix))

dfs <- list()
for (cm in c("ACM", "DCM", "HCM", "Controls")) {

    dfs[[cm]] <- read.delim(paste0(file.path(path, cm), suffix), 
                            sep = "\t", header = T) %>%
                select(!starts_with("X"))
    dfs[[cm]]$CM <- cm
    dfs[[cm]]$f.eid <- as.character(dfs[[cm]]$f.eid)

}


## Equalizing headers ----------------------------------------------------------

message("Equalizing headers")

# Iterate over all CM dataframes
for (l in 1:(length(dfs) - 1)) {

    col <- names(dfs[[l]])
    # Iterate over the last three dfs
    for (n in 2:length(dfs)) {

        con <- names(dfs[[n]])
        # Check if columns from df2 should be added to df1
        if (length(con[!con %in% col]) > 0) {
            # Save the columns not present in df1
            new <- con[!con %in% col]
            # Add the columns, filling them with NA
            for (c in new) {
                dfs[[l]][c] <- NA
            } # End iteration new columns
        } # End check columns1

        # Check if columns from df1 should be added to df2
        if (length(col[!col %in% con]) > 0) {
            # Save the columns not present in df2
            new <- col[!col %in% con]
            # Add the columns, filling them with NA
            for (c in new) {
                dfs[[n]][c] <- NA
            } # End iteration new columns
        } # End check columns2

        # Correct the order of the columns
        dfs[[n]] <- dfs[[n]] %>% select(order(colnames(.))) %>%
            select(f.eid, everything())
    } # End iteration next dfs

    # Correct the order of the columns
    dfs[[l]] <- dfs[[l]] %>% select(order(colnames(.))) %>%
        select(f.eid, everything())
} # End iteration first dfs

# Combine all dfs 
df <- do.call(bind_rows, dfs)


## Saving data -----------------------------------------------------------------

message(paste0("Saving data in data/raw/", output, ".rds"))
saveRDS(df, paste0("data/raw/", output, ".rds"))
