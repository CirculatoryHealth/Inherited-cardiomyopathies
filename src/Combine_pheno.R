## Script information ---------------------------
##
## Script name: Combine_pheno.R
##
## Purpose of script: Combine the several phenotype files of controls and CMs
##
## Author: M. van Vugt
##
## Date Created: 2021-07-14
##
## Copyright (c) M. van Vugt, 2021
## Email: m.vanvugt-2@umcutrecht.nl
##

## Options ---------------------------

options(scipen = 6, digits = 4) # view outputs in non-scientific notation


# Loading arguments -------------------------------------------------------

# Arguments expected:
#     #1 -- Path to the phenotype files,
#           for example /hpc/dhl_ec/mvanvugt/UKBB
#     #2 -- Suffix of the phenotype files,
#           for example /hpc/dhl_ec/mvanvugt/Software/UKB-pipeline-Utrecht
#     #3 -- Output file name, for example MCM_clean_full

args = commandArgs(trailingOnly = TRUE)
path = args[1] # "data/temp/"
suffix = args[2] # "_raw.txt"
output = args[3] # "MCM_raw_full"


## Loading packages -------------------------------------------------------

library(dplyr)


# Loading data ------------------------------------------------------------

message("Loading data")

dfs <- list()
for (cm in c("ACM", "DCM", "HCM", "Controls")) {

  dfs[[cm]] <- read.delim(paste0(path, cm, suffix), sep = "\t", header = T) %>%
    select(!starts_with("X"))
  dfs[[cm]]$CM <- cm
  dfs[[cm]]$f.eid <- as.character(dfs[[cm]]$f.eid)

}


# Equalizing headers ------------------------------------------------------

message("Equalizing headers")

for (l in 1:(length(dfs) - 1)) {


  col <- names(dfs[[l]])
  for (n in 2:length(dfs)) {

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
    dfs[[n]] <- dfs[[n]] %>% select(order(colnames(.))) %>%
      select(f.eid, everything())
  } # End iteration next dfs
  dfs[[l]] <- dfs[[l]] %>% select(order(colnames(.))) %>%
    select(f.eid, everything())
} # End iteration first dfs
df <- do.call(bind_rows, dfs)


# Saving data -------------------------------------------------------------

message("Saving data")

saveRDS(df, paste0("data/raw/", output, ".rds"))
