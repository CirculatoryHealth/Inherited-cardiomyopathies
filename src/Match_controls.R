#!/usr/bin/env Rscript
##
## Script information ----------------------------------------------------------
##
## Script name: Match_controls.R
##
## Purpose of script: Visualize the UKB population so far included and pick 
##                    controls (matched) based on these characteristics
##
## Author: M. van Vugt
##
## Date Created: 2021/06/01
##
## Copyright (c) M. van Vugt, 2022
## Email: m.vanvugt-2@umcutrecht.nl
##

## Options ---------------------------------------------------------------------

options(scipen = 6, digits = 4) # view outputs in non-scientific notation
set.seed(194)
rm(list=ls())


## Loading packages ------------------------------------------------------------

suppressMessages(library(viridis))


## Loading arguments -----------------------------------------------------------

# Arguments expected:
#     #1 -- Path to the phenotype file of all CMs
#           for example /hpc/dhl_ec/mvanvugt/UKBB/data/raw
#     #2 -- Suffix of the phenotype files of all CMs
#           for example _full.tsv
#     #3 -- Path to the figure directory
#           for example /hpc/dhl_ec/mvanvugt/UKBB/results/figures
#     #4 -- Number of controls to be picked per case
#           for example 4

args   = commandArgs(trailingOnly = TRUE)
path   = args[1] 
suffix = args[2]
figs   = args[3]
mul    = as.numeric(args[4])


## Loading functions -----------------------------------------------------------

source("src/functions.R")


## Loading files ---------------------------------------------------------------

message("Loading data")
# Read the phenotype files per CM
dfs <- list()
dfs$acm <- readr::read_delim(paste0(file.path(path, "ACM"), suffix), "\t", 
                             escape_double = FALSE, trim_ws = TRUE)
dfs$dcm <- readr::read_delim(paste0(file.path(path, "DCM"), suffix), "\t", 
                             escape_double = FALSE, trim_ws = TRUE)
dfs$hcm <- readr::read_delim(paste0(file.path(path, "HCM"), suffix), "\t", 
                             escape_double = FALSE, trim_ws = TRUE)
dfs$acm$f.eid <- as.character(dfs$acm$f.eid)
dfs$dcm$f.eid <- as.character(dfs$dcm$f.eid)
dfs$hcm$f.eid <- as.character(dfs$hcm$f.eid)

# Read the WES phenotype file and CMR ID list
wes <- read.table("data/temp/MCM_WES_pheno.tsv", header = TRUE,
                  stringsAsFactors = FALSE)
wes$f.eid <- as.character(wes$f.eid)
mri <- read.table("data/temp/CMR_IDs.txt", header = TRUE,
                  stringsAsFactors = FALSE)


## Loading files ---------------------------------------------------------------

message("Remove cases from potential controls list")
cmi <- rbind(dfs[[1]][,1], dfs[[2]][,1], dfs[[3]][,1])
wes <- anti_join(wes, cmi, by = "f.eid")

# Add information on MRI
mri$MRI <- "Yes"
wes <- merge(wes, mri, by = "f.eid", all.x = TRUE)
rm(cmi, mri)


## Control data preparation ----------------------------------------------------

message("Characterizing controls")
# Fill up one columns NA's with non-missing values from other time points
wes$f.21000.0.0[is.na(wes$f.21000.0.0)] <- wes$f.21000.1.0[is.na(wes$f.21000.0.0)]
wes$f.21000.0.0[is.na(wes$f.21000.0.0)] <- wes$f.21000.2.0[is.na(wes$f.21000.0.0)]

# Select only necessary columns and rename
wes <- wes %>% dplyr::select(f.eid, f.31.0.0, f.34.0.0, f.21000.0.0,
                             f.21003.0.0, MRI)
names(wes) <- c("f.eid", "Sex", "Year_of_birth", "Ethnic_background",
                "Age_at_recruitment", "MRI")

# Clean up coding
wes$Sex[wes$Sex == 1] <- "Male"
wes$Sex[wes$Sex == 0] <- "Female"
wes$Sex <- as.factor(wes$Sex)

# Make categories of ethnicity
wes$Ethnicity <- NA
wes$Ethnicity[wes$Ethnic_background == "1" | wes$Ethnic_background == "1001" | wes$Ethnic_background == "1002" | wes$Ethnic_background == "1003"] <- "White"
wes$Ethnicity[wes$Ethnic_background == "2" | wes$Ethnic_background == "2001" | wes$Ethnic_background == "2002" | wes$Ethnic_background == "2003" | wes$Ethnic_background == "2004"] <- "Mixed"
wes$Ethnicity[wes$Ethnic_background == "3" | wes$Ethnic_background == "3001" | wes$Ethnic_background == "3002" | wes$Ethnic_background == "3003" | wes$Ethnic_background == "3004"] <- "Asian"
wes$Ethnicity[wes$Ethnic_background == "4" | wes$Ethnic_background == "4001" | wes$Ethnic_background == "4002" | wes$Ethnic_background == "4003"] <- "Black"
wes$Ethnicity[wes$Ethnic_background == "5"] <- "Chinese"
wes$Ethnicity[wes$Ethnic_background == "6"] <- "Other"
wes$Ethnicity <- as.factor(wes$Ethnicity)

# Make categories of age to pick controls in similar age categories
wes$Age_cat <- NA
wes$Age_cat[wes$Age_at_recruitment < 50] <- "<50"
wes$Age_cat[wes$Age_at_recruitment >= 50 & wes$Age_at_recruitment <= 60] <- "50-59"
wes$Age_cat[wes$Age_at_recruitment > 60] <- ">60"
wes$Age_cat <- as.factor(wes$Age_cat)

# Clean up coding for CMR availability
wes$MRI[is.na(wes$MRI)] <- "No"
wes$MRI <- as.factor(wes$MRI)

# Select and make right order of necessary columns
wes <- wes[, c(1, 2, 6, 7, 8)]


## Subset data to pick controls from -------------------------------------------
suc <- list()
suc$m_no <- subset(wes, Sex == "Male" & MRI == "No")
suc$f_no <- subset(wes, Sex == "Female" & MRI == "No")
suc$m_yes <- subset(wes, Sex == "Male" & MRI == "Yes")
suc$f_yes <- subset(wes, Sex == "Female" & MRI == "Yes")
rm(wes)


## Case data preparation -------------------------------------------------------

message("Characterizing cases")
cases <- data.frame()
controls <- data.frame()
# Iterate over the phenotype dfs per CM
for (q in 1:length(dfs)) {

    # Temporarily save CM-specific df in df
    df <- dfs[[q]]
    # Fill up one columns NA's with non-missing values from other time points
    df$f.21000.0.0[is.na(df$f.21000.0.0)] <- df$f.21000.1.0[is.na(df$f.21000.0.0)]
    df$f.21000.0.0[is.na(df$f.21000.0.0)] <- df$f.21000.2.0[is.na(df$f.21000.0.0)]

    # Select only necessary columns and rename
    df <- df %>% dplyr::select(f.eid, f.31.0.0, f.34.0.0, f.21000.0.0,
                               f.21003.0.0, Age_MIR)
    names(df) <- c("f.eid", "Sex", "Year_of_birth", "Ethnic_background",
                   "Age_at_recruitment", "Age_MRI")

    # Clean up coding for CMR availability
    df$MRI <- NA
    df$MRI[is.na(df$Age_MRI)] <- "No"
    df$MRI[!is.na(df$Age_MRI)] <- "Yes"
    df$MRI <- as.factor(df$MRI)

    # Clean up coding
    df$Sex[df$Sex == 1] <- "Male"
    df$Sex[df$Sex == 0] <- "Female"
    df$Sex <- as.factor(df$Sex)

    # Make categories of ethnicity
    df$Ethnicity <- NA
    df$Ethnicity[df$Ethnic_background == "1" | df$Ethnic_background == "1001" | df$Ethnic_background == "1002" | df$Ethnic_background == "1003"] <- "White"
    df$Ethnicity[df$Ethnic_background == "2" | df$Ethnic_background == "2001" | df$Ethnic_background == "2002" | df$Ethnic_background == "2003" | df$Ethnic_background == "2004"] <- "Mixed"
    df$Ethnicity[df$Ethnic_background == "3" | df$Ethnic_background == "3001" | df$Ethnic_background == "3002" | df$Ethnic_background == "3003" | df$Ethnic_background == "3004"] <- "Asian"
    df$Ethnicity[df$Ethnic_background == "4" | df$Ethnic_background == "4001" | df$Ethnic_background == "4002" | df$Ethnic_background == "4003"] <- "Black"
    df$Ethnicity[df$Ethnic_background == "5"] <- "Chinese"
    df$Ethnicity[df$Ethnic_background == "6"] <- "Other"
    df$Ethnicity <- as.factor(df$Ethnicity)

    # Make categories of age to pick controls in similar age categories
    df$Age_cat <- NA
    df$Age_cat[df$Age_at_recruitment < 50] <- "<50"
    df$Age_cat[df$Age_at_recruitment >= 50 & df$Age_at_recruitment <= 60] <- "50-59"
    df$Age_cat[df$Age_at_recruitment > 60] <- ">60"
    df$Age_cat <- as.factor(df$Age_cat)

    # Select and make right order of necessary columns
    df <- df[, c(1, 2, 5, 7, 8, 9)]

    # Subset data 
    mri <- subset(df, MRI == "Yes")
    rest <- subset(df, MRI == "No")
    cm <- names(dfs)[q]


    # Plotting ----------------------------------------------------------------

    message("Making some nice graphs")
    p <- rbind(perc_var(mri, "Sex", "MRI"), perc_var(rest, "Sex", "non-MRI")) %>%
        ggplot(aes(fill = value, x = name, y = perc)) +
        geom_col(position = position_dodge(width = 0.7),
                 color = "black", width = 0.7) +
        labs(x = "Group", y = "Percentage of cases") +
        ggtitle(paste0("Sex distribution ", toupper(cm))) +
        geom_text(aes(label = count), vjust = -0.7,
                  position = position_dodge(width = 0.7)) +
        scale_fill_manual(name = "Sex", labels = c("Female", "Male"),
                  values = c("palevioletred1", "skyblue1"))
    ggsave(paste0(file.path(figs, "Bar_sex_all_"), toupper(cm), ".png"), 
           plot = p)

    p <- rbind(perc_var(mri, "Ethnicity", "MRI"), 
               perc_var(rest, "Ethnicity", "non-MRI")) %>%
        ggplot(aes(fill = value, x = name, y = perc)) +
        geom_col(position = position_dodge(width = 0.7),
                 color = "black", width = 0.7) +
        labs(x = "", y = "Percentage of cases") +
        ggtitle(paste0("Ethnicities ", toupper(cm))) +
        geom_text(aes(label = count), vjust = -0.7,
                  position = position_dodge(width = 0.7)) +
        scale_fill_viridis(discrete = TRUE, name = "Ethnic group")
    ggsave(paste0(file.path(figs, "Bar_ethnic_all_"), toupper(cm), ".png"), 
           plot = p)

    p <- ggplot(mri, aes(Age_at_recruitment)) +
        geom_histogram(color = "black", fill = "red3",
                       aes(y = (..count..)/sum(..count..)), binwidth = 10) +
        stat_bin(aes(y=(..count..)/sum(..count..),
                     label=(..count..)),
                 geom="text", size=4, binwidth = 10, vjust=-.7) +
        facet_wrap(~Sex) +
        labs(x = "Age (in years)", y = "Percentage") +
        ggtitle(paste0("Age at recruitment - MRI ", toupper(cm)))
    ggsave(paste0(file.path(figs, "Histogram_age_sex_MRI_"), toupper(cm), ".png"), 
           plot = p)

    p <- ggplot(rest, aes(Age_at_recruitment)) +
        geom_histogram(color = "black", fill = "red3",
                       aes(y = (..count..)/sum(..count..)), binwidth = 10) +
        stat_bin(aes(y=(..count..)/sum(..count..),
                     label=(..count..)),
                 geom="text", size=4, binwidth = 10, vjust=-.7) +
        facet_wrap(~Sex) +
        labs(x = "Age (in years)", y = "Percentage") +
        ggtitle(paste0("Age at recruitment - nonMRI ", toupper(cm)))
    ggsave(paste0(file.path(figs, "Histogram_age_sex_nonMRI_"), toupper(cm), ".png"),
           plot = p)

    p <- ggplot(df, aes(Age_at_recruitment)) +
        geom_histogram(color = "black", fill = "red3",
                       aes(y = (..count..)/sum(..count..)), binwidth = 10) +
        stat_bin(aes(y=(..count..)/sum(..count..),
                     label=(..count..)),
                 geom="text", size=4, binwidth = 10, vjust=-.7) +
        facet_wrap(~Sex) +
        labs(x = "Age (in years)", y = "Percentage") +
        ggtitle(paste0("Age at recruitment - all ", toupper(cm)))
    ggsave(paste0(file.path(figs, "Histogram_age_sex_all_"), toupper(cm), ".png"),
           plot = p)

    # Make column with CM and bind to case df
    df$CM <- toupper(names(dfs)[q])
    cases <- rbind(cases, df)
} # end for-loop iteration CMs

cases$CM <- as.factor(cases$CM)


## Subset case data to determine categories ------------------------------------
message("Picking controls")
sub <- list()
sub$m_no <- subset(cases, Sex == "Male" & MRI == "No")
sub$f_no <- subset(cases, Sex == "Female" & MRI == "No")
sub$m_yes <- subset(cases, Sex == "Male" & MRI == "Yes")
sub$f_yes <- subset(cases, Sex == "Female" & MRI == "Yes")


## Pick controls based on case categories --------------------------------------

for (s in 1:length(sub)) {
    ok <- sub[[s]]
    for (race in levels(ok$Ethnicity)) {
        for (age in levels(ok$Age_cat)) {
            # Determine number of cases to be extracted
            n <- ok %>% filter(Ethnicity == race) %>%
                filter(Age_cat == age) %>% nrow() * mul
            con <- suc[[s]] %>% filter(Ethnicity == race) %>%
                filter(Age_cat == age)
            rows <- sample(1:nrow(con), n)
            for (t in rows) {
                id <- con[t, 1]
                ind <- subset(con, f.eid == id)
                controls <- rbind(controls, ind)
            } # end for-loop picking controls
        } # end for-loop ages
    } # end for-loop races
    nas <- subset(ok, is.na(Ethnicity))
    for (age in levels(nas$Age_cat)) {
        m <- nas %>% filter(Age_cat == age) %>% nrow() * mul
        nac <- suc[[s]] %>% filter(is.na(Ethnicity)) %>%
            filter(Age_cat == age)
        row <- sample(1:nrow(nac), m)
        for (r in row) {
            id <- nac[r, 1]
            ind <- subset(nac, f.eid == id)
            controls <- rbind(controls, ind)
        } # end for-loop picking controls
    } # end for-loop nas-ages
} # end for-loop subset
rm(ok, race, age, n, con, rows, ind, id, t, m, row, nac, nas, s, r)


## Save controls ID list -------------------------------------------------------
message("Saving controls ID-list")
cids <- as.data.frame(controls$f.eid)
write.table(cids, file.path(path, "Control_IDs.tsv"), col.names = "f.eid",
            row.names = FALSE, quote = FALSE)



