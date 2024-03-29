#!/usr/bin/env Rscript
##
## Script information ----------------------------------------------------------
##
## Script name: MCM_pheno_clean.R
##
## Purpose of script: Clean up the phenotype file
##
## Author: Arjen Cupido
## Modified by: M. van Vugt
##
## Date Created: 2021/07/13
##
## Copyright (c) M. van Vugt, 2022
## Email: m.vanvugt-2@umcutrecht.nl
##

## Options ---------------------------------------------------------------------

options(scipen = 6, digits = 4) # view outputs in non-scientific notation


## Loading packages ------------------------------------------------------------

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(stringr))


## Loading arguments -----------------------------------------------------------

# Arguments expected:
#     #1 -- Path to and name of the phenotype file,
#           for example /hpc/dhl_ec/mvanvugt/UKBB/data/raw/MCM_raw_full.rds
#     #2 -- Directory to the helper files,
#           for example /hpc/dhl_ec/mvanvugt/UKBB/data/raw
#     #3 -- Output directory, 
#           for example /hpc/dhl_ec/mvanvugt/UKBB/data/processed
#     #4 -- Prefix of the output files, for example MCM

args   = commandArgs(trailingOnly = TRUE)
pheno  = args[1] 
files  = args[2]
output = args[3] 
prefix = args[4]


## Loading files ---------------------------------------------------------------

message(paste0("Loading data from ", pheno))
d <- data.table(readRDS(pheno))
d$CM <- as.factor(d$CM)
cols <- colnames(d)


## Preparing helper files ------------------------------------------------------

message("Preparing helper files...")
desc <- read.table(file.path(files, "Dictionary_Field_ID_desc.tsv"),
                   header = TRUE, sep = "\t", quote = "",
                   stringsAsFactors = FALSE)
# Modify this dictionary, so replacements are smoother
for (line in 1:nrow(desc)) {
    x <- desc[line, 1]
    desc[line, 1] <- paste0("f.", x, ".")
    desc[line, 2] <- paste0(desc[line, 2], ".")
}
rm(x, line)

summ <- read.table(file.path(files, "Summary_outcomes.txt"), sep = "\t",
                   header = TRUE, stringsAsFactors = FALSE, quote = "")
rep <- read.table(file.path(files, "UKB_repeat_measurements.txt"),
                  header = FALSE, stringsAsFactors = FALSE, quote = "")
nas <- read.table(file.path(files, "UKB_missing_merge.txt"),
                  header = FALSE, stringsAsFactors = FALSE, quote = "")


## Calculate mean height, wheight and BSA --------------------------------------

message("Preparing some measurements")
temp <- d %>% select(f.eid, starts_with("f.50")) %>% as.data.frame
temp$Standing_height <- rowMeans(temp[, 2:ncol(temp)], na.rm = TRUE)
temp <- temp[, c(1, ncol(temp))]
d <- merge(d, temp, all.x = TRUE)
temp <- d %>% select(f.eid, starts_with("f.21002")) %>% as.data.frame
temp$UKB_weight <- rowMeans(temp[, 2:ncol(temp)], na.rm = TRUE)
temp <- temp[, c(1, ncol(temp))]
d <- merge(d, temp, all.x = TRUE)
d$BMI <- d$UKB_weight / ((d$Standing_height / 100) * (d$Standing_height / 100))

# BSA calculations
d$BSA <- 0.20247 * (d$UKB_weight ^ 0.425) * ((d$Standing_height / 100 ) ^ 0.725)


## Defining fields -------------------------------------------------------------

# Specify vectors that contain the columns in which certain data is stored
# This makes it easier to search in these columns later on
message("Defining fields...")
main_date <- grep("41280", cols)
sec_date <- grep("41262", cols)

death_date <- grep("40000", cols)

ICD9 <- grep("41271|41203|41205", cols)
ICD9date <-grep("41281|41263", cols, value = T)

proc <- grep("41272", cols)
procdate <- grep("41282", cols)

death_main <- grep("f.40001", cols)
death_oth <- grep("f.40002|f.41201", cols)
medication <- grep("f.20003.0", cols)
med_spec <- grep("6177", cols)

PC <- grep("f.22009", cols) # 109-148

self_rep <- grep("f.20002", cols) # 4 - 105

death <- c(death_main, death_oth)


## Construc personally defined outcomes ----------------------------------------

message("Construct personally defined outcomes...")
pbar <- txtProgressBar(min = 0, max = nrow(summ), style = 3)

for (line in 1:nrow(summ)) {

    name <- summ[line, 1]

    # Something strange happened for this outcome, so do that separately
    if (name == "Family_Heart_Disease_sum") {

        ids <- d %>% select(f.eid)
        # Save and extract columns where family heart disease is saved in df
        other <- grep(summ[line, 9], cols)
        o <- d %>% select(f.eid, any_of(names(d)[other])) %>% as.data.frame
        odf <- data.frame()
        # Iterate over columns in which family heart disease info is stored
        for (col in 2:ncol(o)) {
            # Add people with family heart disease to df
            new <- subset(o, o[, col] == 1)
            new <- new %>% select(f.eid)
            odf <- rbind(odf, new)
        } # End iteration columns
        # Create column indicating family heart disease
        tmp.name <- gsub("sum", "other", name)
        odf[, tmp.name] <- 1
        # Merge with the full df
        odf <- merge(unique(ids), unique(odf), by = "f.eid", all.x = TRUE)
        d <- merge(d, odf, by = "f.eid", all.x = TRUE)

    } else {

        if (!is.na(summ[line, 2])) {
            # Save columns in which medication information is stored
            meds <- grep(summ[line, 2], cols)
            tmp.name <- gsub("sum", "meds", name)
            # Determine whether people have medication for the specified outcome
            d[, tmp.name] <- as.integer(apply(d[, ..meds], 1, function(r) any(r %in% c(grep(summ[line, 3], r, value = T)))))
        } # End if medication

        if (!is.na(summ[line, 4])) {
            # Save columns in which ICD10 codes are stored
            icd10 <- grep(summ[line, 4], cols)
            tmp.name <- gsub("sum", "icd10", name)
            # Determine if people have ICD10 code for specified outcome
            d[, tmp.name] <- as.integer(apply(d[, ..icd10], 1, function(r) any(r %in% c(grep(summ[line, 5], r, value = T)))))
        } # End if ICD10

        if (!is.na(summ[line, 6])) {
            # Save columns in which cause of death is stored
            death <- grep(summ[line, 6], cols)
            tmp.name <- gsub("sum", "death", name)
            # Determine if people died due to specified outcome
            d[, tmp.name] <- as.integer(apply(d[, ..death], 1, function(r) any(r %in% c(grep(summ[line, 5], r, value = T)))))
        } # End if death

        if (!is.na(summ[line, 7])) {
            # Save columns in which self-reported diseases are stored
            sr <- grep(summ[line, 7], cols)
            tmp.name <- gsub("sum", "sr", name)
            # Determine if people reported they suffer from  specified outcome
            d[, tmp.name] <- as.integer(apply(d[, ..sr], 1, function(r) any(r %in% c(grep(summ[line, 8], r, value = T)))))
        } # End if self-reported

        if (!is.na(summ[line, 9])) {
            # Save column in which other information required for outcome is stored
            other <- grep(summ[line, 9], cols)
            tmp.name <- gsub("sum", "other", name)
            # Determine if people fulfill the criteria for this outcome
            d[, tmp.name] <- as.integer(apply(d[, ..other], 1, function(r) any(r %in% c(grep(summ[line, 10], r, value = T)))))
        } # End if other

    } # Check for Family_Heart_Disease

    # Summarize all information per outcome in one column
    tmp <- d %>% select(f.eid, starts_with(gsub("sum", "", name)))
    tmp[, name] <- ifelse(rowSums(tmp[, 2:ncol(tmp)], na.rm = TRUE) >= 1, "Yes", "No")
    tmp <- tmp %>% select(f.eid, all_of(name))
    # Add summarizing column to df
    d <- merge(d, unique(tmp), by = "f.eid", all.x = TRUE)
    names(d)[ncol(d)] <- name

    setTxtProgressBar(pbar, line)

} # End for-loop defining outcomes
close(pbar)


## Add date of diagnosis -------------------------------------------------------

message("Adding diagnosis dates...")
pbar <- txtProgressBar(min = 0, max = nrow(summ), style = 3)

for (line in 1:nrow(summ)) {

    if (!is.na(summ[line, 4])) {

        name <- summ[line, 1]
        # Select people that have diagnosis for specified outcome
        f <- d[get(name) == "Yes"]

        if (nrow(f) > 0) {

            datelist <- list()
            # Iterate over diagnosed individuals
            for (i in 1:nrow(f)) {

                # Selecting row
                x <- f[i, ]

                datecolumn <- vector()
                # Select columns for which date is available
                icd <- c(base::strsplit(summ[line, 4], "[|]")[[1]], base::strsplit(summ[line, 6], "[|]")[[1]])
                icd <- icd[!icd %in% c("41201", "41204", "f.40002", "f.41201")]

                # Iterate over the columns for which date is available
                for (code in icd) {
                    if (code == "41270") {
                        # Select columns that contain diagnosis
                        sel <- grep(code, cols)
                        # Determine location of data columns
                        daf <- c(apply(x[,..sel], 1, function(r) which(r %in% grep(pattern = summ[line, 5], r, value = T))))
                        datecolumn <- c(datecolumn, main_date[daf])

                    } else if (code == "41202") {
                        # Select columns that contain diagnosis
                        sel <- grep(code, cols)
                        # Determine location of data columns
                        daf <- c(apply(x[,..sel], 1, function(r) which(r %in% grep(pattern = summ[line, 5], r, value = T))))
                        datecolumn <- c(datecolumn, sec_date[daf])

                    } else if (code == "f.40001") {
                        # Select columns that contain diagnosis
                        sel <- grep(code, cols)
                        # Determine location of data columns
                        daf <- c(apply(x[,..sel], 1, function(r) which(r %in% grep(pattern = summ[line, 5], r, value = T))))
                        datecolumn <- c(datecolumn, death_date[daf])

                    } # End check for available date

                } # End iteration diagnoses

                # Select first date
                datesframe <- as.data.frame(subset(x, select = datecolumn))
                datesframe[names(datesframe)] <- lapply(datesframe[names(datesframe)], as.Date)

                if (nrow(datesframe) > 0) {

                    r <- c(apply(datesframe, 1, function(r) min(r)))
                    colnames(datesframe) <- as.character(1:ncol(datesframe))

                    # Save first date
                    date.name <- gsub("sum", "FirstDate", name)
                    datelist[[i]] <- data.frame("f.eid" = x$f.eid, date.name = r, datesframe)
                    names(datelist[[i]])[2] <- date.name

                } # End check empty datesframe

            } # End iteration rows subset

            # Save all other dates
            dateresults <- do.call(bind_rows, datelist)
            colnames(dateresults) <- gsub("X", gsub("sum", "Date", name),
                                          colnames(dateresults))
            # Add date columns to full df
            d <- left_join(d, unique(dateresults), by = "f.eid")

        } # End if subsetting exists

        setTxtProgressBar(pbar, line)

    } # Check if defined by diagnoses

} # End iteration rows outcomes
close(pbar)


## Compare dates ---------------------------------------------------------------

# To define people with non-ischemic HF, we compare the first date of diagnosis
# of both outcomes. People that already had HF diagnosis before ischemia, are
# included as non-ischemic HF anyways

# Select people with both HF and ischemia dianosis
f <- d %>% filter(Heart_Failure_sum == "Yes" & Chronic_ischaemic_heart_disease_sum == "Yes") %>%
    select(f.eid, ends_with("_sum"), ends_with("FirstDate"))
f$Chronic_ischaemic_heart_disease_FirstDate <- as.Date(f$Chronic_ischaemic_heart_disease_FirstDate)
f$Heart_Failure_FirstDate <- as.Date(f$Heart_Failure_FirstDate)
# Calculate differences between the diagnosis dates
f$HF_isch_diff <- as.numeric(f$Chronic_ischaemic_heart_disease_FirstDate - f$Heart_Failure_FirstDate)
# Save dates
write.table(f, file.path(output, "Diagnosis_dates.tsv"), sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)


## Renaming columns ------------------------------------------------------------

message("Renaming columns...")
pbar <- txtProgressBar(min = 0, max = nrow(desc), style = 3)

# The code of the column names are replaced by the name given in the Data
# Dictionary for that column, but the instance and array-numbers are kept
# For example f.31.0.0 becomes Sex.0.0
for (l in 1:nrow(desc)) {
    # Save fieldID and description of this fieldID
    id <- desc[l, 1] 
    nam <- desc[l, 2] 
    # Replace fieldID number for its description
    colnames(d) <- gsub(id, nam, colnames(d), fixed = TRUE)

    setTxtProgressBar(pbar, l)
} # end for-loop total
rm(id, nam, l)
close(pbar)


## Averaging repeated measurements ---------------------------------------------

message("Averaging repeated measurements...")
pbar <- txtProgressBar(min = 0, max = nrow(rep), style = 3)

# There are some repeat measurements, and using this loop, the average will
# be taken. 
# For example, all measurements of Systolic_blood_pressure_manual_reading at
# instance 0 (column name: Systolic_blood_pressure_manual_reading.0.* where *
# means any array number) will be averaged and the average column will have
# the name: Systolic_blood_pressure_manual_reading.0_mean.
for (line in 1:nrow(rep)) {
    # To calculate per instance, we iterate over the instances
    for (inst in (0:20)) {
        # Select columns with this measurement
        x <- paste0(rep[line, ], ".", inst)
        temp <- dplyr::select(d, c(f.eid, dplyr::starts_with(x)))
        if (ncol(temp) > 1) {
            # Calculate mean over all these columns
            temp[, paste0(x, "_mean")] <- rowMeans(temp[, 2:ncol(temp)], na.rm = TRUE)
            # Add mean to full df
            temp <- temp %>% select(f.eid, ends_with("mean"))
            d <- merge(d, unique(temp), by = "f.eid", all = TRUE)
        } else {
            break
        } # end ifelse check for instance loop
    } # end for-instance loop

    setTxtProgressBar(pbar, line)

} # end-averaging loop
rm(line, inst, x, temp)
close(pbar)


## Manipulating some columns ---------------------------------------------------

message("Manipulating some columns")
d <- as.data.frame(d)
cols <- names(d)
for (line in 1:nrow(nas)) {
    sels <- strsplit(nas[line, 2], ";")[[1]]
    for (s in sels) {
        # Select columns for which NAs can be filled with values of other timepoints
        if (s == sels[1]) {
            temp <- dplyr::select(d, c(f.eid, starts_with(paste0(s, "."))))
        } else {
            temp1 <- dplyr::select(d, c(f.eid, starts_with(s)))
            temp <- merge(temp, unique(temp1), all = TRUE, by = "f.eid")
        } # end ifelse-loop for initiating new temporary df
    } # end for-loop iterating over dplyr::selecting columns

    if (ncol(temp) > 1) {
        temp[nas[line, 1]] <- temp[, 2]
        # Fill up NAs with values from other time points
        for (c in 3:(ncol(temp)-1)) {
            temp[, ncol(temp)][is.na(temp[, ncol(temp)])] <- temp[, c][is.na(temp[, ncol(temp)])]
        } # end for-loop iterating over columns
        # Select only f.eid and newly created column and add to full df
        temp <- temp[, c(1, ncol(temp))]
        d <- merge(d, unique(temp), by = "f.eid", all = TRUE)
        names(d)[ncol(d)] <- nas[line, 1]
    } else {
        message(paste0("Columns for this column don't exist:", nas[line, n]))
    } # end ifelse loop existing columns
} # end for-loop
rm(line, sels, s, temp, temp1, c, cols)

# In this loop we will check whether Ethnic_background.0.0 is present in the
# df and if so, the information will be used to summarize. First, NA's will
# be replaced by any value in other instances of the same column and then
# They will be summarized in less specific categories to make it more
# readable.
if ("Ethnic_background.0.0" %in% names(d)) {
    # Fill up NAs with values from other instances
    d$Ethnic_background.0.0[is.na(d$Ethnic_background.0.0)] <- d$Ethnic_background.1.0[is.na(d$Ethnic_background.0.0)]
    d$Ethnic_background.0.0[is.na(d$Ethnic_background.0.0)] <- d$Ethnic_background.2.0[is.na(d$Ethnic_background.0.0)]

    # Categorise Ethnicity
    d$Ethnicity[d$Ethnic_background.0.0 == "1" | d$Ethnic_background.0.0 == "1001" | d$Ethnic_background.0.0 == "1002" | d$Ethnic_background.0.0 == "1003"] <- "White"
    d$Ethnicity[d$Ethnic_background.0.0 == "2" | d$Ethnic_background.0.0 == "2001" | d$Ethnic_background.0.0 == "2002" | d$Ethnic_background.0.0 == "2003" | d$Ethnic_background.0.0 == "2004"] <- "Mixed"
    d$Ethnicity[d$Ethnic_background.0.0 == "3" | d$Ethnic_background.0.0 == "3001" | d$Ethnic_background.0.0 == "3002" | d$Ethnic_background.0.0 == "3003" | d$Ethnic_background.0.0 == "3004"] <- "Asian"
    d$Ethnicity[d$Ethnic_background.0.0 == "4" | d$Ethnic_background.0.0 == "4001" | d$Ethnic_background.0.0 == "4002" | d$Ethnic_background.0.0 == "4003"] <- "Black"
    d$Ethnicity[d$Ethnic_background.0.0 == "5"] <- "Chinese"
    d$Ethnicity[d$Ethnic_background.0.0 == "6"] <- "Other"
    d$Ethnicity[d$Ethnic_background.0.0 <= 0] <- NA
    d$Ethnicity[d$Ethnic_background.0.0 == "0"] <- NA
    d$Ethnicity <- as.factor(d$Ethnicity)
} # end if-loop make ethnicity column if exists

# The Sex-column, if present, is decoded into Male/Female instead of 1/0
if("Sex.0.0" %in% names(d)) {
    d$Sex[d$Sex.0.0 == 1] <- "Male"
    d$Sex[d$Sex.0.0 == 0] <- "Female"
    d$Sex <- as.factor(d$Sex)
} # end if-loop categories sex if exists

# BMI is calcluated based on, if existing, the Height and Weight column.
if ("Height" %in% names(d) && "Weight" %in% names(d)) {
    d$BMI <- d$Weight / ((d$Height/100) * (d$Height/100))
} # end if-loop calculating BMI if columns exist

# Take mean of automated and manual BP reading
temp <- d %>%
    select(f.eid, starts_with("Systolic_blood_pressure")) %>%
    select(f.eid, ends_with("0_mean"))
temp$Systolic_blood_pressure_mean <- rowMeans(temp[, 2:ncol(temp)], na.rm = TRUE)
d <- merge(d, unique(temp))
temp <- d %>%
    select(f.eid, starts_with("Diastolic_blood_pressure")) %>%
    select(f.eid, ends_with("0_mean"))
temp$Diastolic_blood_pressure_mean <- rowMeans(temp[, 2:ncol(temp)], na.rm = TRUE)
d <- merge(d, unique(temp))


## Options ---------------------------------------------------------------------

message(paste0("Saving data in ", file.path(output, prefix), "_cleaned.tsv"))

saveRDS(d, paste0(file.path(output, prefix), "_cleaned.rds"))
write.table(d, file = paste0(file.path(output, prefix), "_cleaned.tsv"), 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


