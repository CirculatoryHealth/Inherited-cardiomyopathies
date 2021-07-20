## Script information ---------------------------
##
## Script name: UKB_MCM_controls.R
##
## Purpose of script: Summarize (graphically) the UKB population for UKB MCM-
## project to pick controls.
##
## Author: M. van Vugt
##
## Date Created: 2021-06-01
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
path = args[1] # "data/temp"
suffix = args[2] # "_full.tsv"
figs = args[3] # "Results/figures/UKB_MCM_Summary.pptx"


## Loading packages ---------------------------

library(viridis)


## Loading functions ---------------------------

rm(list=setdiff(ls(), "all"))
rm(list=ls())
source("src/functions.R")


# Loading data ------------------------------------------------------------

dfs <- list()
dfs$hcm <- readr::read_delim(paste0(path, "/HCM", suffix), "\t", escape_double = FALSE, trim_ws = TRUE)
dfs$dcm <- readr::read_delim(paste0(path, "/DCM", suffix), "\t", escape_double = FALSE, trim_ws = TRUE)
dfs$acm <- readr::read_delim(paste0(path, "/ACM", suffix), "\t", escape_double = FALSE, trim_ws = TRUE)

PptPath <- figs
mul <- 4 # Number of controls to be picked per case

wes <- read.table(paste0(path, "/MCM_WES_pheno.tsv"), header = TRUE,
                  stringsAsFactors = FALSE)
mri <- read.table(paste0(path, "/CMR_IDs.txt"), header = TRUE,
                  stringsAsFactors = FALSE)

# Filter out cases
cmi <- rbind(dfs[[1]][,1], dfs[[2]][,1], dfs[[3]][,1])
wes <- anti_join(wes, cmi, by = "f.eid")

# Add information on MRI
mri$MRI <- "Yes"
wes <- merge(wes, mri, by = "f.eid", all.x = TRUE)
rm(cmi, mri)


# Control data preparation ------------------------------------------------

wes$f.21000.0.0[is.na(wes$f.21000.0.0)] <- wes$f.21000.1.0[is.na(wes$f.21000.0.0)]
wes$f.21000.0.0[is.na(wes$f.21000.0.0)] <- wes$f.21000.2.0[is.na(wes$f.21000.0.0)]

wes <- wes %>% dplyr::select(f.eid, f.31.0.0, f.34.0.0, f.21000.0.0,
                             f.21003.0.0, MRI)
names(wes) <- c("f.eid", "Sex", "Year_of_birth", "Ethnic_background",
                "Age_at_recruitment", "MRI")

wes$Sex[wes$Sex == 1] <- "Male"
wes$Sex[wes$Sex == 0] <- "Female"
wes$Sex <- as.factor(wes$Sex)

wes$Ethnicity <- NA
wes$Ethnicity[wes$Ethnic_background == "1" | wes$Ethnic_background == "1001" | wes$Ethnic_background == "1002" | wes$Ethnic_background == "1003"] <- "White"
wes$Ethnicity[wes$Ethnic_background == "2" | wes$Ethnic_background == "2001" | wes$Ethnic_background == "2002" | wes$Ethnic_background == "2003" | wes$Ethnic_background == "2004"] <- "Mixed"
wes$Ethnicity[wes$Ethnic_background == "3" | wes$Ethnic_background == "3001" | wes$Ethnic_background == "3002" | wes$Ethnic_background == "3003" | wes$Ethnic_background == "3004"] <- "Asian"
wes$Ethnicity[wes$Ethnic_background == "4" | wes$Ethnic_background == "4001" | wes$Ethnic_background == "4002" | wes$Ethnic_background == "4003"] <- "Black"
wes$Ethnicity[wes$Ethnic_background == "5"] <- "Chinese"
wes$Ethnicity[wes$Ethnic_background == "6"] <- "Other"
wes$Ethnicity <- as.factor(wes$Ethnicity)

wes$Age_cat <- NA
wes$Age_cat[wes$Age_at_recruitment < 50] <- "<50"
wes$Age_cat[wes$Age_at_recruitment >= 50 & wes$Age_at_recruitment <= 60] <- "50-59"
wes$Age_cat[wes$Age_at_recruitment > 60] <- ">60"
wes$Age_cat <- as.factor(wes$Age_cat)

wes$MRI[is.na(wes$MRI)] <- "No"
wes$MRI <- as.factor(wes$MRI)

wes <- wes[, c(1, 2, 6, 7, 8)]

# Subset based on MRI and Sex
suc <- list()
suc$m_no <- subset(wes, Sex == "Male" & MRI == "No")
suc$f_no <- subset(wes, Sex == "Female" & MRI == "No")
suc$m_yes <- subset(wes, Sex == "Male" & MRI == "Yes")
suc$f_yes <- subset(wes, Sex == "Female" & MRI == "Yes")
rm(wes)

cases <- data.frame()
controls <- data.frame()
for (q in 1:length(dfs)) {

  # Case Data Preparation ---------------------------------------------------

  df <- dfs[[q]]
  df$f.21000.0.0[is.na(df$f.21000.0.0)] <- df$f.21000.1.0[is.na(df$f.21000.0.0)]
  df$f.21000.0.0[is.na(df$f.21000.0.0)] <- df$f.21000.2.0[is.na(df$f.21000.0.0)]

  df <- df %>% dplyr::select(f.eid, f.31.0.0, f.34.0.0, f.21000.0.0,
                               f.21003.0.0, Age_MRI)
  names(df) <- c("f.eid", "Sex", "Year_of_birth", "Ethnic_background",
                  "Age_at_recruitment", "Age_MRI")
  df$MRI <- NA
  df$MRI[is.na(df$Age_MRI)] <- "No"
  df$MRI[!is.na(df$Age_MRI)] <- "Yes"
  df$MRI <- as.factor(df$MRI)

  df$Sex[df$Sex == 1] <- "Male"
  df$Sex[df$Sex == 0] <- "Female"
  df$Sex <- as.factor(df$Sex)

  df$Ethnicity <- NA
  df$Ethnicity[df$Ethnic_background == "1" | df$Ethnic_background == "1001" | df$Ethnic_background == "1002" | df$Ethnic_background == "1003"] <- "White"
  df$Ethnicity[df$Ethnic_background == "2" | df$Ethnic_background == "2001" | df$Ethnic_background == "2002" | df$Ethnic_background == "2003" | df$Ethnic_background == "2004"] <- "Mixed"
  df$Ethnicity[df$Ethnic_background == "3" | df$Ethnic_background == "3001" | df$Ethnic_background == "3002" | df$Ethnic_background == "3003" | df$Ethnic_background == "3004"] <- "Asian"
  df$Ethnicity[df$Ethnic_background == "4" | df$Ethnic_background == "4001" | df$Ethnic_background == "4002" | df$Ethnic_background == "4003"] <- "Black"
  df$Ethnicity[df$Ethnic_background == "5"] <- "Chinese"
  df$Ethnicity[df$Ethnic_background == "6"] <- "Other"
  df$Ethnicity <- as.factor(df$Ethnicity)

  df$Age_cat <- NA
  df$Age_cat[df$Age_at_recruitment < 50] <- "<50"
  df$Age_cat[df$Age_at_recruitment >= 50 & df$Age_at_recruitment <= 60] <- "50-59"
  df$Age_cat[df$Age_at_recruitment > 60] <- ">60"
  df$Age_cat <- as.factor(df$Age_cat)

  df <- df[, c(1, 2, 5, 7, 8, 9)]

  mri <- subset(df, MRI == "Yes")
  rest <- subset(df, MRI == "No")
  cm <- names(dfs)[q]


  # Plotting ----------------------------------------------------------------

  p <- rbind(perc_var(mri, "Sex", "MRI"), perc_var(rest, "Sex", "non-MRI")) %>%
    ggplot(aes(fill = value, x = name, y = perc)) +
    geom_col(position = position_dodge(width = 0.7),
             color = "black", width = 0.7) +
    labs(x = "Group", y = "Percentage of cases") +
    ggtitle(paste0("Sex distribution ", toupper(cm))) +
    geom_text(aes(label = count), vjust = -0.7,
              position = position_dodge(width = 0.7)) +
    scale_fill_manual(name = "Sex", labels = c("Female", "Male"),
                      values = c("palevioletred1", "skyblue1")) +
    my_theme()
  ggsave(paste0("Results/figures/Bar_sex_all_",
                toupper(cm), ".png"),
         plot = p)
  create_pptx(p, PptPath)


  p <- rbind(perc_var(mri, "Ethnicity", "MRI"), perc_var(rest, "Ethnicity", "non-MRI")) %>%
    ggplot(aes(fill = value, x = name, y = perc)) +
    geom_col(position = position_dodge(width = 0.7),
             color = "black", width = 0.7) +
    labs(x = "", y = "Percentage of cases") +
    ggtitle(paste0("Ethnicities ", toupper(cm))) +
    geom_text(aes(label = count), vjust = -0.7,
              position = position_dodge(width = 0.7)) +
    scale_fill_viridis(discrete = TRUE, name = "Ethnic group") +
    my_theme()
  ggsave(paste0("Results/figures/Bar_ethnic_all_",
                toupper(cm), ".png"),
         plot = p)
  create_pptx(p, PptPath)

  p <- ggplot(mri, aes(Age_at_recruitment)) +
    geom_histogram(color = "black", fill = "red3",
                   aes(y = (..count..)/sum(..count..)), binwidth = 10) +
    stat_bin(aes(y=(..count..)/sum(..count..),
                 label=(..count..)),
             geom="text", size=4, binwidth = 10, vjust=-.7) +
    facet_wrap(~Sex) +
    labs(x = "Age (in years)", y = "Percentage") +
    ggtitle(paste0("Age at recruitment - MRI ", toupper(cm))) +
    my_theme()
  ggsave(paste0("Results/figures/Histogram_age_sex_MRI_",
                toupper(cm), ".png"),
         plot = p)
  create_pptx(p, PptPath)

  p <- ggplot(rest, aes(Age_at_recruitment)) +
    geom_histogram(color = "black", fill = "red3",
                   aes(y = (..count..)/sum(..count..)), binwidth = 10) +
    stat_bin(aes(y=(..count..)/sum(..count..),
                 label=(..count..)),
             geom="text", size=4, binwidth = 10, vjust=-.7) +
    facet_wrap(~Sex) +
    labs(x = "Age (in years)", y = "Percentage") +
    ggtitle(paste0("Age at recruitment - nonMRI ", toupper(cm))) +
    my_theme()
  ggsave(paste0("Results/figures/Histogram_age_sex_nonMRI_",
                toupper(cm), ".png"),
         plot = p)
  create_pptx(p, PptPath)

  p <- ggplot(df, aes(Age_at_recruitment)) +
    geom_histogram(color = "black", fill = "red3",
                   aes(y = (..count..)/sum(..count..)), binwidth = 10) +
    stat_bin(aes(y=(..count..)/sum(..count..),
                 label=(..count..)),
             geom="text", size=4, binwidth = 10, vjust=-.7) +
    facet_wrap(~Sex) +
    labs(x = "Age (in years)", y = "Percentage") +
    ggtitle(paste0("Age at recruitment - all ", toupper(cm))) +
    my_theme()
  ggsave(paste0("Results/figures/Histogram_age_sex_all_",
                toupper(cm), ".png"),
         plot = p)
  create_pptx(p, PptPath)

  df$CM <- toupper(names(dfs)[q])
  cases <- rbind(cases, df)
} # end for-loop iteration CMs

cases$CM <- as.factor(cases$CM)

sub <- list()
sub$m_no <- subset(cases, Sex == "Male" & MRI == "No")
sub$f_no <- subset(cases, Sex == "Female" & MRI == "No")
sub$m_yes <- subset(cases, Sex == "Male" & MRI == "Yes")
sub$f_yes <- subset(cases, Sex == "Female" & MRI == "Yes")

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

cids <- as.data.frame(controls$f.eid)
write.table(cids, paste0(path, "/Control_IDs.txt"), col.names = "f.eid",
            row.names = FALSE, quote = FALSE)
