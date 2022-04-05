#!/usr/bin/env Rscript
##
## Script information ----------------------------------------------------------
##
## Script name: Create_Figure4_S2_S3.R
##
## Purpose of script: Create the figures 4, S2 and S3
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
pix <- 0.393700787


## Loading packages ------------------------------------------------------------

suppressMessages(library(ggpubr))
suppressMessages(library(forestplot))
suppressMessages(library(VennDiagram))
suppressMessages(library(svglite))
source("src/functions.R")


## Loading arguments -----------------------------------------------------------

# Arguments expected:
#     #1 -- Path to and name of the file with CM-stratified differences
#           for example /hpc/dhl_ec/mvanvugt/UKBB/results/output/TableS6_7.tsv
#     #2 -- Path to and name of the file with gene-stratified differences,
#           for example /hpc/dhl_ec/mvanvugt/UKBB/results/output/TableS9.tsv
#     #3 -- Path to figure directory,
#           for example /hpc/dhl_ec/mvanvugt/UKBB/results/figures

args   = commandArgs(trailingOnly = TRUE)
cmdif  = args[1]
gendif = args[2]
figs   = args[3]


## General options -------------------------------------------------------------

pheno <- c("Outcomes", "Heart_Failure", "Cardiomyopathy", "HFCM", "Pheno",
           "Ventricular_arrhythmias", "Atrial_arrhythmias", "Heart_Arrhythmia", 
           "Chronic_ischaemic_heart_disease", "Angina", 
           "Cardiovascular_Death", "All_cause_mortality")


##  Figure S2 - all CMs --------------------------------------------------------

# Read data
message(paste0("Read data from ", cmdif))
dat <- read.delim(cmdif) %>%
    filter(CM %in% c("ARVC G+ vs G-", "DCM G+ vs G-", "HCM G+ vs G-", "strict HCM G+ vs G-")) %>%
    filter(Phenotype %in% pheno) %>% arrange(factor(Phenotype, levels = pheno))

# Some label changes
dat$CM <- gsub(" vs G-", "", dat$CM)
dat$CM[dat$CM == "strict HCM G+"] <- "HCM* G+"

# Some preparations
dat[, 2:4] <- lapply(dat[, 2:4], as.numeric)
dat <- dat[order(dat[, 1], dat[, 7]),] %>% 
    arrange(factor(Phenotype, levels = pheno)) 
names(dat) <- c("Phenotype", "OR", "LCI", "UCI", "pvalue", "test", "CM")

# Create dataframe with text for forestplot
message("Creating Figure S2")
text <- data.frame()
# Iterate over the phenotypes
for (i in pheno) {

    # Create the header
    if (i == "Outcomes") {
        t <- c(i, "OR (95% CI)")
    } else {
        temp <- dat %>% filter(Phenotype == i)
        for (g in 1:nrow(temp)) {
            # Paste together the OR and 95% CI into text
            s <- paste0(format(round(temp[g, 2], 2), nsmall = 2), " (", format(round(temp[g, 3], 2), nsmall = 2), ";", format(round(temp[g, 4], 2), nsmall = 2), ")")
            if (g == 1) {
                # Create variable 
                p <- s
            } else {
                # Paste with newline to variable p
                p <- paste0(p, "\n", s)
            } # end if first gene
            t <- c(i, p)
        } # iterate over CMs  
    }
    text <- rbind(text, t)
    text[1:2] <- lapply(text[1:2], as.character)
}

# Some more textual changes
text[, 1] <- gsub("_", " ", text[, 1])
text[,1][text[, 1] == "Heart Failure"] <- "Heart failure"
text[,1][text[, 1] == "HFCM"] <- "Heart failure + cardiomyopathy"
text[,1][text[, 1] == "Pheno"] <- "Phenotype positive"
text[,1][text[, 1] == "Angina"] <- "Angina pectoris"
text[,1][text[, 1] == "Chronic ischaemic heart disease"] <- "Chronic ischemic heart disease"
text[,1][text[, 1] == "Heart Arrhythmia"] <- "Self-reported heart arrhythmias"
text[,1][text[, 1] == "Cardiovascular Death"] <- "Cardiovascular death"
text[,1][text[, 1] == "All cause mortality"] <- "All-cause mortality"

# Create and save figure S2
message(paste0("Saving Figure S2 in ", file.path(figs, "FigureS2.pdf")))
pdf(file.path(figs, "FigureS2.pdf"), paper = "a4", width = 10 * pix, 
    height = 12 * pix)
dat %>% group_by(CM) %>% 
    forestplot(labeltext = text, mean = OR, lower = LCI, upper = UCI, 
               zero = 1, fn.ci_norm = c(fpDrawNormalCI, fpDrawDiamondCI, 
                                        fpDrawCircleCI, fpDrawCircleCI),
               shapes_gp = fpShapesGp(box = c("#ffd167", "#168ab2", "#ef476f", "black") %>% 
                                      lapply(function(x) gpar(fill = x, col = x)),
                                  line = c("#ffd167", "#168ab2", "#ef476f", "black") %>% 
                                      lapply(function(x) gpar(fill = x, col = x)),
                                  default = gpar(vertices = TRUE)),
               boxsize = .1, graph.pos = 2, align = "r", lwd.xaxis = .5,
               lwd.zero = .5, lwd.ci = .5, xlog = TRUE,
               txt_gp = fpTxtGp(label = gpar(cex = .3), xlab = gpar(cex = .35),
                                ticks = gpar(cex= .3), legend = gpar(cex = .3)),
               is.summary = c(TRUE, rep(FALSE, 24)), xlab = "Odds ratio (95% CI)")
dev.off()


## Figure 4 - Only CMs ---------------------------------------------------------

# Filter out strict HCM group
dat <- dat %>% filter(CM != "HCM* G+")

message("Creating Figure 4")
# Create text df for forestplot
text <- data.frame()
# Iterate over phenotypes
for (i in pheno) {

    # Create header
    if (i == "Outcomes") {
        t <- c(i, "OR (95% CI)")
    } else {
        temp <- dat %>% filter(Phenotype == i)
        for (g in 1:nrow(temp)) {
            # Paste together the OR and 95% CI into text
            s <- paste0(format(round(temp[g, 2], 2), nsmall = 2), " (", format(round(temp[g, 3], 2), nsmall = 2), ";", format(round(temp[g, 4], 2), nsmall = 2), ")")
            if (g == 1) {
                # Create variable 
                p <- s
            } else {
                # Paste with newline to variable p
                p <- paste0(p, "\n", s)
            } # end if first gene
            t <- c(i, p)
        } # iterate over CMs  
    }
    text <- rbind(text, t)
    text[1:2] <- lapply(text[1:2], as.character)
}

# Some textual changes 
text[, 1] <- gsub("_", " ", text[, 1])
text[,1][text[, 1] == "Heart Failure"] <- "Heart failure"
text[,1][text[, 1] == "HFCM"] <- "Heart failure + cardiomyopathy"
text[,1][text[, 1] == "Pheno"] <- "Phenotype positive"
text[,1][text[, 1] == "Angina"] <- "Angina pectoris"
text[,1][text[, 1] == "Chronic ischaemic heart disease"] <- "Chronic ischemic heart disease"
text[,1][text[, 1] == "Heart Arrhythmia"] <- "Self-reported heart arrhythmias"
text[,1][text[, 1] == "Cardiovascular Death"] <- "Cardiovascular death"
text[,1][text[, 1] == "All cause mortality"] <- "All-cause mortality"

# Create and save forestplot Figure 4
message(paste0("Saving Figure S4 in ", file.path(figs, "FigureS4.pdf")))
pdf(file.path(figs, "Figure4.pdf"), paper = "a4", width = 10 * pix, 
              height = 9 * pix)
dat %>% group_by(CM) %>% 
    forestplot(labeltext = text, mean = OR, lower = LCI, upper = UCI, zero = 1, 
               fn.ci_norm = c(fpDrawNormalCI, fpDrawDiamondCI, fpDrawCircleCI),
               shapes_gp = fpShapesGp(box = c("#ffd167", "#168ab2", "#ef476f") %>% 
                                      lapply(function(x) gpar(fill = x, col = x)),
                                  line = c("#ffd167", "#168ab2", "#ef476f") %>% 
                                      lapply(function(x) gpar(fill = x, col = x)),
                                  default = gpar(vertices = TRUE)),
               boxsize = .1, graph.pos = 2, align = "r", lwd.xaxis = .5,
               lwd.zero = .5, lwd.ci = .5, xlog = TRUE,
               txt_gp = fpTxtGp(label = gpar(cex = .3), xlab = gpar(cex = .35),
                                ticks = gpar(cex= .3), legend = gpar(cex = .3)),
               is.summary = c(TRUE, rep(FALSE, 18)), xlab = "Odds ratio (95% CI)")
dev.off()
rm(dat, temp, text, i, t)


## Figure S3 - Stratified by genes ---------------------------------------------

data <- read.delim(gendif) %>%
    filter(Phenotype %in% pheno) %>% arrange(factor(Phenotype, levels = pheno)) 
data$CM <- gsub(" vs G-", "", data$CM)

# Iterate over CMs
message("Creating figure S3")
for (cm in levels(as.factor(data$CM))) {

    # Select only correct genes for this CM
    if (cm == "ARVC G+") {
        dat <- data %>% filter(Gene %in% c("DES", "PKP2"))
        # define necessary height of final pdf
        h <- 8
    } else if (cm == "DCM G+") {
        dat <- data %>% filter(Gene %in% c("BAG3", "DES", "FLNC", "LMNA", 
                                           "MYH7", "TTN"))
        # define necessary height of final pdf
        h <- 18
    } else if (cm == "HCM G+") {
        dat <- data %>% filter(Gene %in% c("JPH2", "MYBPC3", "MYH7"))
        # define necessary height of final pdf
        h <- 12
    }

    # Include some extra rows for the forest plot header
    ex <- unique(subset(dat, CM == cm)$Gene)
    ex <- data.frame(Phenotype = "Outcomes", OR = NA, new = NA, old = NA, 
                     pvalue = NA, CM = cm, Gene = ex)
    names(ex) <- names(dat)
    dat <- rbind(dat, ex) %>% arrange(factor(Phenotype, levels = pheno))

    # Clean up df
    dat[, 2:4] <- lapply(dat[, 2:4], as.numeric)
    dat <- dat[order(dat[, 1], dat[, 6]),] %>% arrange(factor(Phenotype, levels = pheno)) 
    names(dat) <- c("Phenotype", "mean", "lower", "upper", "pvalue", "CM", "Gene")

    # Create text df for forestplot
    text <- data.frame()
    for (i in pheno) {

        # Create header
        if (i == "Outcomes") {
            t <- c(i, "OR (95% CI)")
        } else {
            temp <- dat %>% filter(Phenotype == i) %>% filter(CM == cm)
            for (g in 1:nrow(temp)) {
                # Paste together the OR and 95% CI into text
                s <- paste0(format(round(temp[g, 2], 2), nsmall = 2), " (", format(round(temp[g, 3], 2), nsmall = 2), ";", format(round(temp[g, 4], 2), nsmall = 2), ")")
                if (g == 1) {
                    # Create variable 
                    p <- s
                } else {
                    # Paste with newline to variable p
                    p <- paste0(p, "\n", s)
                } # end if first gene
                t <- c(i, p)
            } # iterate over genes
        } # end if phenotype is outcome or not
        text <- rbind(text, t)
        text[1:2] <- lapply(text[1:2], as.character)
    } # Iterate over phenotypes

    # Some textual changes
    text[, 1] <- gsub("_", " ", text[, 1])
    text[,1][text[, 1] == "Heart Failure"] <- "Heart failure"
    text[,1][text[, 1] == "HFCM"] <- "Heart failure + cardiomyopathy"
    text[,1][text[, 1] == "Pheno"] <- "Phenotype positive"
    text[,1][text[, 1] == "Angina"] <- "Angina pectoris"
    text[,1][text[, 1] == "Chronic ischaemic heart disease"] <- "Chronic ischemic heart disease"
    text[,1][text[, 1] == "Heart Arrhythmia"] <- "Self-reported heart arrhythmias"
    text[,1][text[, 1] == "Cardiovascular Death"] <- "Cardiovascular death"

    # Remove strange values that raise errors in plotting
    dat[dat == Inf] <- NA
    dat[dat == 0] <- NA

    # prepare some settings for the forestplot function
    shapes <- c(fpDrawNormalCI, fpDrawDiamondCI, fpDrawCircleCI)
    colors <- c("#519872", "#ffd167", "#168ab2", "#ef476f", "black", "#034C3C", "grey")

    message(paste0("Saving figure in ", file.path(figs, "FigureS3_"), cm, ".pdf"))
    pdf(paste0(file.path(figs, "FigureS3_"), cm, ".pdf"), paper = "a4",
        width = 10 * pix, height = h * pix)
    plot <- dat %>% filter(CM == cm) %>% group_by(Gene) %>% 
        forestplot(labeltext = text, zero = 1, 
                   fn.ci_norm = rep(shapes, ceiling(nrow(temp)/3))[1:nrow(temp)],
                   shapes_gp = fpShapesGp(box = colors[1:nrow(temp)] %>%
                                          lapply(function(x) gpar(fill = x, col = x)),
                                      line = colors[1:nrow(temp)] %>%
                                          lapply(function(x) gpar(fill = x, col = x)),
                                      default = gpar(vertices = TRUE)),
                   boxsize = .1, graph.pos = 2, align = "r", lwd.xaxis = .5,
                   lwd.zero = .5, lwd.ci = .5, xlog = TRUE, 
                   txt_gp = fpTxtGp(label = gpar(cex = .3), xlab = gpar(cex = .35),
                                    ticks = gpar(cex= .3), legend = gpar(cex = .3)),
                   is.summary = c(TRUE, rep(FALSE, 24)), xlab = "Odds ratio (95% CI)")
    print(plot)
    dev.off()

} # iterate over CMs

