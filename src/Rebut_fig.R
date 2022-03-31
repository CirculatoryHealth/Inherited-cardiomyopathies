library(dplyr)
library(ggplot2)
library(ggpubr)
library(forestplot)
source("src/functions.R")
pix <- 0.393700787


# Forest plot -------------------------------------------------------------

pheno <- c("Outcomes", "Heart_Failure", "Cardiomyopathy", "HFCM", "Pheno",
           "Ventricular_arrhythmias", "Atrial_arrhythmias", "Heart_Arrhythmia", 
           "Chronic_ischaemic_heart_disease", "Angina", 
           "Cardiovascular_Death", "All_cause_mortality")
data <- read.delim("results/output/Reb_Differences_genes.tsv") %>%
    filter(Phenotype %in% pheno) %>% arrange(factor(Phenotype, levels = pheno)) 
# dat <- read.delim("~/surfdrive/Mendelian_CM_WES_UKB/Difference_testing_v4.tsv") %>%
#   filter(CM %in% c("ACM", "DCM", "HCM", "strict HCM")) %>% 
#   filter(Phenotype %in% pheno) %>% arrange(factor(Phenotype, levels = pheno)) 
# dat$CM <- paste0(dat$CM, " G+")
# dat$CM[dat$CM == "strict HCM G+"] <- "HCM* G+"

for (cm in c("ACM", "DCM", "HCM")) {

    # Select only correct genes for this CM
    if (cm == "ACM") {
        dat <- data %>% filter(Gene %in% c("DES", "PKP2"))
        # define necessary height of final pdf
        h <- 8
    } else if (cm == "DCM") {
        dat <- data %>% filter(Gene %in% c("BAG3", "DES", "FLNC", "LMNA", 
                                          "MYH7", "TTN"))
        # define necessary height of final pdf
        h <- 18
    } else if (cm == "HCM") {
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

    dat[, 2:4] <- lapply(dat[, 2:4], as.numeric)
    dat <- dat[order(dat[, 1], dat[, 6]),] %>% arrange(factor(Phenotype, levels = pheno)) 
    names(dat) <- c("Phenotype", "mean", "lower", "upper", "pvalue", "CM", "Gene")

    text <- data.frame()
    for (i in pheno) {

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

    pdf(paste0("results/figures/FigureS2_", cm, ".pdf"), paper = "a4",
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

