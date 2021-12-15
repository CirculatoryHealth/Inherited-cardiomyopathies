#!/usr/bin/env Rscript

############################################################
# Loading packages
############################################################

library(argparse)
library(fs)

############################################################
# Argparse
############################################################

# create parser object
parser <- ArgumentParser()

parser$add_argument('--main',
    type='character',
    help='Path to main dir with .csv files from the extraction_results.py, optional mapper files with columns: sel, xlim, at, xlab, title',
    nargs='?',
    default='current'
)

parser$add_argument('--coderepos',
    type='character',
    help='Path to the coderepos dir',
    nargs='?',
    default='/home/amand/google_drive/Research/CodeRepos/'
)

# getting arguments
args <- parser$parse_args()

############################################################
# PATHS
############################################################

# getting root dir
argv <- commandArgs(trailingOnly = FALSE)[4]
base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
base_dir <- path_dir(base_dir)
print(base_dir)
# base_dir = args$coderepos

# sourcing custom stuff
source(paste0(base_dir,'/R/PointAndCI.R'))
source(paste0(base_dir,'/R/Preamble_standard.R')) # need a more specific file
source(paste0(base_dir,'/R/MR_forestplots.R'))

# normal packages
library(dplyr)
library(forestplot)

############################################################
# LOADING DATA
############################################################

# for debugging
# args$main = '/home/amand/google_drive/Research/UMCU/Cupido/project_cupido/results_100kb'
if(args$main == 'current'){
    wd <- getwd()
} else{
    wd <- args$main
}
print(paste0('The current wd: ', wd))

# the relevant data
setwd(wd)
content <- list.files(pattern = '.csv')
# optinal selection
# content <- content[!grepl('_all.csv', content)]

# find mapper data, dealing with quotes
mapper <- read.csv(paste0(wd, '/mapper.txt'), sep ='\t')
if(!any(colnames(mapper) %in% 'ensemblid')){
    mapper <- read.csv(paste0(wd, '/mapper.txt'), sep ='\t', quote="'")
}

############################################################
# CONSTANTS
############################################################

# alpha
a <- 0.05
# figures sizes
w <- 16 * 0.393700787
h <- 14 * 0.393700787
# other
h_ci <- 0
size_ci <- 0.2

############################################################
# Generating the figures
############################################################

for(i in content){

    results <- read.csv(i)
    # out filename
    outfile <- sub('.csv','.pdf',i)
    split1 <- unlist(strsplit(i, '_|\\.'))
    ensemblid <- split1[1]
    gene <- split1[2]
    # file specific settings
    settings = mapper[mapper$ensemblid == ensemblid,]
    # labels
    if( any(colnames(settings) %in% 'xlab') ){
        xlab1 <- as.character(settings$xlab)
        title1 <- as.character(settings$title)
    } else {
        xlab1 <- 'Increasing exposure'
        title1 <- gene
    }
    # selecting columns
    if( any(colnames(settings) %in% 'noncvd')){
        r_sel2 <- as.character(settings$noncvd)
        r_sel2 <- unlist(strsplit(r_sel2, "[|]"))
    } else {
        r_sel2 <- c('dmmr', 'ckdmr', 'asthmamr', 'ibsmr', 'chronmr', 'ucmr','msmr', 'alzheimersmr')
    }
    if( any(colnames(settings) %in% 'cvd')){
        r_sel1 <- as.character(settings$cvd)
        r_sel1 <- unlist(strsplit(r_sel1, "[|]"))
    } else {
        r_sel1 <- c('chdmr', 'asmr', 'aismr',
                'lasmr', 'cesmr', 'svsmr', 'hfmr', 'afmr')
    }
    # limits
    if( any(colnames(settings) %in% 'xlimstart')){
        xlim1 <- c(settings$xlimstart, settings$xlimend)
        at1 <- c(settings$xlimstart, settings$middle, settings$xlimend)
    } else {
        xlim1 = c(0.5, 2); at1 = c(0.5, 1.0, 2.0)
    }
    if(nrow(results) == 0 | mean(is.na(results$point)) == 1){
    pdf(paste0(outfile,'_empty'), width =w , height = h)
    dev.off()
    # skipp next part
    next
    }
    # if not empyt do the following
    # sub-setting
    if(any(colnames(results) %in% 'outcome.1')){ results$selrows <- as.character(results$outcome); results$outcome <- as.character(results$outcome.1)
    }else if(any(colnames(results) %in% 'outcome.2')){ results$selrows <- as.character(results$outcome.1); results$outcome <- as.character(results$outcome.2)
    } else {
        results$selrows <- as.character(results[,1])
        results$outcome <- as.character(results$outcome)
    }
    bin <- results[results$unit == 'logOR',]
    input <- bin[,c('outcome', 'selrows', 'point','se','pvalue','no_cases', 'samplesize')]
    # rownames(input) <- as.character(input$outcome)
    # scaling [optional]
    if ( any(colnames(settings) %in% 'scale')){
        scale <- as.numeric(settings$scale)
    } else{
            scale <- 1}
    input$point <- scale * input$point
    input$se2 <- abs(scale) * input$se
    input$lb <- input$point - qnorm(1-a/2) * input$se2
    input$ub <- input$point + qnorm(1-a/2) * input$se2

    # order tables
    fig1_p <- input[input$selrows %in%  r_sel1,]
    fig1_p <- fig1_p[match(r_sel1,fig1_p$selrows),]
    fig2_p <- input[input$selrows %in%  r_sel2 ,]
    fig2_p <- fig2_p[match(r_sel2,fig2_p$selrows),]

    # plots
    pdf(outfile, width =w , height = h, compress = F)
    forest(format(fig1_p,fig2_p), xlab = xlab1, xlim = xlim1, at = at1)
    grid.text(title1, just = 'left', .03, .95, gp=gpar(cex=1.5, fontface='bold'))
    dev.off()

}

