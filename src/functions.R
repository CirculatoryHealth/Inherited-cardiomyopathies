## Script information ---------------------------
##
## Script name: functions.R
##
## Purpose of script: Define functions necessary for this project
##
## Author: M. van Vugt
##
## Date Created: 2021-07-16
##
## Copyright (c) M. van Vugt, 2021
## Email: m.vanvugt-2@umcutrecht.nl
##

## Options ---------------------------

options(scipen = 6, digits = 4) # view outputs in non-scientific notation


## Loading packages ---------------------------

library(dplyr)
library(ggplot2)
library("officer")
library("rvg")


# Start functions ---------------------------------------------------------

my_theme <- function() {
  theme(
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.key.size = unit(.2, "cm"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.line = element_line(colour = "black"),
    panel.grid.major.y = element_line(colour = "gray90"),
    strip.background = element_rect(fill = "transparent", colour = "black"),
    text = element_text(size = 4),
    line = element_line(size = .2),
  )
}

pc_theme <- function() {
  theme(
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.line = element_line(colour = "black"),
    panel.grid.major.y = element_line(colour = "gray90"),
    strip.background = element_rect(fill = "transparent", colour = "black"),
    text = element_text(size = 15)
  )
}


perc_var <- function(full, vars, desc = vars) {
  #' Summarizing counts and percentages of factor levels.
  #' 
  #' @description This function summarizes the numbers and frequencies of 
  #' appearances per level of all specified factors.
  #' 
  #' @param full Dataframe containing the factors.
  #' @param vars Vector of all column names that should be treated. Should 
  #' be factors.
  #' @param desc Vector with description of all column names specified in vars.
  #' Should be same length and in the same order as vars.
  #' 
  #' @return A dataframe with four columns (value, count, perc, name). Value 
  #' is the level of the factor, count is the number of times this level 
  #' occurred for this factor, percentage is the percentage of this level for 
  #' this factor and name is the description as specified by the user in desc.
  # if (all(vars != colnames(full))) {
  #   temp <- full %>% dplyr::select(all_of(vars))
  # } else {
  #   temp <- full
  # }
  list <- list()
  for (i in 1:length(vars)) {
    temp1 <- full %>%
      group_by(get(vars[i])) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(perc = (count / sum(count)) * 100)
    temp1$name <- desc[i]
    colnames(temp1) <- c("value", "count", "perc", "name")
    list[i] <- list(temp1)
  }
  df <- do.call("rbind", list)
}


# Save figures to PPT -----------------------------------------------------

create_pptx <- function(plt = last_plot(), path = file.choose()) {
  if (!file.exists(path)) {
    out <- read_pptx()
  } else {
    out <- read_pptx(path)
  }
  
  out %>%
    add_slide(layout = "4_Aangepaste indeling", master = "8_Office-thema") %>%
    ph_with(value = dml(ggobj = plt), location = ph_location_fullsize()) %>%
    print(target = path)
}

