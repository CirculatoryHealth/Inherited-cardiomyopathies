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

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(officer))
suppressMessages(library(rvg))
suppressMessages(library(docstring))


# ggplot2 themes ----------------------------------------------------------

# theme ready for publication
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

# theme for presentations (screens)
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


# Summarizing function ----------------------------------------------------

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


# Create incidence matrix graph -------------------------------------------

inc_mat <- function(data, sig = 0.05, sig1 = NULL, pdif = "color", 
                    xas = names(data)[ncol(data)], legend = "none", 
                    odif = NULL, oname = NULL, cat = NULL) {

    #' Create incidence matrix graphs
    #' 
    #' @description This function creates an incidence matrix graph, showing 
    #' whether tests were significant by displaying different colors
    #' 
    #' @param data Dataframe with header of all exposures tested and p-values 
    #' in the dataframe. One of the columns should be the outcomes, column 
    #' name should be specified by `xas`
    #' @param sig Numerical, significance level for which p-values below this 
    #' level will have lay-out indicating significance. Default is 0.05
    #' @param sig1 Either NULL or numerical Second significance level, should 
    #' be lower than `sig`. Lay-out will be different for the values. 
    #' Default is NULL 
    #' @param pdif Character, type of differentiation for the p-value, choose 
    #' from ("color", "colour", "shape")
    #' @param xas Character, name of the column in which outcomes for all tests
    #' are given. Default is last column name of `data`
    #' @param legend Character, the position of legends ("none", "left", 
    #' "right", "bottom", "top", or two-element numeric vector). 
    #' Default is "none"
    #' @param odif Character, type of differentiation for other parameter, 
    #' which is displayed by the p-value being positive or negative
    #' choose from ("color", "colour", "shape"). Default is NULL
    #' @param oname Character vector, containing the name of the legend for 
    #' the other parameter, the name of a negative pvalue and a positive one.
    #' Default is NULL
    #' @param cat Either NULL or a character vector, in latter case with length
    #' equal to number of exposures, indicating categories corresponding to the
    #' categories of the exposures. These will be displayed on the right side of 
    #' the graph and seperated by horizontal lines. Default is NULL
    #' 
    #' @return A ggplot object with the incidence matrix graph

    ### Check parameters ###
    if (!all(c(odif, pdif) %in% c("color", "colour", "shape"))) {
        stop("Argument pdif is invalid, should be color or shape")
    }
    if (!xas %in% names(data)) {
        stop("Argument xas is invalid, should be a column name of data")
    }
    if (odif == pdif) {
        stop("Arguments odif and pdif should be different!")
    }


    ### Prepare data ###
    data <- as_tibble(data)

    # Get number of outcomes and exposures
    row <- nrow(unique(data[, xas]))
    new <- data %>% select(-any_of(xas))
    col <- ncol(new)

    # Create grid for graph
    xy <- data.frame(x = rep(c(1:row), col),
                     y = sort(rep(c(1:col), row)),
                     p = unlist(c(new[, 1:ncol(new)])))

    # Create p-value categories
    xy$pcat[abs(xy$p) >= sig] <- paste0("p > ", sig)
    xy$pcat[abs(xy$p) < sig] <- paste0("p < ", sig)
    if (!is.null(sig1)) {
        # Check pvalue significance
        if (sig1 >= sig) stop("sig1 should be smaller than sig!")
        # Correct format of significance value
        if (sig1 < 0.001) sigf <- format(sig1, scientific = TRUE)
        xy$pcat[abs(xy$p) < sig1] <- paste0("p < ", sigf)
        lev <- c(paste0("p > ", sig), paste0("p < ", sig), paste0("p < ", sigf))

        # Define colors or shapes depending on pvalue
        if (pdif %in% c("color", "colour")) {
            cols <- c("gray50", "orange2", "red3")
            siz = c(0.1, 0.25, 0.35)
        } else {
            shap <- c(17, 16, 15)
            siz = c(0.4, 0.8, 0.9)
        }
    } else {
        lev <- c(paste0("p > ", sig), paste0("p < ", sig))
        # Define colors or shapes depending on pvalue
        if (pdif %in% c("color", "colour")) {
            cols <- c("gray50", "red3")
            size = c(0.1, 0.3)
        } else {
            shap <- c(17, 16)
            siz = c(0.4, 1)
        }
    }
    xy$pcat <- as.factor(xy$pcat)

    if (!is.null(odif)) {
        # Define categories
        xy$change[xy$p < 0] <- oname[2]
        xy$change[xy$p >= 0] <- oname[3]
        change <- unique(xy$change)
        xy$change <- as.factor(xy$change)

        # Define colors or shapes depending on categories
        if (odif == "shape") {
            shap <- c(17, 16)
        } else {
            cols <- c("red3", "green3")
        }
    }

    # Add categories
    if (!is.null(cat)) {
        if (length(cat) != ncol(new)) stop("Number of categories and exposures are different, please check!")
        ycat <- c()
        ytxt = c()
        ycum = 0
        for (c in unique(cat)) {
            ycat <- c(ycat, length(grep(c, cat)) + 0.5 + ycum)
            ytxt <- c(ytxt, (ycum + ycat[length(ycat)]) / 2)
            ycum = ycum + length(grep(c, cat))
        }
        xtxt = row + 1
        ycat <- ycat[1:(length(ycat) - 1)]
        txt = unique(cat)
    } else {
        ycat = 0
        xtxt = 0
        ytxt = 0
        txt = ""
    }

    ### Make plot ###
    if (is.null(odif) & pdif %in% c("color", "colour")) {
        p <- ggplot(xy, aes(x = x, y = y, color = pcat, size = pcat)) +
            geom_point() +
            geom_hline(yintercept = ycat, size = .2) +
            annotate(geom = "text", x = xtxt, y = ytxt, label = txt, 
                     angle = 90, size = 1, fontface = 2) +j
            scale_x_continuous(breaks = unique(xy$x), 
                               labels = pull(data, xas), 
                               name = NULL, expand = c(0.05, 0.05)) +
            scale_y_continuous(breaks = unique(xy$y), labels = names(new), 
                               name = NULL, expand = c(0.005, 0.005)) +
            my_theme() +
            theme(panel.grid.major.x = element_line(colour = "gray90"),
                  axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, 
                                                    hjust=1),
                  axis.line.x.top =  element_line(color = "black"),
                  axis.text.x.top = element_text(angle = -90, vjust = 0.5, 
                                                 hjust=1),
                  axis.line.y.right =  element_line(color = "black"),
                  axis.ticks.y.right = element_blank(),
                  axis.text.y.right = element_blank(),
                  legend.position = legend) +
            guides(x.sec = "axis", y.sec = "axis") +
            scale_color_manual(breaks = lev, values = cols, 
                               na.translate = FALSE, name = "Significance") +
            scale_size_manual(breaks = lev, values = siz, 
                              na.translate = FALSE, name = "Significance") +
            coord_cartesian(xlim = c(1, row), ylim = c(0, col + 1), 
                            clip = "off")
    } else if (is.null(odif) & pdif == "shape") {

        p <- ggplot(xy, aes(x = x, y = y, shape = pcat, size = pcat)) +
            geom_point() +
            geom_hline(yintercept = ycat, size = .2) +
            annotate(geom = "text", x = xtxt, y = ytxt, label = txt, 
                     angle = 90, size = 1, fontface = 2) +
            scale_x_continuous(breaks = unique(xy$x), 
                               labels = pull(data, xas), 
                               name = NULL, expand = c(0.05, 0.05)) +
            scale_y_continuous(breaks = unique(xy$y), labels = names(new), 
                               name = NULL, expand = c(0.005, 0.005)) +
            my_theme() +
            theme(panel.grid.major.x = element_line(colour = "gray90"),
                  axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, 
                                                    hjust=1),
                  axis.line.x.top =  element_line(color = "black"),
                  axis.text.x.top = element_text(angle = -90, vjust = 0.5, 
                                                 hjust=1),
                  axis.line.y.right =  element_line(color = "black"),
                  axis.ticks.y.right = element_blank(),
                  axis.text.y.right = element_blank(),
                  legend.position = legend) +
            guides(x.sec = "axis", y.sec = "axis") +
            scale_shape_manual(breaks = lev, values = cols, 
                               na.translate = FALSE, name = "Significance") +
            scale_size_manual(breaks = lev, values = siz, 
                              na.translate = FALSE, name = "Significance") +
            coord_cartesian(xlim = c(1, row), ylim = c(0, col + 1), 
                            clip = "off")

    } else if (odif == "shape" & pdif %in% c("color", "colour")) {

        p <- ggplot(xy, aes(x = x, y = y, color = pcat, size = pcat, 
                            shape = change)) +
            geom_point() +
            geom_hline(yintercept = ycat, size = .2) +
            annotate(geom = "text", x = xtxt, y = ytxt, label = txt, 
                     angle = 90, size = 1, fontface = 2) +
            scale_x_continuous(breaks = unique(xy$x), 
                               labels = pull(data, xas), 
                               name = NULL, expand = c(0.05, 0.05)) +
            scale_y_continuous(breaks = unique(xy$y), labels = names(new), 
                               name = NULL, expand = c(0.005, 0.005)) +
            my_theme() +
            theme(panel.grid.major.x = element_line(colour = "gray90"),
                  axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, 
                                                    hjust=1),
                  axis.line.x.top =  element_line(color = "black"),
                  axis.text.x.top = element_text(angle = -90, vjust = 0.5, 
                                                 hjust=1),
                  axis.line.y.right =  element_line(color = "black"),
                  axis.ticks.y.right = element_blank(),
                  axis.text.y.right = element_blank(),
                  legend.position = legend) +
            guides(x.sec = "axis", y.sec = "axis") +
            scale_color_manual(breaks = lev, values = cols, 
                               na.translate = FALSE, name = "Significance") +
            scale_size_manual(breaks = lev, values = siz, 
                              na.translate = FALSE, name = "Significance") +
            scale_shape_manual(breaks = change, values = shap, name = oname[1],
                               na.translate = FALSE) +
            coord_cartesian(xlim = c(1, row), ylim = c(0, col + 1), 
                            clip = "off")
    } else if (odif %in% c("color", "colour") & pdif == "shape") {

        p <- ggplot(xy, aes(x = x, y = y, color = change, size = pcat,
                            shape = pcat)) +
            geom_point() +
            geom_hline(yintercept = ycat, size = .2) +
            annotate(geom = "text", x = xtxt, y = ytxt, label = txt, 
                     angle = 90, size = 1, fontface = 2) +
            scale_x_continuous(breaks = unique(xy$x), 
                               labels = pull(data, xas), 
                               name = NULL, expand = c(0.05, 0.05)) +
            scale_y_continuous(breaks = unique(xy$y), labels = names(new), 
                               name = NULL, expand = c(0.005, 0.005)) +
            my_theme() +
            theme(panel.grid.major.x = element_line(colour = "gray90"),
                  axis.text.x.bottom = element_text(angle = 90, vjust = 0.5, 
                                                    hjust=1),
                  axis.line.x.top =  element_line(color = "black"),
                  axis.text.x.top = element_text(angle = -90, vjust = 0.5, 
                                                 hjust=1),
                  axis.line.y.right =  element_line(color = "black"),
                  axis.ticks.y.right = element_blank(),
                  axis.text.y.right = element_blank(),
                  legend.position = legend) +
            guides(x.sec = "axis", y.sec = "axis") +
            scale_color_manual(breaks = change, values = cols, 
                               na.translate = FALSE, name = oname[1]) +
            scale_size_manual(breaks = lev, values = siz, 
                              na.translate = FALSE, name = "Significance") +
            scale_shape_manual(breaks = lev, values = shap, 
                               na.translate = FALSE, name = "Significance") +
            coord_cartesian(xlim = c(1, row), ylim = c(0, col + 1), 
                            clip = "off")
    }

return(p)
    }
