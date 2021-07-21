library(plotly)

m <- list(l = 50, r = 50, b = 300, t = 100, pad = 4)

pie <- plot_ly(tmp, 
               labels = paste0(tmp$Gene, "\n", tmp$N, " (", round(tmp$prop), "%)"),
               values = ~N, type = "pie", showlegend = FALSE,
               textposition = ifelse(tmp$prop < 4, "outside", "inside"), 
               textinfo = "label",
               insidetextfont = list(color = "black"),
               marker = list(colors = viridis(nrow(tmp), alpha = 0.8),
                             line = list(color = "#FFFFFF", width = 1)))
pie <- pie %>% layout(title = paste0("Mutated genes in ", cm),
                      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                      autosize = FALSE, margin = m)
pie
