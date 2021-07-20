args = commandArgs(trailingOnly = TRUE)
input = args[1]
genes = args[2]
vector = args[3]
output = args[4]

library(data.table)
library(dplyr)

cat("Loading the data\n")
df <- data.table(read.table(input, header = TRUE, stringsAsFactors = FALSE))
gen <- read.table(genes, header = FALSE, stringsAsFactors = FALSE)
gen <- na.omit(gen)
vec <- read.table(vector, header = FALSE, stringsAsFactors=FALSE)[, 1]

ind <- data.frame(ID = df$ID, LV = df$LV, RV = df$RV)

cat("Starting the loop of genes\n")
for (i in 1:nrow(gen)) {
  x <- gen[i,]
  if (x %in% vec) {
    
    cat(paste0("\nWorking on ", x, ", gene ", i, " / ", nrow(gen), "\n"))
    lok <- names(df)[grep(x, df)]
    ind[, paste0(x, "_variant")] <- as.integer(apply(df[, ..lok], 1, function(r) any(r %in% c(grep("1", r, value = T)))))
    
  } else {
    cat(paste0("\n Gene ", x, " not found among the predefined genes, going to the next\n"))
    next
  }
}

  
ind$Total_variants <- rowSums(ind[, 4:ncol(ind)], na.rm = TRUE)
ind <- slice_head(ind, n = (nrow(ind)-1))

cat("\nWriting output")
write.table(ind, output, col.names = TRUE, row.names = FALSE, quote = FALSE)
