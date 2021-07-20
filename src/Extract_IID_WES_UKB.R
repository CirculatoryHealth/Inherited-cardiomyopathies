args = commandArgs(trailingOnly = TRUE)
input = args[1]
genes = args[2]
vector = args[3]
output = args[4]

cat("Loading the data\n")
df <- read.table(input, header = TRUE, stringsAsFactors = FALSE)
gen <- read.table(genes, header = FALSE, stringsAsFactors = FALSE)
gen <- na.omit(gen)
vec <- read.table(vector, header = FALSE, stringsAsFactors=FALSE)[, 1]

ind <- data.frame(cbind(df$ID, df$LV, df$RV))[1 : (nrow(df) - 1),]
names(ind) <- c("ID", "LV", "RV")
df <- df[, colSums(is.na(df)) == 0]

cat("Starting the loop of genes\n")
for (i in 1:nrow(gen)) {
  x <- gen[i,]
  if (x %in% vec) {
    cat(paste0("\nWorking on ", x, ", gene ", i, " / ", nrow(gen), "\n"))
    df_gen <- df[min(grep(x, df)) : max(grep(x, df))]
    snps <- names(df_gen)
    df_gen <- data.frame(df_gen[1 : (nrow(df_gen) - 1),])
    df_gen <- data.frame(apply(df_gen, 2, as.numeric))
    names(df_gen) <- snps
    df_gen[, x] <- NA
    if (ncol(df_gen) == 2) {
      df_gen[, x] <- df_gen[, 1]
    } else {
      df_gen[, x] <- rowSums(df_gen[, 1 : (ncol(df_gen) - 1)])
    }
    ind <- cbind(ind, df_gen[, x])
    names(ind)[ncol(ind)] <- x
    df_gen <- cbind(ind$ID, df_gen)
    names(df_gen)[1] <- "ID"
    try <- reshape2::melt(df_gen, id.vars = c("ID", x))
    try <- subset(try, value != 0, select = c(ID, variable))
    if (any(duplicated(try$ID))) {
      cat("Duplicates present! Adding an extra column")
      names(try) <- c("ID", paste0(x, "_variant.1"))
      dup <- try$ID[duplicated(try$ID)]
      for (i in dup) {
        sub <- subset(try, ID == i)
        try <- subset(try, ID != i)
        sub$num <- c(1:nrow(sub))
        names(sub)[2] <- paste0(x, "_variant")
        sub <- reshape(sub, idvar = "ID", timevar = "num", direction = "wide")
        if (i == dup[1]) {
          subs <- sub
        } else {
          subs <- merge(subs, sub, all = TRUE, by = names(subs))
        }
      }
    try <- merge(try, subs, all = TRUE, by = names(try))
    } else {
      names(try) <- c("ID", paste0(x, "_variant"))
    }
    ind <- merge(ind, try, all = TRUE, by = "ID")
  } else {
    cat(paste0("\n Gene ", x, " not found among the predefined genes, going to the next\n"))
    next
  }

}


cat("\nCalculating total number of variants\n")
sum <- dplyr::select(ind, !contains("variant"))
num <- names(sum)[4:ncol(sum)]
sum[num] <- lapply(sum[num], as.numeric)
sum$Total_variants <- rowSums(sum[, 4:ncol(sum)], na.rm = TRUE)
ind <- merge(ind, sum[, c(1, ncol(sum))], all= TRUE, by = "ID")
ind <- dplyr::select(ind, c(1:3, contains("variant")))
ind <- subset(ind, Total_variants != 0)

cat("\nWriting output")
write.table(ind, output, col.names = TRUE, row.names = FALSE, quote = FALSE)
