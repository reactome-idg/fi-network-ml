read.file <- function(file.name) {
    cor.values <- read.delim(file.name, header = T, sep = "\t")
    names(cor.values) <- c("GenePair", "Downloaded", "Calculated")
    cor.values
}

compare.correlations <- function(cor.values, title = "Correlation") {
    t.test.results <- t.test(cor.values$Downloaded, cor.values$Calculated)
    print(t.test.results)
    boxplot(cor.values[2 : length(cor.values)], xlab = title)
}

compare.two.correlations <- function(cor.values1, cor.values2) {
    t.test.results <- t.test(cor.values1$Download, cor.values2$Download)
    print("T-test for two downloded results:")
    print(t.test.results)
    t.test.results <- t.test(cor.values1$Calculated, cor.values2$Calculated)
    print("T-test for two calculated results:")
    print(t.test.results)
}

dir <- "/Users/wug/Documents/wgm/work/reactome-idg/archs4"
file.name1 <- "correlations_comparison.tsv"
file.name2 <- "correlations_comparison_time_filtered.tsv"
cor.values.all <- read.file(paste(dir, file.name1, sep = "/"))
cor.values.filtered <- read.file(paste(dir, file.name2, sep = "/"))

# Plot two results together
par(mfrow = c(1, 2))
print("Check all correlations:")
compare.correlations(cor.values.all, "Without Filtering")
print("Check filtered correlations:")
compare.correlations(cor.values.filtered, "With Filtering")

# Compare two downloaded values
print("Compare the same correlations:")
compare.two.correlations(cor.values.all, cor.values.filtered)


