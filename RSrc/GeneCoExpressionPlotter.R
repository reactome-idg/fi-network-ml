library("ggplot2")
# Get parameters from the command
args <- commandArgs(TRUE);
if (length(args) < 3) {
    stop("Three arguments (coexpression.file.name, out.file.name, title) should be provided in order to run this script!");
}
# dir.name <- "/Users/wug/git/reactome-idg/idg-pairwise/examples"
# file.name <- paste(dir.name, "test.txt", sep = "/")
# tile <- "test"
file.name <- args[1]
out.file.name <- args[2]
title <- args[3]
data <- read.delim(file.name, header = TRUE)
plot <- ggplot(data, aes(x = Coexpression, y = ..density..)) + 
        geom_histogram(color = "grey", fill = "lightgrey", binwidth = 0.001) +
        geom_density(color = "blue") +
        labs(title = title) +
        theme(plot.title = element_text(hjust = 0.50))
# print(plot)
ggsave(out.file.name, plot = plot)
