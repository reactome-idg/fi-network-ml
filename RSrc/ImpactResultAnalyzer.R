# The following function is used to draw a heatmap for pathway impact scores and
# genes using package "ComplexHeatmap". The code is based on this web post:
# https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/

# Make sure the package is installed. Otherwise, install it
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
    BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
# For colors
library(circlize)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(ggpubr)

# Need to get the same color scheme for top pathways
set.seed(1234)

plot.heatmap <- function(file.name, 
                         score.name, 
                         gene.family.file,
                         top.pathway.file,
                         need.dev.level = FALSE) {
    df <- read.delim(file.name, header = TRUE, sep = "\t", check.names = FALSE)
    df.tmp <- df[, 2:length(df)]
    row.names(df.tmp) <- df[, 1]
    df.matrix <- as.matrix(df.tmp)
    df.matrix.dim <- dim(df.matrix)
    print(paste("Genes and pathways:", df.matrix.dim[1], df.matrix.dim[2], sep = " "))
    # For protein family annotation
    family.df <- load.gene.families(gene.family.file, df.matrix)
    # df.matrix <- scale(df.matrix)
    # Color genes based on their family
    if (need.dev.level) {
        la <- rowAnnotation(Family = family.df$Family,
                            DevLevel = family.df$DevLevel,
                            col = list(Family = c("GPCR" = "red", "IC" = "green", "Kinase" = "blue"),
                                       DevLevel = c("Tdark" = rgb(244/255, 67/255, 54/255), 
                                                    "Tbio" = rgb(255/255, 178/255, 89/255),
                                                    "Tchem" = rgb(91/255, 192/255, 222/255),
                                                    "Tclin" = rgb(51/255, 122/255, 183/255))))
    }
    else {
        la <- rowAnnotation(Family = family.df$Family,
                            col = list(Family = c("GPCR" = "red", "IC" = "green", "Kinase" = "blue")))
    }
    
    # Color pathways based on their topics
    top.df <- load.top.pathways(top.pathway.file, df.matrix)
    if (is.null(top.colors))
        top.colors <<- create.top.colors(top.pathway.file)
    # View(top.colors)
    ta <- columnAnnotation(Topics = top.df$Topic,
                           col = list(Topics = top.colors))
    # View(top.df)
    # View(colnames(df.matrix))
    score.max = max(df.matrix)
    score.min = min(df.matrix)
    score.mid = score.min + (score.max - score.min) / 10
    mycolors <- colorRamp2(breaks = c(score.min, score.mid, score.max), colors = c("lightgrey", "white", "red"))
    map <- Heatmap(df.matrix,
                   name = score.name,
                   column_title = "Pathways",
                   show_row_names = dim(df.matrix)[1] < 100,
                   show_column_names = dim(df.matrix)[2] < 500,
                   row_title = "Genes",
                   col = mycolors,
                   #col = colorRamp2(c(16, 3, 0), brewer.pal(n=3, name="RdYlBu")),
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 2.5),
                   clustering_method_rows = "ward.D2",
                   clustering_method_columns = "ward.D2",
                   left_annotation = la,
                   top_annotation = ta)
    # Call this in the function context
    print(map)
}

load.top.pathways <- function(top.pathway.file, score.df.matrix) {
    top.df <- read.delim(top.pathway.file, header = T, sep = "\t", check.names = FALSE)
    which <- top.df$Pathway %in% colnames(score.df.matrix)
    top.df <- top.df[which, ]
    # Need to sort it 
    order <- match(colnames(score.df.matrix), top.df$Pathway)
    top.df <- top.df[order, ]
    return (top.df)
}

load.gene.families <- function(gene.family.file, score.df.matrix) {
    # For protein family annotation
    family.df <- read.delim(gene.family.file, header = TRUE, sep = "\t", check.names = FALSE)
    which <- family.df$GeneName %in% rownames(score.df.matrix)
    family.df <- family.df[which, ]
    order <- match(rownames(score.df.matrix), family.df$GeneName)
    family.df <- family.df[order, ]
    return(family.df)
}

create.top.colors <- function(top.pathway.file) {
    # Color pathways based on their topics
    top.df <- read.delim(top.pathway.file, header = T, sep = "\t", check.names = FALSE)
    # Create colors
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    topics <- unique(top.df$Topic)
    top.colors <- sample(col_vector, length(topics))
    names(top.colors) <- topics
    return (top.colors)
}

plot.gene.pathways <- function(score.file.name, gene) {
    plot.genes.pathways(score.file.name,
                        "Gene",
                        gene,
                        "PathwayName",
                        "Reactome Pathway")
}

plot.pathway.genes <- function(score.file.name, pathway, dark.target.only = FALSE) {
    plot.genes.pathways(score.file.name, 
                        "PathwayName",
                        pathway,
                        "Gene",
                        "Gene",
                        dark.target.only)
}

plot.genes.pathways <- function(score.file.name, 
                                col.name,
                                key,
                                x.col.name,
                                x.lab,
                                dark.target.only = FALSE) {
    df <- read.delim(score.file.name, header = TRUE, sep = "\t", check.names = FALSE)
    key.which <- df[, col.name] %in% key
    key.df <- df[key.which, ]
    # Process FDR
    key.df$FDR <- -log10(key.df$FDR)
    key.df$FDR <- (key.df$FDR - min(key.df$FDR)) / (max(key.df$FDR - min(key.df$FDR)))
    key.df$Average_Activation <- (key.df$Average_Activation - min(key.df$Average_Activation)) / (max(key.df$Average_Activation - min(key.df$Average_Activation)))
    key.df$Average_Inhibition <- (key.df$Average_Inhibition - min(key.df$Average_Inhibition)) / (max(key.df$Average_Inhibition - min(key.df$Average_Inhibition)))
    key.df$Sum_Score <- key.df$FDR + key.df$Average_Activation + key.df$Average_Inhibition
    # Do another filter here is needed so that we can put dark proteins in all proteins' scale.
    if (dark.target.only) {
        which <- key.df$DevLevel %in% "Tdark"
        key.df <- key.df[which, ]
    }
    # Pick the top 10 pathways only
    order <- order(key.df$Sum_Score, decreasing = TRUE)
    key.df <- key.df[order, ]
    # Top 20 pathways
    if (dim(key.df)[1] > 20)
        key.df <- key.df[1:20, ]
    # Reshape the dataframe
    key.df.melt <- melt(key.df, measure.vars = c("Average_Activation", "Average_Inhibition", "FDR"),
                         variable.name = "Method", value.name = "Relative_Score")
    plot <- ggplot(key.df.melt,
                   aes(x = reorder(key.df.melt[, x.col.name], -Sum_Score), y = Relative_Score, fill = Method)) +  # Need to order it
        geom_bar(stat = "identity", position = position_dodge()) + 
        xlab(x.lab) +
        ylab("Relative Score") +
        ggtitle(key) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    print(plot)
}

check.diff.annot.genes <- function(score.file.name) {
    df <- read.delim(score.file.name, header = TRUE, sep = "\t", check.names = FALSE)
    # Use negative log for FDR
    df$FDR <- -log10(df$FDR)
    # Scale it between 0 and 1
    df$FDR <- (df$FDR - min(df$FDR)) / (max(df$FDR) - min(df$FDR))
    col.names <- c("Average_Activation", "Average_Inhibition", "FDR")
    for (col.name in col.names) {
        t.test.result <- t.test(df[, col.name] ~ IsAnnotated, data = df)
        print(paste("Check ", col.name, sep = ""))
        print(t.test.result)
    }
    df.melt <- melt(df,
                    measure.vars = col.names,
                    variable.name = "Method", 
                    value.name = "Relative_Score")
    plot <- ggplot(df.melt, aes(x = Method, y = Relative_Score, fill = IsAnnotated)) +
            ylab("Relative Score") +
            geom_violin(position = position_dodge(0.9)) + 
            geom_boxplot(position = position_dodge(0.9), width = 0.02) +
            xlim(c("Average_Activation", "Average_Inhibition")) + ylim(c(0, 0.2)) + # Focus on these two methods
            stat_compare_means(label = "p.format", label.x = 1.5, label.y = 0.2)
    print(plot)
}

dir.name <- "/Users/wug/git/reactome-idg/fi-network-ml/results/impact_analysis"
# FDR for enrichment based analysis: FDR < 0.05
file.name <- "impact_analysis_092121_df_fdr_100721.txt"
score.name <- "FDR (-Log)"
need.dev.level <- FALSE
# Without filter
file.name <- "impact_analysis_092121_df_fdr_no_filter_100821.txt"
score.name <- "FDR (-Log)"
need.dev.level <- FALSE
# # Without filter: three protein familes for all dev levels
file.name <- "impact_analysis_092121_df_fdr_no_filter_three_families_all_100821.txt"
score.name <- "FDR (-Log)"
need.dev.level <- TRUE
# # Avearge activation from fuzzy logic simulation
file.name <- "impact_analysis_092121_df_activation_100821.txt"
score.name <- "Activation"
need.dev.level <- FALSE
# # Three protein families for all dev levels
file.name <- "impact_analysis_092121_df_activation_three_families_all_100821.txt"
score.name <- "Activation"
need.dev.level <- TRUE
# For inhibition
file.name <- "impact_analysis_092121_df_inhibition_100821.txt"
score.name <- "Inhibition"
need.dev.level <- FALSE
# # Three protein families for all dev levels
file.name <- "impact_analysis_092121_df_inhibition_three_families_all_100821.txt"
score.name <- "Inhibition"
need.dev.level = TRUE

# Some meta information
gene.family.file <- "/Users/wug/git/reactome-idg/fi-network-ml/src/main/resources/UniProtGeneDevLevelsTypes_100721.txt"
top.pathway.file <- paste(dir.name, "pathway2topic_100721.txt", sep = "/")
# When first run this, eter top.colors <- NULL, in the R consol first
if (exists("top.colors"))
    top.colors <- create.top.colors(top.pathway.file)

file.name <- paste(dir.name, file.name, sep = "/")
# plot.heatmap(file.name, score.name, gene.family.file, top.pathway.file, need.dev.level)

# Score file name
score.file.name <- paste(dir.name, "impact_analysis_092121_with_enrichment_092921_dev_level_100521.txt", sep = "/")
gene <- "C1QL1"
# gene <- "TANC1"
# plot.gene.pathways(score.file.name, gene)

# pathway <- "Neurexins and neuroligins"
# dark.target.only <- TRUE
# plot.pathway.genes(score.file.name, pathway, dark.target.only)

check.diff.annot.genes(score.file.name)
