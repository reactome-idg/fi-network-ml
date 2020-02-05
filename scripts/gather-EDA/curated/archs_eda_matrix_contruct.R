usePackage <- function(p)
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("limma")
# BiocManager::install("rhdf5")
# BiocManager::install("preprocessCore")
# 
# usePackage("tidyverse")
# usePackage("doParalle")
# usePackage("tools")
# usePackage("rentrez")
# usePackage("foreach")
# usePackage("reshape")
# usePackage("dply")
# usePackage("plotly")
# usePackage("plyr")

require("rhdf5")
require("tools")
require("preprocessCore")
require("rentrez")
require("GEOquery")
require("doParallel")
require("limma")
require("reshape")
require("tidyverse")
require("plotly")
require("plyr")
require("foreach")
# -----------------------------------------------
options(stringsAsFactors = FALSE)
set.seed(1234)
# registerDoParallel(detectCores() - 4)
# -----------------------------------------------
# or use this higher spec 
# -----------------------------------------------
# n <- 10
# no_cores <- detectCores() - n
# cl <- makeCluster(no_cores) 
# registerDoParallel(cl)

# -----------------------------------------------
source("tissue_gsms.R") 
source("celline_gsms.R") 

names(celllines.metadata)[1] <- "tissue"
# tissues.metadata <- celllines.metadata # cell line only usecase 

# for aggregating both cellline and tissue metadata 
both.dat <- as.data.frame(rbind(tissues.metadata, celllines.metadata))

# remove duplicates between cell line and tissue metadata ------
gsm.duplicated <- unique(both.dat$gsm[duplicated(both.dat$gsm)])

both.dat$duplicated <- ifelse(both.dat$gsm %in% gsm.duplicated, 1, 0)
both.dat1 <- both.dat[which(both.dat$duplicated == 1), ]
both.dat0 <- both.dat[which(both.dat$duplicated == 0), ]

# keep annotations with cell line tags (may be biologically simpler)
both.dat1 <- both.dat1 %>% 
  group_by(gsm) %>% 
  filter(grepl("cl",tissue)) %>% 
  as.data.frame()

both.dat1 <- distinct(both.dat1, gsm, .keep_all= TRUE)
both.dat<- as.data.frame(rbind(both.dat0, both.dat1))
tissues.metadata <- both.dat[order(both.dat$tissue), ]

sub.sample <- T

if (sub.sample) {
  samp <- sample(tissues.metadata$gsm, 1000)
  tissues.metadata <- tissues.metadata[which(tissues.metadata$gsm %in% samp), ]
}

# -----------------------------------------------
# absolute paths to input and output files
input.path <- "input/"
output.path <- "output/"

# -----------------------------------------------
source_file <- paste0(input.path, 'human_matrix.h5')

if(!file.exists(source_file)){
  print(paste("Downloading compressed gene expression matrix to ", source_file))
  url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.rda"
  download.file(url, source_file, quiet = FALSE)
} else {
  print("Local file already exists.")
}

h5ls(source_file)

samples <- h5read(source_file, "meta/Sample_geo_accession")
tissue <- h5read(source_file, "meta/Sample_source_name_ch1")
genes <- h5read(source_file, "meta/genes")
platform_id <- h5read(source_file, "meta/Sample_platform_id")
series <- h5read(source_file, "meta/Sample_series_id")

# which samples do we have curated tissue/cellline metadta 
msk <- which(samples %in% tissues.metadata$gsm)

samples <- samples[msk]
tissue <- tissue[msk]
platform_id <- platform_id[msk]
series <- series[msk]

tissues.metadata <- tissues.metadata[which(tissues.metadata$gsm %in% samples), ]
tissues.metadata <- tissues.metadata[match(samples, tissues.metadata$gsm), ]
tissues.metadata$tissue_info <- tissue
tissues.metadata$platform_id <- platform_id
tissues.metadata$series <- series

# ---------------------------------------------------------------------
# filter samples by tissue or cell line of interest (in focused studies or local machine)
# NOTE: in large scale studies, only filter by samples with curated metadata on tissue/cellline
# ---------------------------------------------------------------------
# Unit of measure of these expressions are gene counts 
expression <- h5read(source_file, "data/expression", index=list(1:length(genes), msk))
H5close()

# --------------------------------------------------------------------
# normalize samples and correct for differences in 'gene count distribution'
# --------------------------------------------------------------------
expression <- log2(expression + 1)
expression <- normalize.quantiles(expression)

rownames(expression) <- genes
colnames(expression) <- samples

# batch effects in gene expression
batchid <- match(series, unique(series))
# --------------------------------------------------------------------
# gene counts present in at leaset 20% of the probes-data
# --------------------------------------------------------------------
gene.counts <- rownames(expression)[apply(expression, 1, function(x)  length(which(x > 0))) > ((20 * dim(expression)[2]) / 100)]
expression.norm <- expression[which(rownames(expression) %in% gene.counts),]

# ---------------------------------------------------------------------
# regress out batch effect and compute PCA 
# ---------------------------------------------------------------------
adjusted.expression <- limma::removeBatchEffect(expression.norm, batch=tissues.metadata$series, batch2=tissues.metadata$platform_id, 
                                                 design = model.matrix(~ tissue, data = tissues.metadata))

# save adjusted expression as csv 
ae.dat.df <- as.data.frame(adjusted.expression)
write.csv(ae.dat.df, paste0(output.path, "adjusted_tissue_expression.csv"), row.names = T)

# preserved genes/features between tissues should be above blue region on smooth-scatter
exp.mad <- apply(adjusted.expression, 1, mad)
exp.mean <- apply(adjusted.expression, 1, mean)
# smoothScatter(x = exp.mean, y = exp.mad)

# evaluate PCs 
PC <- prcomp(adjusted.expression, scale.=T, center = T)
# Plot first 2 PCs
plotdata <- data.frame(SampleID=rownames(PC$rotation), 
                       PC1=PC$rotation[,1], 
                       PC2=PC$rotation[,2])

plotdata$tissue <- tissues.metadata$tissue

# ---------------------------------------------------------------------
# detect sample outlier via PC1 for downstream analysis 
# ---------------------------------------------------------------------
tag.pca.outlier <- function(my.vector, sampleqc){
  # PCA outlier tagging 
  pop.sd <- sd(my.vector) * sqrt((length(my.vector) - 1)/(length(my.vector)))
  pop.mean <- mean(my.vector)
  pop.median <- median(my.vector)
  pop.mad <- mad(my.vector)
  
  # zscore of PC1s ----------------------------------------------
  z.score <- lapply(seq_along(my.vector), function(i) {
    data.point <- my.vector[i]
    z <- (data.point - pop.mean) / pop.sd
    if (is.finite(z)) {
      ifelse((abs(z) < 3), "non-outlier", "outlier")
    } else {
      "non-outlier"
    }
  })
  vec <- as.vector(do.call(cbind, z.score))
  
  samplepca.outliers <- as.data.frame(cbind(sampleqc, vec))
  names(samplepca.outliers) <- c("SampleID", "pca_outlier")
  samplepca.outliers[samplepca.outliers=="outlier"] <- 1
  samplepca.outliers[samplepca.outliers=="non-outlier"] <- 0
  
  return(samplepca.outliers)
}
samplepca.outlier <- tag.pca.outlier(plotdata$PC1, sampleqc = rownames(PC$rotation))
write.csv(samplepca.outlier, paste0(output.path, "sample_outliers_tissue.csv"), row.names = F)

# ---------------------------------------------------------------------
# plot and save 
# ---------------------------------------------------------------------
# batchid <- as.factor(batchid)
# 
# g1 <- ggplot(plotdata, aes(x=PC1, y=PC2, color=batchid)) + 
#   geom_point(aes( 
#     text = paste("SampleID: ", rownames(PC$rotation), 
#                  "</br> outlier_tag: ", samplepca.outlier$pca_outlier))) + 
#   labs(title = "")
# p1 <- plotly::ggplotly(g1)
# g2 <- ggplot(plotdata, aes(x=PC1, y=PC2, color=tissue)) + 
#   geom_point(aes( 
#     text = paste("SampleID: ", rownames(PC$rotation), 
#                  "</br> outlier_tag: ", samplepca.outlier$pca_outlier))) + 
#   labs(title = "")
# p2 <- plotly::ggplotly(g2)
# 
# htmlwidgets::saveWidget(as_widget(p1), paste0(output.path, "batchid_expression_pca.html"))
# htmlwidgets::saveWidget(as_widget(p2), paste0(output.path,"tissue_expression_pca.html"))

# PLOT PER TISSUE --------------------------------------------------
tissue.ls <- unique(tissues.metadata$tissue)

eda.dat <- lapply(tissue.ls, function(x){
  dat.sliced <- tissues.metadata[which(tissues.metadata$tissue %in% x), ]
  dat <- adjusted.expression[ ,which(colnames(adjusted.expression) %in% dat.sliced$gsm)]
  # evaluate PCs
  PC <- prcomp(dat, scale.=T, center = T)
  # Plot first 2 PCs
  plotdata <- data.frame(SampleID=rownames(PC$rotation),
                         PC1=PC$rotation[,1],
                         PC2=PC$rotation[,2])
  plotdata$tissue <- x
  plotdata$batchid <- as.factor(dat.sliced$series)

  samplepca.outlier <- tag.pca.outlier(plotdata$PC1, sampleqc = rownames(PC$rotation))
  g1 <- ggplot(plotdata, aes(x=PC1, y=PC2, color=batchid)) +
    geom_point(aes(
      text = paste("SampleID: ", rownames(PC$rotation),
                   "</br> outlier_tag: ", samplepca.outlier$pca_outlier))) +
    labs(title = "ARCHS4 Curated Data PCA", subtitle = x)
  p1 <- plotly::ggplotly(g1)

  write.csv(samplepca.outlier, paste0(output.path, paste0(x, "-sample-outliers-tissue.csv")), row.names = F)
  htmlwidgets::saveWidget(as_widget(p1), paste0(output.path, paste0(x, "-batchid-expression-pca.html")), selfcontained = TRUE)
})

# -----------------------------------------------------------------------
# Conduct univariate analysis to see if a signal exists by tissue (pheno)
# tissue ~ gene_i - on a gene X sample should evaluate a long tail distribution
# -----------------------------------------------------------------------
tissues.metadata$tissue_factors <- as.numeric(as.factor(tissues.metadata$tissue))

pvals <- adply(adjusted.expression, 1, function(x){
  dat <- cbind(x, tissues.metadata$tissue_factors)
  dat <- as.data.frame(dat)
  colnames(dat) <- c("gene",  "tissue")
  dat$gene <- as.numeric(dat$gene)
  dat$tissue <- as.numeric(dat$tissue)
  f <- as.formula("~ gene")
  fit <- eBayes(lmFit(dat$tissue, model.matrix(f, data = dat)))
  pval <- fit$p.value[2]
  return(pval)
}, .parallel = T, .expand = F)
colnames(pvals) <- c("gene", "pval")

# save non-ggplot or plotly plots in pdf file 
pdf(file = paste0(output.path, "test.pdf"),  wi = 10, he = 8)
hist(pvals$pval, xlab='P-Values of tissue ~ gene_i', breaks=200, col = 'lightblue', main = "ARCHS tissue metadata ~ each gene's gene expression")
smoothScatter(x = exp.mean, y = exp.mad, xlab = "Mean", ylab = "Median Absolute Deviation (MAD)")
dev.off()
