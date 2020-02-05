# set workspace
setwd("~/Documents/wgm/work/reactome-idg/archs4")

# The following code is copied directly from archs4's help page:
# R script to download selected samples
# Copy code and run on a local machine to initiate download
# Check for dependencies and install if missing
packages <- c("rhdf5", "preprocessCore", "sva")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    print("Install required packages")
    source("https://bioconductor.org/biocLite.R")
    biocLite("rhdf5")
    biocLite("preprocessCore")
    biocLite("sva")
}
library("rhdf5")
library("preprocessCore")
library("sva")
destination_file = "human_matrix.h5"
extracted_expression_file = "example_expression_matrix.tsv"

# Check if gene expression file was already downloaded, if not in current directory download file form repository
if(!file.exists(destination_file)){
    print("Downloading compressed gene expression matrix.")
    url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"
    download.file(url, destination_file, quiet = FALSE)
} else{
    print("Local file already exists.")
}

# Selected samples to be extracted
samp = c("GSM1224927","GSM1066120","GSM1224923","GSM1224929","GSM1224924","GSM1066118","GSM1066119","GSM1224925","GSM1224930","GSM1872071","GSM2282084","GSM1872064","GSM1872067","GSM1704845")

# Retrieve information from compressed data
samples = h5read(destination_file, "meta/Sample_geo_accession")
# Identify columns to be extracted
sample_locations = which(samples %in% samp)

tissue = h5read(destination_file, "meta/Sample_source_name_ch1")
genes = h5read(destination_file, "meta/genes")
series = h5read(destination_file, "meta/Sample_series_id")
series = series[sample_locations]

print(paste("series", series, sep = ": "))

# extract gene expression from compressed data
expression = h5read(destination_file, "data/expression", index=list(1:length(genes), sample_locations))
# expression = h5read(destination_file, "data/expression", index=list(1:length(genes), 1:length(series)))
H5close()

# normalize samples and correct for differences in gene count distribution
expression = log2(expression+1)
expression = normalize.quantiles(expression)

rownames(expression) = genes
colnames(expression) = samples[sample_locations]

# correct batch effects in gene expression
batchid = match(series, unique(series))

print(paste("batchid", batchid, sep = ": "))

correctedExpression <- ComBat(dat=expression, batch=batchid, par.prior=TRUE, prior.plots=FALSE)
