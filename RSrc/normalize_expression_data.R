# Perform normalization (using normalize.quantiles) on the expression data, save to new HDF
# This variation also includes batch correction (using ComBat).
#
# Author: sshorser
###############################################################################

# You will need to install these packages: BiocManager, preprocessCore, rhdf5, sva.

library("sva")
library("rhdf5")
library("preprocessCore")

# args from CLI
args <- commandArgs(TRUE)

print(args)

if (length(args) < 2)
{
    print("Two arguments are required: 1) path to source file 2) path to output file")
    stopifnot(length(args) < 2)
}

source_file <- args[1]
# source_file <- '/tmp/human_matrix.h5'

output_file <- args[2]
#output_file <- '/tmp/normalized_human_matrix.h5'

# Get the input RDA file.
if(!file.exists(source_file))
{
    print(paste("Downloading compressed gene expression matrix to ", source_file))
    url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.rda"
    download.file(url, source_file, quiet = FALSE)
} else {
    print("Local file already exists.")
}
print(paste("Loading ", source_file))
# load the RDA file
load(source_file)
print("Loading is complete!")
class(meta)
class(expression)

samples <- meta["Sample_geo_accession"]
class(samples)
genes <- meta["genes"]
class(genes)

length(genes)

print("correcting bad value")
# The value at 7045,166150 is -9, for some reason. No other negative numbers, so we'll just set this to 0.
expression[7045,166150] <- 0

# prepare for batch correction (these are used later)
series <- meta$Sample_series_id
batchid <- match(series, unique(series))

# normalized <- normalize.quantiles(expression[,1:100])
print("applying log2 transformation")
expression <- log2(expression+1)
print("normalizing expression values")
normalized <- normalize.quantiles(expression)
# Now, do batch correction
print("performing batch correction")
batchCorrectedNormalized <- ComBat(dat=normalized, batch=batchid, par.prior=TRUE, prior.plots=FALSE)

# cleanup before attempting to create new file.
if(file.exists(output_file))
{
  print("Removing pre-existing output file.")
  file.remove(output_file)
}

h5createFile(output_file)
h5createGroup(output_file,"data")
h5createDataset(output_file, "data/normalized_expression", dim(batchCorrectedNormalized), storage.mode=storage.mode(batchCorrectedNormalized), fillValue=0.0, chunk=c(200,200), level=6)
h5write(batchCorrectedNormalized, file=output_file, "data/normalized_expression")
h5closeAll()
