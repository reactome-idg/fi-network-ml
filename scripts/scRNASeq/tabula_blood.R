BiocManager::install('iaconogi/bigSCale2')
library(bigSCale)

# Load filtered and converted tabula data
blood <- readRDS('TS_Blood_filtered.rds')
counts <- blood@assays$RNA@counts
genes <- rownames(blood)
rm(blood)

# Running bigSCale2 on tabula blood dataset
results_blood <- compute.network(expr.data=counts, gene.names=genes, speed.preset='fast')
corr <- results_blood$correlations
saveRDS(results_blood, 'results_blood.rds')
write.csv(corr, 'tabula_blood_corr_recursive.csv')

