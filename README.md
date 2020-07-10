# idg-fi-network-ml
This project is used to collect features from different sources and generate files subject to machine learning approaches to predict functional interactions among dark proteins and annotated proteins in Reactome. The final output actually is an expanded functional interaction network based on a much larger feature space than the one used to construct the Reactome FI network.

## Features used for the machine learning

### Protein-protein interactions
Protein-protein interactions are collected from StringDB, BioGrid, and BioPlex. For the actual files used, see the configuration file, application.properties, in folder src/main/resource. Protein-protein interactions in non-human species were mapped to human by using the orthologous mappings downloaded from panther for fly, worm, and yeast or from ensembl for mouse. The choice of the mapping sources was based on the numbers of mappable proteins from non-human species to human. For example, using the ensebml mapping for mouse produced more interaction pairs than using the panther mapping.

### Gene ontology annotation sharing
Following the procedure established for constructing the Reactome FI network, a  GO biological sharing was used as a feature. The human gene ontology annotation file was downloaded from the gene ontology web site. Gene pairs were checked to see whether or not the two genes in the pairs shared at least one GO BP term. 

### Domain-domain interaction
Following the procedure used for constructing the Reactome FI network, pFam domain-domain interactions are used to check if two proteins could interact via their domain annotation and domain-domain interactions.

### Harmonizome gene similarity
The gene similarity files were downloaded from the harmonizome web site: http://amp.pharm.mssm.edu/Harmonizome/download. For example, the following link was used for archilles: https://amp.pharm.mssm.edu/static/hdfs/harmonizome/data/achilles/gene_similarity_matrix_cosine.txt.gz. We manually curated all collected data sets in the harmonizome (Harmonizome_datasets_annotations_062819.txt in src/main/resources), chose data sets focusing on pair-wise relationships, excluding data sets collected from Reactome and other pathway data sources. Furthermore, we conducted an odds ratio analysis using the FIs extracted from Reactome pathways based on the Java code: [FeatureChecker](https://github.com/reactome/fi_network_build/blob/master/src/org/reactome/fi/util/FeatureChecker.java). We chose datasets having odds ratio >= 4.98, excluding two files having similarity cutoff = 0.0 (clinvar and omim). The final threshold values used for similarities were also based on numbers of pairs that could be categorized as positive. The actual code uses a manually tuned file, harmonizome_selected_files.txt, in the src/main/resources folder. In total, 20 features were selected from this data source.

### GTEx tissue and TCGA cancer-specific gene co-expression
The original gene expression samples were downloaded from the GTEx and TCGA GDC data portals. The [reactome-idg/gather-app](https://github.com/reactome-idg/gather-app) pre-processed, QA/QC'ed (PCA, MCL clustering, pathway enrichment analysis), and constructed the final spearman correlation adjacency of all normal and cancer specific tissues. To choose the threshold values for the calculated Spearman correlations, we chose an adaptive cutoff as described in paper by Iacono G et al: [Single-cell transcriptomics unveils gene regulatory network plasticity](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1713-4). As in the original paper, we chose top 0.1% correlations as positive features for individual tissue- or cancer-specific correlations. Furthermore, we also ensured that the odds ratio for each selected correlation feature should be >= 5.0. With these odds ratio threshold, two GTEx data, Brain-Putamen-basalganglia_Spearman_Adj (odd ratio: 4.58 +- 0.35) and Brain-AnteriorcingulatecortexBA24_Spearman_Adj (odd ratio: 3.57 +- 0.19), and one TCGA file, TCGA-UVM_Spearman_Adj (odds ratio = 4.22 +- 0.16), were excluded as features.

---

In total, we got 106 features: 5 PPIs + 1 domain interaction + 1 BP sharing + 20 harmonizime gene similarities + 31 TCGA coexpression + 48 GTEx coexpression.


