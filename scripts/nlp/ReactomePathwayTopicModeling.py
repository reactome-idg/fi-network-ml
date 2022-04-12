# This script is used to conduct some topic modeling using text annotated for Reactome events using BERT.
import math
import os
import pickle
import random
import time
from typing import Tuple, List, Union, Dict

import bertopic
import numpy as np
import pandas as pd
import scanpy as sc
from sentence_transformers import SentenceTransformer
import logging

import UniProtHandler as uph
import PubmedHandler as ph
import TextEmbedder as embedder
import TopicModelingResultAnalyzer as resultAnalyzer

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# DIR = '/ssd/d0/ml/nlp/'
# At mac pro
DIR = "../../results/impact_analysis/nlp_files/"
TOP_PMID_NUMBER = 1000 # This is rather arbitray. Most for performance.
# DIR = '/Volumes/ssd/results/reactome-idg/fi-network-ml/impact_analysis/nlp_files/'
# PATHWAY_TEXT_FILE = DIR + "PathwayText_111921.txt"
PATHWAY_TEXT_FILE = DIR + "PathwayText_120721.txt"
# The lanugage model to be used
# BERT_LANGUAGE_MODEL = 'microsoft/BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext'
# For some reason, the aboe lanaguage model throws error. Use the folowing instead.
BERT_LANGUAGE_MODEL = 'dmis-lab/biobert-v1.1'
# PATHWAY_EMBEDDING_FILE = DIR + "pathway2embedding_biobert_112421.pkl"
# OUT_FILE = DIR + "pathway2embedding_112321.pkl"
# The default model: tried this based model. The performance of this model performs better than the above.
# Very weird!!!
# BERT_LANGUAGE_MODEL = 'bert-base-uncased'
# PATHWAY_EMBEDDING_FILE = DIR + 'pathway2embedding_bert_base_uncased_112421.pkl'
# BERT_LANGUAGE_MODEL = 'bert-large-cased'
# PATHWAY_EMBEDDING_FILE = DIR + 'pathway2embedding_bert_large_cased_112421.pkl'
# BERT_LANGUAGE_MODEL = 'roberta-base'
# PATHWAY_EMBEDDING_FILE = DIR + 'pathway2embedding_robert_base_112421.pkl'
UMAP_CLUSTER_FILE = DIR + "UMAP_Pathways_cluster_112421.html"
UMAP_TOPIC_FILE = DIR + "UMAP_Pathways_topic_112421.html"
PATHWAY_2_TOPIC_FILE = DIR + "../pathway2topic_100721.txt"
# This file was downloaded from pubmed after searching for TANC1 and then select abstract
TANC1_PUBMED_ABSTRACT_FILE = DIR + "abstract-TANC1-set.txt"
TANC1_ABSTRACT_COSINE_SIMILIARTY_FILE = DIR + "abstract_TANC1_Pathway_Cosine_Similarity_112421.tsv"
IMPACT_ANALYSIS_FILE = DIR + "../impact_analysis_092121_with_enrichment_092921.txt"
LAYERS = '-1'
# The batch size to be used at rws00061
# BATCH_SIZE = 50
# BATCH_SIZE = 3
# Minimum score to remove something that is not counted
MINIMUM_SCORE = 1.0e-22
SENTENCE_TRANSFORMER_MODEL = 'all-MiniLM-L6-v2'
MAX_SENTENCE_LENGTH = 512
PATHWAY_EMBEDDING_FILE = DIR + "pathway2embedding_" + SENTENCE_TRANSFORMER_MODEL + "_120721.pkl"


def load_pathway_document(file_name: str = PATHWAY_TEXT_FILE) -> Dict[str, List[str]]:
    """
    Load the summation text for individual pathway into a dict
    :param file_name:
    :return:
    """
    file = open(file_name, "r")
    pathway_name = None
    pathway_text_list = None
    text = None
    pathway2text = dict()
    for line in file:
        line = line.strip("\n")  # Get rid of new line first
        if line.startswith("###"):
            if pathway_name is not None:
                pathway2text[pathway_name] = pathway_text_list
            # Need to get rid of the new line char at the end for some reason
            pathway_name = line.strip("###")
            pathway_text_list = list()
            text = ""
        elif line == "//":
            pathway_text_list.append(text)
            text = ""
        else:
            text += line
    file.close()
    return pathway2text


def generate_topic_sentence_embedding(pathway2doc: Union[dict, str] = PATHWAY_TEXT_FILE,
                                      save: str = None) -> Dict[str, List[np.ndarray]]:
    if isinstance(pathway2doc, str):
        pathway2doc = load_pathway_document(PATHWAY_TEXT_FILE)
    pathway2embedding = embedder.generate_sentence_embedding(pathway2doc)
    if save is not None:
        file = open(save, 'wb')
        pickle.dump(pathway2embedding, file)
        file.close()
    return pathway2embedding


# To be run at the server
# generate_topic_sentence_embedding(save=PATHWAY_EMBEDDING_FILE)
def generate_topic_embedding(pathway2doc: Union[dict, str] = PATHWAY_TEXT_FILE,
                             language_model: str = BERT_LANGUAGE_MODEL,
                             save: str = PATHWAY_EMBEDDING_FILE) -> dict:
    """
    Generate the pathway to embedding using Flair
    :param pathway2doc
    :param language_model: the pre-trained language model to be used. See the documents
    for possible model that can be used in the Flair API.
    :param save: If anything is provided, the output will be saved
    :return:
    """
    if isinstance(pathway2doc, str):
        pathway2doc = load_pathway_document(PATHWAY_TEXT_FILE)
    pathway2embedding = embedder.generate_bert_embedding(pathway2doc)
    if save is not None:
        file = open(save, 'wb')
        pickle.dump(pathway2embedding, file)
        file.close()
    return pathway2embedding


def run_bert_topic() -> Tuple[bertopic.BERTopic, List[int], np.ndarray]:
    pathway2topic = load_pathway2topic()
    docs = list(pathway2topic.values())
    topic_model = bertopic.BERTopic()
    topics, probs = topic_model.fit_transform(docs)
    return topic_model, topics, probs


def run_topic_modeling(pathway2embedding: dict) -> Tuple[bertopic.BERTopic,
                                                         Union[np.ndarray, None]]:
    """
    Construct a BERTopic object for visualization.
    :param pathway2embedding:
    :return:
    """
    # Convert embedding from tenson to numpy as required by BERTopic
    # Pick a value
    embedding_matrix, pathway_list = embedder.convert_to_matrix(pathway2embedding)
    topic_model = bertopic.BERTopic()
    topics, probs = topic_model.fit_transform(pathway_list, embedding_matrix)
    return topic_model, probs


def run_umap(pathway2embedding: Union[dict, str] = PATHWAY_EMBEDDING_FILE,
             pathway2topic: dict = None) -> sc.AnnData:
    if pathway2topic is None:
        pathway2topic = load_pathway2topic(PATHWAY_2_TOPIC_FILE)
    if isinstance(pathway2embedding, str):
        # pathway2embedding = load_pathway2embedding(pathway2embedding)
        # Switch to use the sentencetransformer
        pathway2doc = load_pathway_document()
        pathway2embedding = generate_topic_sentence_embedding(pathway2doc)
    # Convert to a matrix
    embedding_matrix, pathway_list = embedder.convert_to_matrix(pathway2embedding)
    return resultAnalyzer.run_umap(embedding_matrix,
                                   pathway_list,
                                   pathway2topic,
                                   UMAP_CLUSTER_FILE,
                                   UMAP_TOPIC_FILE)


def load_pathway2embedding(file: str = PATHWAY_EMBEDDING_FILE) -> dict:
    logger.info("Loading pre-generated pathway2embedding...")
    file = open(file, "rb")
    pathway2embedding = pickle.load(file)
    file.close()
    logger.info("Done loading. The size of pathway2embeding: {}".format(len(pathway2embedding)))
    return pathway2embedding


def load_pathway2topic(file: str = PATHWAY_2_TOPIC_FILE) -> dict:
    file = open(file, 'r')
    pathway2topic = dict()
    for line in file:
        tokens = line.split("\t")
        pathway2topic[tokens[0]] = tokens[1]
    file.close()
    return pathway2topic


def load_impact_results(gene: str = None,
                        file: str = IMPACT_ANALYSIS_FILE) -> pd.DataFrame:
    df = pd.read_csv(file, sep='\t')
    if gene is not None:
        df = df.loc[df['Gene'] == gene, :]
    return df


def calculate_cor_impact_cosine_all() -> Tuple[float, float]:
    """
    Run all calculation in one single function
    :return:
    """
    print('loading pathway2embedding..')
    pathway2embedding = load_pathway2embedding()
    abstract_file_name = TANC1_PUBMED_ABSTRACT_FILE
    # abstract_file_name = DIR + 'abstract-LRFN1-set.txt'
    # abstract_file_name = DIR + 'abstract-NLGN1-set.txt'
    abstract = ph.load_pubmed_abstract(DIR, abstract_file_name)
    print('calculating embedding...')
    abstract_embedding = embedder.embed(abstract)
    print('calculating cosine similarity...')
    pathway2cosine = resultAnalyzer.calculate_cosine_similiarity(abstract_embedding,
                                                                 pathway2embedding,
                                                                 None)
    impact_results = load_impact_results('TANC1')
    print("Correlation for FDR:")
    results = resultAnalyzer.calculate_cor_impact_cosine(pathway2cosine,
                                                         impact_results,
                                                         'FDR',
                                                         True)
    print(results)
    print("Correation for Average_Activation:")
    results = resultAnalyzer.calculate_cor_impact_cosine(pathway2cosine,
                                                         impact_results,
                                                         'Average_Activation',
                                                         True)
    print(results)
    print("Correation for Average_Inhibition:")
    results = resultAnalyzer.calculate_cor_impact_cosine(pathway2cosine,
                                                         impact_results,
                                                         'Average_Inhibition',
                                                         True)
    print(results)


def calculate_cor_impact_cosine_via_sentence_transformer(gene: str,  # Three genes: TANC1, LRFN1, NLGN1
                                                         pathway2embedding: dict = None) -> Tuple[float, float]:
    """
    Run all calculation in one single function
    :return:
    """
    pathway2embedding, abstract_embedding, pathway2cosine, impact_results = _prepare_sent_cor_analysis(gene,
                                                                                                       pathway2embedding)
    print("Correlation for FDR:")
    results = resultAnalyzer.calculate_cor_impact_cosine(pathway2cosine,
                                                         impact_results,
                                                         'FDR',
                                                         True)
    print(results)
    print("Correation for Average_Activation:")
    results = resultAnalyzer.calculate_cor_impact_cosine(pathway2cosine,
                                                         impact_results,
                                                         'Average_Activation',
                                                         True)
    print(results)
    print("Correation for Average_Inhibition:")
    results = resultAnalyzer.calculate_cor_impact_cosine(pathway2cosine,
                                                         impact_results,
                                                         'Average_Inhibition',
                                                         True)
    print(results)
    return pathway2embedding


def search_abstracts_for_all_genes():
    """
    Perform batch search for all analyzed genes.
    :return:
    """
    logger.info("Load impact scores...")
    # Load impact result as a DataFrame
    impact_df = pd.read_csv(IMPACT_ANALYSIS_FILE, sep="\t")
    logger.info("The shape of impact_df: {}.".format(impact_df.shape))
    # Apparently there is a null gene included in the analysis. Need to exclude rows having for this!
    impact_df = impact_df.dropna()
    logger.info("Cleaned impact_df: {}.".format(impact_df.shape))
    genes = impact_df['Gene'].unique()
    logger.info("Total genes: {}.".format(len(genes)))
    ph.search_abstracts_for_all_via_ray(genes)


def batch_analyze_cor_impact_cosine():
    """
    Perform batch analysis using all downloaded pubmed abstracts.
    :return:
    """
    logger.info("Loading pathway embedding...")
    pathway2embedding = load_pathway2embedding()
    logger.info("Size of pathway2emebdding: {}.".format(len(pathway2embedding)))
    logger.info("Loading abstract emebdding...")
    pmid2emebdding = ph.load_pmid2embedding()
    logger.info("Size of abstract2embedding: {}.".format(len(pmid2emebdding)))
    ph.log_mem(logger)
    logger.info("Loading genes2pmids...")
    gene2pmids = ph.load_gene2pmids()
    logger.info("Size of gene2pmids: {}.".format(len(gene2pmids)))
    ph.log_mem(logger)
    logger.info("Load impact scores...")
    # Load impact result as a DataFrame
    impact_df = pd.read_csv(IMPACT_ANALYSIS_FILE, sep="\t")
    logger.info("The shape of impact_df: {}.".format(impact_df.shape))
    # Apparently there is a null gene included in the analysis. Need to exclude rows having for this!
    impact_df = impact_df.dropna()
    logger.info("Cleaned impact_df: {}.".format(impact_df.shape))
    genes = impact_df['Gene'].unique()
    logger.info("Total genes: {}.".format(len(genes)))
    # Do a filtering for easy handling
    which = np.isin(genes, list(gene2pmids.keys()))
    genes = genes[which]
    logger.info("Total genes in gene2pmids: {}.".format(len(genes)))
    # Used to select top genes
    pmid_reactome_sim_df = load_pmid2reactome_similarity_df()
    # Only need this sorted index
    pmid_reactome_sim_df_index = pmid_reactome_sim_df.index
    logger.info("Load pmid2reactome_sim_df.shape: {}.".format(pmid_reactome_sim_df.shape))
    # Keep all results in this DataFrame
    cols = ["Gene", "Pearson", "Peason_PValue", "Spearman", "Spearman_PValue", "Impacted_Pathways"]
    # Three types of impact scores
    impact_score_types = ["FDR", "Average_Activation", "Average_Inhibition"]
    impact_rows = [[], [], []]
    counter = 0
    time0 = time.time()
    # For local test
    genes = random.sample(genes.tolist(), 100)
    genes = ['DLG4', 'NLGN1', 'LRFN1', 'TANC1']
    logger.info("Genes subject to analysis: {}.".format(len(genes)))
    for gene in genes:
        logger.info("Handling gene: {}...".format(gene))
        logger.info("Searching pubmed abstracts")
        gene_pmids = gene2pmids[gene]
        if gene_pmids is None or len(gene_pmids) == 0:
            logger.info("Cannot find pmids for {}.".format(gene))
            continue
        logger.info("Found pmids: {}.".format(len(gene_pmids)))
        # Pick top pmids if there are too many pmids
        if len(gene_pmids) > TOP_PMID_NUMBER:
            # Don't use np.isin. Panda's index isin should be much faster
            # which = np.isin(pmid_reactome_sim_df.index, gene_pmids)
            which = pmid_reactome_sim_df_index.isin(gene_pmids)
            gene_pmids = pmid_reactome_sim_df_index[which].to_list()[:TOP_PMID_NUMBER]
            logger.info("Selected top {} PMIDS.".format(TOP_PMID_NUMBER))
        gene_pmid2embedding = {pmid: pmid2emebdding[pmid] for pmid in gene_pmids if pmid in pmid2emebdding.keys()}
        if len(gene_pmid2embedding) == 0:
            logger.info("No abstract emebedding for {}.".format(gene))
            continue
        gene_pathway2cosine = resultAnalyzer.calculate_cosine_similiarity(list(gene_pmid2embedding.values()),
                                                                          pathway2embedding)
        gene_impact_scores = impact_df.loc[impact_df['Gene'] == gene, ]
        for i in range(len(impact_score_types)):
            gene_cors = resultAnalyzer.calculate_cor_impact_cosine(gene_pathway2cosine,
                                                                   gene_impact_scores,
                                                                   impact_score_types[i])
            if gene_cors is None:
                logger.info("No enough impacted pathways for analysis!")
                continue;
            gene_row = {cols[0]: gene,
                        cols[1]: gene_cors[0][0],
                        cols[2]: gene_cors[0][1],
                        cols[3]: gene_cors[1].correlation,
                        cols[4]: gene_cors[1].pvalue,
                        cols[5]: gene_cors[2]}
            impact_rows[i].append(gene_row)
        counter += 1
        # if counter == 2:
        #     break;
    logger.info("Total genes analyzed: {}.".format(counter))
    time1 = time.time()
    logger.info("Total time used: {} seconds.".format(time1 - time0))
    ph.log_mem(logger)
    type2df = {impact_score_types[i]: pd.DataFrame(impact_rows[i], columns=cols) for i in range(len(impact_score_types))}
    # Save the files
    file_name = DIR + "{}_impact_pubmed_score_cor.txt"
    for type, df in type2df.items():
        type_file_name = file_name.format(type)
        df.to_csv(type_file_name, sep = '\t', index=False)
    return type2df


def _prepare_sent_cor_analysis(gene: str,
                               pathway2embedding: dict):
    if pathway2embedding is None:
        print('calculating pathway2embedding..')
        pathway2doc = load_pathway_document()
        pathway2embedding = generate_topic_sentence_embedding(pathway2doc)
    abstract = ph.load_pubmed_abstract(DIR, gene)
    print('calculating abstract embedding...')
    embedding_approach = SentenceTransformer(SENTENCE_TRANSFORMER_MODEL)
    embedding_approach.max_seq_length = MAX_SENTENCE_LENGTH
    abstract_embedding = embedding_approach.encode(abstract)
    print('calculating cosine similarity...')
    pathway2cosine = resultAnalyzer.calculate_cosine_similiarity(abstract_embedding,
                                                                 pathway2embedding,
                                                                 None)
    impact_results = load_impact_results(gene)
    return pathway2embedding, abstract_embedding, pathway2cosine, impact_results


def plot_sentence_transformer_cor_analysis(gene: str,
                                           pathway2embedding: dict = None) -> None:
    """
    Plot the results with plotly.
    :param gene:
    :param pathway2embedding:
    :return:
    """
    pathway2embedding, abstract_embedding, pathway2cosine, impact_results = _prepare_sent_cor_analysis(gene,
                                                                                                       pathway2embedding)
    return resultAnalyzer.plot_cor_analysis(gene, impact_results, pathway2cosine, DIR)


def plot_luna_nlp_cor_analysis(gene: str = 'UGT1A5') -> None:
    """
    This function is used to calculate and plot correlation for Augustin topic modeling results based on
    MeSH term usage.
    :param gene:
    :return:
    """
    dir_name = DIR + "../luna_nlp_results/"
    file_name = dir_name + "idg_{}tiab_and_free_full_textsb__topic_papers_pathway_scores.txt".format(gene.lower())
    print(file_name)
    # Have to manually parse the lines because of the format issue in the original file
    file = open(file_name, 'r')
    pathway2cosine = {}
    n = -1
    for line in file:
        n += 1
        if n == 0:
            continue
        tokens = line.split(',')
        pathway = tokens[0]
        # Get the DB_ID only and we need integer
        pathway = int(pathway.split('-')[2])
        pathway2cosine[pathway] = float(tokens[1])
    file.close()
    # Load the gene impact scores
    impact_results = load_impact_results(gene)
    # Make sure there is result
    if impact_results.empty:
        print("No impact results for gene {}.".format(gene))
        return
    resultAnalyzer.plot_cor_analysis(gene,
                                     impact_results,
                                     pathway2cosine,
                                     dir_name,
                                     use_pathway_name=False)
    # Perform correlation analysis here
    results = resultAnalyzer.calculate_cor_impact_cosine(pathway2cosine,
                                                         impact_results,
                                                         'FDR',
                                                         False)
    print(results)
    print("Correation for Average_Activation:")
    results = resultAnalyzer.calculate_cor_impact_cosine(pathway2cosine,
                                                         impact_results,
                                                         'Average_Activation',
                                                         False)
    print(results)
    print("Correation for Average_Inhibition:")
    results = resultAnalyzer.calculate_cor_impact_cosine(pathway2cosine,
                                                         impact_results,
                                                         'Average_Inhibition',
                                                         False)
    print(results)


def run_luna_nlp_cor_analysis(dir_name=DIR + "../luna_nlp_results"):
    # Get the list files for analysis
    genes = []
    for file in os.listdir(dir_name):
        if file.endswith('tiab_and_free_full_textsb__topic_papers_pathway_scores.txt'):
            # Get the gene name in the file
            gene = file.split('_')[1]
            index = gene.index('tiab')
            gene = gene[:index]
            genes.append(gene.upper())
    for gene in genes:
        print("Working on {}...".format(gene))
        plot_luna_nlp_cor_analysis(gene)
        print('Done\n')


def sort_pubmed_abstracts_on_similarity():
    """
    Sort the abstracts based on the similarity to text summartions in Reactome pathways.
    The comparison is based on text embedding. The final score is the average values between
    pathway embeddings and the abstract embedding.
    :return:
    """
    time0 = time.time()
    pathway2embedding = load_pathway2embedding()
    pmid2embedding = ph.load_pmid2embedding()
    time1 = time.time()
    logger.info("Time used to load pathway2embedding and pmid2embedding: {} seconds.".format(time1 - time0))
    ph.log_mem(logger)
    # For test
    # pmid2embedding = ph.sample_pmid2abstract(pmid2embedding, 100)
    logger.info("Total pmid2emedding for similarity: {}.".format(len(pmid2embedding)))
    pmid2similarity = resultAnalyzer.calculate_pathway_abstract_cosine_similarity_via_ray(pmid2embedding,
                                                                                          pathway2embedding)
    time2 = time.time()
    logger.info("Total pmid2similarity: {}.".format(len(pmid2similarity)))
    logger.info("Time for similarity: {} seconds.".format(time2 - time1))
    file = DIR + "pmid2reactome_similarity.pkl"
    file = open(file, 'wb')
    pickle.dump(pmid2similarity, file)
    file.close()

    ph.log_mem(logger)
    df = pd.DataFrame(pmid2similarity.items(), columns=('PMID', 'Similarity'))
    df.sort_values(by='Similarity', inplace=True, ascending=False)
    time3 = time.time()
    logger.info("Time for sorting in DataFrame: {}.".format(time3 - time2))
    ph.log_mem(logger)
    # Need to cache it into a file
    file = DIR + "pmid2reactome_similarity_df.pkl"
    file = open(file, "wb")
    pickle.dump(df, file)
    time4 = time.time()
    logger.info("Time for saving the dataframe: {}.".format(time4 - time3))
    return df


def load_pmid2reactome_similarity_df(file_name: str = DIR + "pmid2reactome_similarity_df.pkl") -> pd.DataFrame:
    file = open(file_name, 'rb')
    pd = pickle.load(file)
    # Want to use pmid as index for performance
    pd.set_index('PMID', inplace=True)
    return pd


if __name__ == '__main__':
    # search_abstracts_for_all_genes()
    results_dfs = batch_analyze_cor_impact_cosine()
    for impact_type, results_df in results_dfs.items():
        print("{}:\n{}".format(impact_type, results_df))
# calculate_cor_impact_cosine_via_sentence_transformer('LRFN1', load_pathway2embedding())
