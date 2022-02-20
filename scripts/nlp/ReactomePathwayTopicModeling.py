# This script is used to conduct some topic modeling using text annotated for Reactome events using BERT.
import gc
import math
import os
import pickle
import re
from typing import Tuple, List, Union, Dict

import bertopic
import numpy as np
import pandas as pd
import scanpy as sc
import plotly.express as px
import scipy.spatial.distance
from scipy import stats
from flair.data import Sentence
from flair.embeddings import TransformerDocumentEmbeddings
from sentence_transformers import SentenceTransformer, util

# DIR = '/ssd/d0/ml/nlp/'
DIR = '/Volumes/ssd/results/reactome-idg/fi-network-ml/impact_analysis/nlp_files/'
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
RANDOM_STATE = 123456
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


def embed(document: str,
          language_model: str = BERT_LANGUAGE_MODEL) -> np.ndarray:
    embedding_approach = TransformerDocumentEmbeddings(language_model, layers=LAYERS, layer_mean=True)
    sentence = Sentence(document)
    embedding_approach.embed(sentence)
    return sentence.embedding.detach().numpy()


def generate_topic_sentence_embedding(pathway2doc: Union[dict, str] = PATHWAY_TEXT_FILE,
                                      save: str = None) -> Dict[str, List[np.ndarray]]:
    if isinstance(pathway2doc, str):
        pathway2doc = load_pathway_document(PATHWAY_TEXT_FILE)
    embedding_approach = SentenceTransformer(SENTENCE_TRANSFORMER_MODEL)
    embedding_approach.max_seq_length = MAX_SENTENCE_LENGTH
    pathway2embedding = dict()
    counter = 1
    for pathway, doc in pathway2doc.items():
        print("{}: {}".format(counter, pathway))
        embedding = embedding_approach.encode(doc)
        pathway2embedding[pathway] = embedding
        counter = counter + 1
        # if counter == 5:
        #     break
    print("The size of pathway2embeding: {}".format(len(pathway2embedding)))
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
    # The same emebeding may be used in multiple batches as long as the sentences are cleaned
    embedding_approach = TransformerDocumentEmbeddings(language_model,
                                                       layer_mean=True,
                                                       layers=LAYERS,
                                                       pooling='cls')  # cls produces the best results
    pathway2embedding = dict()
    # pathway_list = random.sample(list(pathway2doc.keys()), 10)
    # pathway_list = list(pathway2doc.keys())
    # doc_list = [pathway2doc[pathway] for pathway in pathway_list]
    # Want to run embedding sentence by sentence. Tried to use batch, which gives different results.
    # Need to do more reading about how this is implemented.
    counter = 1
    for pathway, doc in pathway2doc.items():
        print("{}: {}".format(counter, pathway))
        sentence = Sentence(doc)
        embedding_approach.embed(sentence)
        pathway2embedding[pathway] = sentence.embedding.detach().numpy()
        # Force gc: https://stackoverflow.com/questions/1316767/how-can-i-explicitly-free-memory-in-python
        del sentence
        gc.collect()
        counter = counter + 1
    print("The size of pathway2embeding: {}".format(len(pathway2embedding)))
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
    embedding_matrix, pathway_list = _convert_to_matrix(pathway2embedding)
    topic_model = bertopic.BERTopic()
    topics, probs = topic_model.fit_transform(pathway_list, embedding_matrix)
    return topic_model, probs


def _convert_to_matrix(pathway2embedding: dict) -> Tuple[np.ndarray, List[str]]:
    # Convert to a matrix
    embedding_matrix = None
    pathway_list = []
    for pathway in pathway2embedding.keys():
        pathway_list.append(pathway)
        embedding_array = pathway2embedding[pathway]
        if embedding_matrix is None:
            embedding_matrix = embedding_array
        else:
            embedding_matrix = np.vstack((embedding_matrix, embedding_array))
    return embedding_matrix, pathway_list


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
    embedding_matrix, pathway_list = _convert_to_matrix(pathway2embedding)
    # Convert pathway2embedding as AnnData object
    adata = sc.AnnData(embedding_matrix, dict(obs_names=pathway_list))
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=10, metric='cosine')  # pp: preprocess
    sc.tl.leiden(adata, random_state=RANDOM_STATE)
    sc.tl.umap(adata, random_state=RANDOM_STATE)
    # Plot
    x_umap = adata.obsm['X_umap']
    df = pd.DataFrame(x_umap, columns=['UMAP_1', 'UMAP_2'])
    # Don't use int for cluster. Otehrwise, a continuous color bar will be used
    df['Cluster'] = [i for i in adata.obs['leiden']]
    # Want to split DB_IDs and Pathway Names into two list
    db_ids = [pathway.split('||')[0] for pathway in pathway_list]
    pathway_names = [pathway.split('||')[1] for pathway in pathway_list]
    df['Pathway'] = pathway_names
    df['DB_ID'] = db_ids
    hove_data = ['Pathway', 'Cluster', 'DB_ID']
    topics = [pathway2topic[pathway.strip(' ')] for pathway in pathway_names]
    df['Topic'] = topics
    hove_data.append('Topic')
    # Do two plots
    color_2_file = {'Cluster':UMAP_CLUSTER_FILE,
                    'Topic':UMAP_TOPIC_FILE}
    for color_key, save in color_2_file.items():
        df.sort_values(by=color_key, inplace=True)
        fig = px.scatter(df,
                         x='UMAP_1',
                         y='UMAP_2',
                         color=color_key,
                         hover_data=hove_data)
        fig.write_html(save)
    return adata


def load_pathway2embedding(file: str = PATHWAY_EMBEDDING_FILE) -> dict:
    file = open(file, "rb")
    pathway2embedding = pickle.load(file)
    file.close()
    return pathway2embedding


def load_pathway2topic(file: str = PATHWAY_2_TOPIC_FILE) -> dict:
    file = open(file, 'r')
    pathway2topic = dict()
    for line in file:
        tokens = line.split("\t")
        pathway2topic[tokens[0]] = tokens[1]
    file.close()
    return pathway2topic


def load_pubmed_abstract(gene: str = 'TANC1') -> List[str]:
    """
    Load the downloaded abstract text from pubmed into a single string text.
    The implmentation may not be robust enough and should be done more test.
    :param file_name:
    :return:
    """
    file_name = DIR + 'abstract-{}-set.txt'.format(gene)
    text_list = []
    after_title = False
    is_in_title = False
    title = ''
    after_title = False
    is_in_abstract = False
    before_abstract = False
    abstract = ''
    file = open(file_name, 'r')
    for line in file:
        line = line.strip(' |\n')
        if re.match("^\\d*\\. \\w*", line):
            # Commit
            if len(title) > 0 and len(abstract) > 0:
                # See https://colab.research.google.com/drive/12hfBveGHRsxhPIUMmJYrll2lFU4fOX06
                # for use [SEP]
                text_list.append(title + "[SEP]" + abstract)
                title = ''
                abstract = ''
                after_title = False
            continue
        if len(line) == 0:
            if before_abstract:
                before_abstract = False
                is_in_abstract = True
            elif is_in_abstract:
                is_in_abstract = False
            elif is_in_title:
                is_in_title = False
                after_title = True
            elif not after_title:
                is_in_title = True
        elif line.startswith('Author information') or line.startswith('Erratum in'):
            before_abstract = True
            is_in_abstract = False
        else:
            if is_in_title and not after_title:
                title = title + " " + line
            elif is_in_abstract:
                abstract = abstract + " " + line
    # Pick up the last one
    if len(title) > 0 and len(abstract) > 0:
        text_list.append(title + "[SEP]" + abstract)
    file.close()
    return text_list


def calculate_cosine_similiarity(abstract_embedding: np.ndarray,
                                 pathway2embedding: dict,
                                 save: str = None) -> dict:
    pathway2cosine = dict()
    for pathway, embedding in pathway2embedding.items():
        # Note we want to have the similarity, not distance
        # cosine = 1.0 - scipy.spatial.distance.cosine(abstract_embedding, embedding)
        cosine = util.cos_sim(abstract_embedding, embedding)
        # Want to use the mean of this
        pathway2cosine[pathway] = cosine.mean().item() # Get the mean
    if save is not None:
        file = open(save, 'w')
        file.write('Pathway\tCosine_Similarity\n')
        for pathway, cosine in pathway2cosine.items():
            file.write('{}\t{}\n'.format(pathway,cosine))
        file.close()
    return pathway2cosine


def load_impact_results(gene: str = None,
                        file: str = IMPACT_ANALYSIS_FILE) -> pd.DataFrame:
    df = pd.read_csv(file, sep = '\t')
    if gene is not None:
        df = df.loc[df['Gene'] == gene, :]
    return df


def calculate_cor_impact_cosine(pathway2cosine: dict,
                                impact_results: pd.DataFrame,
                                col_name: str = 'FDR',
                                use_pathway_name: bool = True) -> Tuple[float, float]:
    cos_list = []
    impact_list = []
    pathway_col_name = 'PathwayName'
    if use_pathway_name is False:
        pathway_col_name = 'DBID'
    for pathway, cosine in pathway2cosine.items():
        if use_pathway_name is True:
            pathway = pathway.split('||')[1]
        # Check if we have impact score
        # This is so complicated in Python
        impact_score = impact_results[col_name][impact_results[pathway_col_name] == pathway]
        if impact_score.empty:
            continue
        impact_score = impact_score.iloc[0]
        if impact_score < MINIMUM_SCORE:
            continue  # Set a minimum
        if col_name == 'FDR':
            impact_score = -math.log10(impact_score)
        impact_list.append(impact_score)
        cos_list.append(cosine)
    print("The size of list for cor: {}".format(len(cos_list)))
    return stats.pearsonr(cos_list, impact_list), stats.spearmanr(cos_list, impact_list)


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
    abstract = load_pubmed_abstract(abstract_file_name)
    print('calculating embedding...')
    abstract_embedding = embed(abstract)
    print('calculating cosine similarity...')
    pathway2cosine = calculate_cosine_similiarity(abstract_embedding, pathway2embedding)
    impact_results = load_impact_results('TANC1')
    print("Correlation for FDR:")
    results = calculate_cor_impact_cosine(pathway2cosine, impact_results, 'FDR')
    print(results)
    print("Correation for Average_Activation:")
    results = calculate_cor_impact_cosine(pathway2cosine, impact_results, 'Average_Activation')
    print(results)
    print("Correation for Average_Inhibition:")
    results = calculate_cor_impact_cosine(pathway2cosine, impact_results, 'Average_Inhibition')
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
    results = calculate_cor_impact_cosine(pathway2cosine, impact_results, 'FDR')
    print(results)
    print("Correation for Average_Activation:")
    results = calculate_cor_impact_cosine(pathway2cosine, impact_results, 'Average_Activation')
    print(results)
    print("Correation for Average_Inhibition:")
    results = calculate_cor_impact_cosine(pathway2cosine, impact_results, 'Average_Inhibition')
    print(results)
    return pathway2embedding


def _prepare_sent_cor_analysis(gene: str,
                               pathway2embedding: dict):
    if pathway2embedding is None:
        print('calculating pathway2embedding..')
        pathway2doc = load_pathway_document()
        pathway2embedding = generate_topic_sentence_embedding(pathway2doc)
    abstract = load_pubmed_abstract(gene)
    print('calculating abstract embedding...')
    embedding_approach = SentenceTransformer(SENTENCE_TRANSFORMER_MODEL)
    embedding_approach.max_seq_length = MAX_SENTENCE_LENGTH
    abstract_embedding = embedding_approach.encode(abstract)
    print('calculating cosine similarity...')
    pathway2cosine = calculate_cosine_similiarity(abstract_embedding, pathway2embedding)
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
    return _plot_cor_analysis(gene, impact_results, pathway2cosine)


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
    _plot_cor_analysis(gene,
                       impact_results,
                       pathway2cosine,
                       use_pathway_name=False,
                       out_dir=dir_name)
    # Perform correlation analysis here
    results = calculate_cor_impact_cosine(pathway2cosine, impact_results, 'FDR', False)
    print(results)
    print("Correation for Average_Activation:")
    results = calculate_cor_impact_cosine(pathway2cosine, impact_results, 'Average_Activation', False)
    print(results)
    print("Correation for Average_Inhibition:")
    results = calculate_cor_impact_cosine(pathway2cosine, impact_results, 'Average_Inhibition', False)
    print(results)


def run_luna_nlp_cor_analysis(dir_name = DIR + "../luna_nlp_results"):
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


def _plot_cor_analysis(gene,
                       impact_results,
                       pathway2cosine,
                       use_pathway_name = True,
                       out_dir = DIR):
    # Do a little bit converting
    if use_pathway_name:
        pathway2cosine = {pathway.split('||')[1]: cosine for pathway, cosine in pathway2cosine.items()}
    # Add the similarity column into the impact_results for plot
    sim_cols = []
    pathway_key = 'PathwayName'
    if use_pathway_name is False:
        pathway_key = 'DBID'
    for pathway in impact_results[pathway_key]:
        if pathway in pathway2cosine.keys():
            similarity = pathway2cosine[pathway]
        else:
            similarity = None
        sim_cols.append(similarity)
    impact_results['Cosine_Similarity'] = sim_cols
    fdr_scores = [-math.log10(fdr) for fdr in impact_results['FDR']]
    # # Need to score it
    impact_results['FDR_Score'] = scale_scores(fdr_scores)
    impact_results['Activation_Score'] = scale_scores(impact_results['Average_Activation'].to_list())
    impact_results['Inhibition_Score'] = scale_scores(impact_results['Average_Inhibition'].to_list())
    fig = px.scatter(impact_results,
                     x='Cosine_Similarity',
                     y=['FDR_Score', 'Activation_Score', 'Inhibition_Score'],  # Support multiple columns
                     hover_data=['PathwayName', 'DBID', 'Gene'],
                     trendline='ols')
    out_file = out_dir + gene + "_Plot.html"
    fig.write_html(out_file)
    # Print out the top 5 pathways
    df = pd.DataFrame.from_dict(pathway2cosine, orient="index", columns=['Cosine_Similarity'])
    df.sort_values(by='Cosine_Similarity', ascending=False, inplace=True)
    print(df.head(5))
    fig = px.histogram(df, x='Cosine_Similarity', nbins=50)
    out_file = out_dir + gene + "_cosine_similiary_histgram.html"
    fig.write_html(out_file)
    return pathway2cosine, impact_results


def scale_scores(scores: list):
    min_score = min(scores)
    if min_score < MINIMUM_SCORE:
        min_score = MINIMUM_SCORE
    max_score = max(scores)
    scaled_score = []
    for score in scores:
        if score < MINIMUM_SCORE:
            scaled_score.append(None)
        else:
            scaled_score.append((score - min_score) / (max_score - min_score))
    return scaled_score
