# This script is used to conduct some topic modeling using text annotated for Reactome events using BERT.
import gc
import pickle
from typing import Tuple, List, Union

import bertopic
import numpy as np
import pandas
import scanpy as sc
import plotly.express as px
from flair.data import Sentence
from flair.embeddings import TransformerDocumentEmbeddings

# DIR = '/ssd/d0/ml/nlp/'
DIR = '/Volumes/ssd/results/reactome-idg/fi-network-ml/impact_analysis/nlp_files/'
PATHWAY_TEXT_FILE = DIR + "PathwayText_111921.txt"
# The lanugage model to be used
# BERT_LANGUAGE_MODEL = 'microsoft/BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext'
# For some reason, the aboe lanaguage model throws error. Use the folowing instead.
BERT_LANGUAGE_MODEL = 'dmis-lab/biobert-v1.1'
OUT_FILE = DIR + "pathway2embedding_112321.pkl"
UMAP_CLUSTER_FILE = DIR + "UMAP_Pathways_cluster_112321.html"
UMAP_TOPIC_FILE = DIR + "UMAP_Pathways_topic_112321.html"
PATHWAY_2_TOPIC_FILE = DIR + "../pathway2topic_100721.txt"
# The batch size to be used at rws00061
BATCH_SIZE = 50
# BATCH_SIZE = 3
RANDOM_STATE = 123456


def load_pathway_document(file_name:str = PATHWAY_TEXT_FILE) -> dict:
    """
    Load the summation text for individual pathway into a dict
    :param file_name:
    :return:
    """
    file = open(file_name, "r")
    pathway_name = None
    text = None
    pathway2text = dict()
    for line in file:
        if line.startswith("###"):
            if pathway_name is not None:
                pathway2text[pathway_name] = text
            # Need to get rid of the new line char at the end for some reason
            pathway_name = line.strip("###\n")
            text = ""
        else:
            text += line
    file.close()
    return pathway2text


def generate_topic_embedding(pathway2doc: dict,
                             language_model: str = BERT_LANGUAGE_MODEL,
                             save: str = None) -> dict:
    """
    Generate the pathway to embedding using Flair
    :param pathway2doc
    :param language_model: the pre-trained language model to be used. See the documents
    for possible model that can be used in the Flair API.
    :param save: If anything is provided, the output will be saved
    :return:
    """
    # The same emebeding may be used in multiple batches as long as the sentences are cleaned
    embedding_approach = TransformerDocumentEmbeddings(language_model)
    pathway2embedding = dict()
    # pathway_list = random.sample(list(pathway2doc.keys()), 10)
    pathway_list = list(pathway2doc.keys())
    doc_list = [pathway2doc[pathway] for pathway in pathway_list]
    # sent_list = [Sentence(doc) for doc in doc_list]
    # To control the memory usage, we need to send bak for process
    start = 0
    end = start
    while end < len(doc_list):
        end = start + BATCH_SIZE
        if end > len(doc_list):
            end = len(doc_list)
        print("Working on {} to {}".format(start, end))
        pathway_list_batch = pathway_list[start:end]
        doc_list_batch = doc_list[start:end]
        sent_list_batch = [Sentence(doc) for doc in doc_list_batch]
        embedding_approach.embed(sent_list_batch)
        start = end
        for i in range(len(sent_list_batch)):
            # Convert a tenson to a np.array for easy manupulate
            pathway2embedding[pathway_list_batch[i]] = sent_list_batch[i].embedding.detach().numpy()
        # Force gc: https://stackoverflow.com/questions/1316767/how-can-i-explicitly-free-memory-in-python
        del sent_list_batch
        gc.collect()
    print("The size of pathway2embeding: {}".format(len(pathway2embedding)))
    if save is not None:
        file = open(save, 'wb')
        pickle.dump(pathway2embedding, file)
        file.close()
    return pathway2embedding


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


def run_umap(pathway2embedding: dict,
             save: str = None,
             pathway2topic: dict = None,
             color_key: str = 'Cluster') -> sc.AnnData:
    # Convert to a matrix
    embedding_matrix, pathway_list = _convert_to_matrix(pathway2embedding)
    # Convert pathway2embedding as AnnData object
    adata = sc.AnnData(embedding_matrix, dict(obs_names=pathway_list))
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=5, metric='cosine')  # pp: preprocess
    sc.tl.leiden(adata, random_state=RANDOM_STATE)
    sc.tl.umap(adata, random_state=RANDOM_STATE)
    if save is not None:
        x_umap = adata.obsm['X_umap']
        df = pandas.DataFrame(x_umap, columns=['UMAP_1', 'UMAP_2'])
        df['Cluster'] = adata.obs['leiden'].to_list()
        # Want to split DB_IDs and Pathway Names into two list
        db_ids = [pathway.split('||')[0] for pathway in pathway_list]
        pathway_names = [pathway.split('||')[1] for pathway in pathway_list]
        df['Pathway'] = pathway_names
        df['DB_ID'] = db_ids
        hove_data = ['Pathway', 'Cluster', 'DB_ID']
        if pathway2topic is not None:
            topics = []
            for pathway in pathway_names:
                topics.append(pathway2topic[pathway.strip(' ')])
            df['Topic'] = topics
            hove_data.append('Topic')
        fig = px.scatter(df,
                         x='UMAP_1',
                         y='UMAP_2',
                         color=color_key,
                         hover_data=hove_data)
        fig.write_html(save)
    return adata


def load_pathway2embedding(file: str = OUT_FILE) -> dict:
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


# pathway2doc = load_pathway_document()
# pathwa2embedding = generate_topic_embedding(pathway2doc, save=OUT_FILE)