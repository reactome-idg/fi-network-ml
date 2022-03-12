# This script is used to embed text using SentenceEmbedder.
from sentence_transformers import SentenceTransformer
from flair.embeddings import TransformerDocumentEmbeddings
from flair.data import Sentence
from typing import Dict, List, Tuple
import numpy as np
import logging
import gc

logger = logging.getLogger(__name__)

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
LAYERS = '-1'

# Minimum score to remove something that is not counted
MINIMUM_SCORE = 1.0e-22
SENTENCE_TRANSFORMER_MODEL = 'all-MiniLM-L6-v2'
MAX_SENTENCE_LENGTH = 512


def sentence_embed(text: str,
                   embedding_approach: SentenceTransformer = None):
    if embedding_approach is None:
        embedding_approach = create_sentence_transformer()
    return embedding_approach.encode(text)


def sentence_embed_pmid(abstract: str,
                        pmid: str,
                        embedding_approach: SentenceTransformer):
    embedding = sentence_embed(abstract)
    return {pmid: embedding}


def create_sentence_transformer() -> SentenceTransformer:
    logger.info("Create a new SentenceTransformer object...")
    embedding_approach = SentenceTransformer(SENTENCE_TRANSFORMER_MODEL)
    embedding_approach.max_seq_length = MAX_SENTENCE_LENGTH
    return embedding_approach


def generate_sentence_embedding(pathway2doc,
                                embedding_approach: SentenceTransformer = None) -> Dict[str, List[np.ndarray]]:
    logger.info("generate_sentence_embedding for {}...".format(len(pathway2doc)))
    if embedding_approach is None:
        embedding_approach = create_sentence_transformer()
    pathway2embedding = dict()
    counter = 1
    for pathway, doc in pathway2doc.items():
        # logger.info("{}: {}".format(counter, pathway))
        embedding = sentence_embed(doc, embedding_approach)
        pathway2embedding[pathway] = embedding
        counter = counter + 1
        # if counter == 5:
        #     break
    logger.info("The size of pathway2embeding: {}".format(len(pathway2embedding)))
    return pathway2embedding


def bert_embed(document: str,
               language_model: str = BERT_LANGUAGE_MODEL) -> np.ndarray:
    embedding_approach = TransformerDocumentEmbeddings(language_model, layers=LAYERS, layer_mean=True)
    sentence = Sentence(document)
    embedding_approach.embed(sentence)
    return sentence.embedding.detach().numpy()


def generate_bert_embedding(pathway2doc: dict,
                            language_model: str = BERT_LANGUAGE_MODEL) -> dict:
    """
    Generate the pathway to embedding using Flair
    :param pathway2doc
    :param language_model: the pre-trained language model to be used. See the documents
    for possible model that can be used in the Flair API.
    :param save: If anything is provided, the output will be saved
    :return:
    """
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
    logging.info("The size of pathway2embeding: {}".format(len(pathway2embedding)))
    return pathway2embedding


def convert_to_matrix(pathway2embedding: dict) -> Tuple[np.ndarray, List[str]]:
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

