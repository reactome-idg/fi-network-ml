# This script is used to handle pubmed related processes, e.g. download, parse, etc. The code is heavily based
# on this web post: https://www.toptal.com/python/beginners-guide-to-concurrency-and-parallelism-in-python
import concurrent.futures
import gzip
import logging
import os
import pickle
import random
import re
import time
import ray
import xml.etree.ElementTree as et
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from pathlib import Path
from typing import List
from urllib.request import urlopen

import psutil

import TextEmbedder as embedder

# A much more efficient parallel computing to avoid using buggy Python version
# See this for more information: https://towardsdatascience.com/10x-faster-parallel-python-without-python-multiprocessing-e5017c93cce1

# file_name = 'PubmedHandler_031722.log'
file_name = None
logging.basicConfig(level=logging.INFO,
                    filename=file_name,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

logger = logging.getLogger(__name__)

DOWNLOAD_URL = "https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/pubmed22n{:04d}.xml.gz"
OUT_DIR = "../../results/impact_analysis/pubmed_baseline"
MAX_ID = 1114 # The largest number of 2022 annual baseline of pubmed
MAX_WORKER = psutil.cpu_count(logical=False)

# Cache loaded pubmed abstracts, which should be big
_pmid2abstract = None

def ensure_out_dir(dir_name: str):
    path = Path(dir_name)
    if not path.exists():
        logger.info("{} doesn't exist. Create...".format(dir_name))
        path.mkdir()


def dowload_link(dir_name: str,
                 link: str):
    """
    Download the linked file into a local directory
    :param dir_name:
    :param link:
    :return:
    """
    logger.info("Download {}...".format(link))
    download_path = Path(dir_name + "/" + os.path.basename(link))
    with urlopen(link) as xml, download_path.open('wb') as f:
        f.write(xml.read())
    logger.info('Download {} Done.'.format(link))


def download():
    ensure_out_dir(OUT_DIR)
    # Get links for download
    links = [DOWNLOAD_URL.format(id) for id in range(1, MAX_ID + 1) if id % 200 == 0]
    # dowload_link(OUT_DIR, links[0])
    with ThreadPoolExecutor(max_workers=MAX_WORKER) as executor:
        # Need to use future in order to catch exceptions
        # See:https://www.digitalocean.com/community/tutorials/how-to-use-threadpoolexecutor-in-python-3
        futures = [executor.submit(dowload_link, OUT_DIR, link) for link in links]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                logger.error("Exception thrown:", exc_info=e)


def extract_all_abstracts(dir_name: str = OUT_DIR) -> dict:
    file = Path(dir_name)
    if not file.exists():
        raise FileNotFoundError("Cannot find directory: {}".format(dir_name))
    mem = psutil.Process().memory_info().rss / (1024 * 1024)
    logger.info('memory used before loading: {} M.'.format(mem))
    time0 = time.time()
    # Try to run in a multi-process way
    with ProcessPoolExecutor(max_workers=MAX_WORKER) as executor:
        files = [file for file in file.iterdir() if '.xml' in file.name]
        results = executor.map(extract_abstract, files)
    time01 = time.time()
    logger.info("Done parsing all files: {} seconds.".format(time01 - time0))
    logger.info("Starting merging dicts...")
    pmid2abstract = {}
    for result in results:
        pmid2abstract.update(result) # Don't use the unpacking way (i.e. **) for merging even though it is faster for
                                     # two dict merging. This should be faster since there is no need to unpack the merged
                                     # dict repeatedly.
    logger.info("Total pmid2abstract: {}".format(len(pmid2abstract)))
    mem = psutil.Process().memory_info().rss / (1024 * 1024)
    logger.info('memory used after loading: {} M.'.format(mem))
    time1 = time.time()
    logger.info("Total time used for merging: {} seconds.".format(time1 - time01))
    # Save the object
    logger.info("Dumping pmid2abstract...")
    file = open(dir_name + "/pmid2abstract.pkl", 'wb')
    pickle.dump(pmid2abstract, file)
    file.close()
    time2 = time.time()
    logger.info("Total time used: {} seconds.".format(time2 - time0))
    return pmid2abstract


def load_pmid2abstract(file_name: str = OUT_DIR + "/pmid2abstract.pkl") -> dict:
    logger.info("Loading saved pmid2abstract...")
    file = open(file_name, "rb")
    pmid2abstract = pickle.load(file)
    logger.info("Total abstracts loaded: {}.".format(len(pmid2abstract)))
    return pmid2abstract


def extract_abstract(file_name: str,
                     dir_name: str = OUT_DIR) -> dict:
    """
    Extract the XML encoded abstract from a zipped file.
    :param file_name:
    :return:
    """
    logger.info("Parsing {}...".format(file_name))
    if isinstance(file_name, Path):
        file = file_name
    else:
        file = Path(dir_name + "/" + file_name)
        if not file.exists():
            raise FileNotFoundError('{} cannot be found.'.format(file_name))
    # Need to open the file as unzip
    if file.name.endswith('.gz'):
        file = gzip.open(file_name, 'r')
    # Select all articles having abstract. But it is difficult to check AbstractText directly.
    pmid2abstract = {}
    root = et.parse(file).getroot()
    for citation in root.iterfind('./PubmedArticle/MedlineCitation'):
        # Escape if no abstract text
        abstractText = citation.find('./Article/Abstract/AbstractText')
        if abstractText is None:
            continue
        pmid = citation.find('PMID').text
        # Use the following method to get rid of tags in the abstract text, which are used for html rendering
        # This should work for text only abstractText.
        pmid2abstract[pmid] = ' '.join(abstractText.itertext())
    logger.info("Done parsing: {}.".format(file_name))
    return pmid2abstract


def search_abstract(gene: str) -> dict:
    """
    Search for abstracts for a gene. This is a simple text match for words and should be improved in the future.
    :param gene:
    :param pmid2abstract:
    :return:
    """
    global _pmid2abstract
    if _pmid2abstract is None:
        _pmid2abstract = load_pmid2abstract()
    found = {}
    gene = gene.lower()
    for pmid in _pmid2abstract.keys():
        abstract = _pmid2abstract[pmid]
        if gene in abstract.lower():
            found[pmid] = abstract
    return found


def load_pubmed_abstract(dir_name: str,
                         gene: str = 'TANC1') -> List[str]:
    """
    Load the downloaded abstract text from pubmed into a single string text.
    The implmentation may not be robust enough and should be done more test.
    :param file_name:
    :return:
    """
    file_name = dir_name + 'abstract-{}-set.txt'.format(gene)
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


def sample_pmid2abstract(pmid2abstract: dict,
                         k: int = 1000) -> dict:
    # For local test
    pmids = random.choices(list(pmid2abstract), k = k)
    pmid2abstract = {pmid: pmid2abstract[pmid] for pmid in pmids}
    logger.info("The size of pmid2abstract: {}".format(len(pmid2abstract)))
    return pmid2abstract


_sentence_transformer = None


def embed_abstract(abstract):
    # logger.info("{}: {}".format(pmid, abstract[0:30]))
    global _sentence_transformer
    if _sentence_transformer is None:
        _sentence_transformer = embedder.create_sentence_transformer()
    return _sentence_transformer.encode(abstract)


def embed_abstracts():
    """
    Embedding all pubmed abstracts
    :return:
    """
    pmid2abstract = load_pmid2abstract()
    pmid2abstract = sample_pmid2abstract(pmid2abstract, 1000)
    mem = psutil.Process().memory_info().rss / (1024 * 1024)
    logger.info("Memory used after loading abstracts but before embedding: {} MB.".format(mem))
    time0 = time.time()
    # Try to use multiple processes
    pmids = list(pmid2abstract.keys())
    abstracts = list(pmid2abstract.values())
    pmid2embedding = {}
    cache_file_name = OUT_DIR + "/pmid2embedding_cache.pkl"
    with ProcessPoolExecutor(max_workers=MAX_WORKER) as executor:
        for pmid, embedding in zip(pmids, executor.map(embed_abstract, abstracts)):
            pmid2embedding[pmid] = embedding
            if (len(pmid2embedding)) % 100 == 0:
                logger.info("Done: {}: ".format(len(pmid2embedding)))
                logger.info("Caching...")
                cache_obj(pmid2embedding, cache_file_name)
    time1 = time.time()
    logger.info("Done embedding: {} seconds.".format(time1 - time0))
    logger.info("Total embedding: {}.".format(len(pmid2embedding)))
    mem = psutil.Process().memory_info().rss / (1024 * 1024)
    logger.info("Memory used: {} MB.".format(mem))
    file = open(OUT_DIR + "/pmid2embedding.pkl", "wb")
    pickle.dump(pmid2embedding, file)


def embed_abstracts_via_ray():
    """
    Embedding all pubmed abstracts
    :return:
    """
    pmid2abstract = load_pmid2abstract()
    pmid2abstract = sample_pmid2abstract(pmid2abstract, 1000)
    mem = psutil.Process().memory_info().rss / (1024 * 1024)
    logger.info("Memory used after loading abstracts but before embedding: {} MB.".format(mem))
    time0 = time.time()
    # Try to use multiple processes via ray
    ray.init(num_cpus=MAX_WORKER)
    logger.info("Initializing {} ray actors...".format(MAX_WORKER))
    embedding_actors = [AbstractEmbedder.remote() for _ in range(MAX_WORKER)]
    counter = 0
    for pmid, abstract in pmid2abstract.items():
        embedding_actors[counter % MAX_WORKER].embed.remote(pmid, abstract)
        counter += 1
    time1 = time.time()
    logger.info("Done embedding: {} seconds.".format(time1 - time0))
    logger.info("Starting merging...")
    pmid2embedding = {}
    for embedding_actor in embedding_actors:
        pmid2embedding.update(ray.get(embedding_actor.get_pmid2embedding.remote()))
    time2 = time.time()
    logger.info("Done merging: {} seconds.".format(time2 - time1))
    logger.info("Total embedding: {}.".format(len(pmid2embedding)))
    mem = psutil.Process().memory_info().rss / (1024 * 1024)
    logger.info("Memory used: {} MB.".format(mem))
    file = open(OUT_DIR + "/pmid2embedding.pkl", "wb")
    pickle.dump(pmid2embedding, file)


def cache_obj(obj, file):
    # First check if the back-up file is there
    path_bak = Path(file + ".bak")
    if path_bak.exists():
        path_bak.unlink() # Delete this file
    path = Path(file)
    if path.exists():
        path.rename(file + ".bak")
    file = open(file, "wb")
    pickle.dump(obj, file)


# For ray worker
@ray.remote
class AbstractEmbedder(object):
    def __init__(self):
        self.sentence_transfomer = embedder.create_sentence_transformer()
        self.pmid2embedding = dict()
        logging.info("Initialized AbstractEmbedder: {}".format(self))

    def embed(self, pmid, abstract):
        self.pmid2embedding[pmid] = embedder.sentence_embed(abstract,
                                                            self.sentence_transfomer)
        if len(self.pmid2embedding) % 100 == 0:
            logging.info("{}: {}.".format(self, len(self.pmid2embedding)))

    def get_pmid2embedding(self):
        self.sentence_transfomer = None
        return self.pmid2embedding


if __name__ == '__main__':
    embed_abstracts_via_ray()
