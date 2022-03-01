# This script is used to handle pubmed related processes, e.g. download, parse, etc. The code is heavily based
# on this web post: https://www.toptal.com/python/beginners-guide-to-concurrency-and-parallelism-in-python
import concurrent.futures
import gzip
import logging
import os
import time
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from pathlib import Path
from urllib.request import urlopen
import xml.etree.ElementTree as et

import psutil

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

DOWNLOAD_URL = "https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/pubmed22n{:04d}.xml.gz"
OUT_DIR = "../../results/impact_analysis/pubmed_baseline"
MAX_ID = 1114 # The largest number of 2022 annual baseline of pubmed
MAX_WORKER = 4


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
    print('memory used before loading: {} M.'.format(mem))
    time0 = time.time()
    # Try to run in a multi-process way
    with ProcessPoolExecutor(max_workers=MAX_WORKER) as executor:
        files = [file for file in file.iterdir() if '.xml' in file.name]
        results = executor.map(extract_abstract, files)
    pmid2abstract = {}
    for result in results:
        pmid2abstract = {**result, **pmid2abstract}
    print("Total pmid2abstract: {}".format(len(pmid2abstract)))
    mem = psutil.Process().memory_info().rss / (1024 * 1024)
    print('memory used after loading: {} M.'.format(mem))
    time1 = time.time()
    print("Total time used: {} seconds.".format(time1 - time0))
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
        pmid2abstract[pmid] = abstractText.text
    return pmid2abstract


if __name__ == '__main__':
    extract_all_abstracts()
