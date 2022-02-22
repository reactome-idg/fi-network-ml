# This script is used to handle pubmed related processes, e.g. download, parse, etc. The code is heavily based
# on this web post: https://www.toptal.com/python/beginners-guide-to-concurrency-and-parallelism-in-python

import logging
import os
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from pathlib import Path
from urllib.request import urlopen, Request

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
    print(download_path)
    with urlopen(link) as xml, download_path.open('wb') as f:
        f.write(xml.read())
    logger.info('Download {} Done.'.format(link))


def download():
    ensure_out_dir(OUT_DIR)
    # Get links for download
    links = [DOWNLOAD_URL.format(id) for id in range(1, MAX_ID + 1) if id % 200 == 0]
    # dowload_link(OUT_DIR, links[0])
    with ThreadPoolExecutor(max_workers=MAX_WORKER) as executor:
        try:
            func = partial(dowload_link, OUT_DIR)
            executor.map(func, links)
        except Exception as e:
            logger.error("Exception thrown: " + e)


if __name__ == '__main__':
    download()