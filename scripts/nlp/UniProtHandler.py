# This script is used to pull gene names and protein names together from UniProt database so that we can do synonym
# based search.
import pandas as pd
import re
import logging

SOURCE_FILE_NAME = "../../results/impact_analysis/uniprot/uniprot_protein_gene_names_030422.tab"

# cache loaded abstracts, which should be big
_gene2names = None


def load_gene2othernames(file_name: str = SOURCE_FILE_NAME) -> dict:
    """
    Load the map from gene symbols to other names, including synonyms and protein names from the
    specified UniProt table file. The table is downloaded from UniProt via its table view.
    :param file_name:
    :return:
    """
    logging.info("Loading gene2othernames...")
    df = pd.read_csv(file_name, sep="\t")
    gene2names = {}
    for i in range(df.shape[0]):
        names = set()
        names.add(df.iloc[i, 3]) # Primary gene names
        # Get gene synonyms
        if not pd.isnull(df.iloc[i, 4]):
            syns = df.iloc[i, 4].split(" ")
            names.update(syns)
        # Get protein names
        if not pd.isnull(df.iloc[i, 2]):
            protein_names = _parse_protein_names_via_split(df.iloc[i, 2])
            names.update(protein_names)
        gene2names[df.iloc[i, 3]] = names
    return gene2names


def _parse_protein_names_via_split(text: str) -> list:
    # Remove any contained names
    match = re.search('\\[.*\\]', text)
    if match:
        text = text[0:match.start() - 1].strip(' ')
    names = list()
    tokens = text.split(" (")
    for token in tokens:
        # Get rid of ")" if any
        token = token.strip(' ')
        if token.endswith(')'):
            token = token[0:len(token) - 1]
        if len(token) == 0:
            continue; # Escape contained name or empty space
        _add_protein_name(names, token)
    return names


def _add_protein_name(names: list,
                      name: str) -> list:
    name = name.strip(' ')
    # Avoid EC number
    if (re.match('^EC (\d*|-)\\.(\d*|-)\\.(\d*|-)\\.(\d*|-)$', name)):
        return names
    names.append(name)
    return names


def get_names(gene: str) -> dict:
    """
    Get all names for a gene as collected in UniProt
    :param gene:
    :return:
    """
    # Want to get the defined variable outside of this function scope
    global _gene2names
    if _gene2names is None:
        _gene2names = load_gene2othernames()
    return _gene2names[gene]

