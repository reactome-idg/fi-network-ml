# This script is used to process the topic analysis results.
import psutil
from sentence_transformers import util
import numpy as np
import pandas as pd
from typing import Tuple
from scipy import stats
import math
import scanpy as sc
import plotly.express as px
import ray
import logging
import time

MINIMUM_SCORE = 1.0e-22
RANDOM_STATE = 123456


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
    if len(cos_list) < 10: # This is arbitray for the time being
        return None
    # print("The size of list for cor: {}".format(len(cos_list)))
    return stats.pearsonr(cos_list, impact_list), stats.spearmanr(cos_list, impact_list), len(cos_list)


def run_umap(embedding_matrix,
             pathway_list,
             pathway2topic,
             umap_cluster_file,
             umap_topic_file) -> sc.AnnData:
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
    color_2_file = {'Cluster':umap_cluster_file,
                    'Topic':umap_topic_file}
    for color_key, save in color_2_file.items():
        df.sort_values(by=color_key, inplace=True)
        fig = px.scatter(df,
                         x='UMAP_1',
                         y='UMAP_2',
                         color=color_key,
                         hover_data=hove_data)
        fig.write_html(save)
    return adata


def plot_cor_analysis(gene,
                       impact_results,
                       pathway2cosine,
                       out_dir,
                       use_pathway_name = True):
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


def calculate_pathway_abstract_cosine_similarity_via_ray(pmid2emebedding,
                                                         pathway2embedding):
    """
    Use ray to parallel the calculation of cosine similiarty between pathways
    and abstracts.
    Note: This method is very similar to embed_abstracts_via_ray() in PubmedHandler.py.
    :param pmid2emebedding:
    :param pathway2embedding:
    :return:
    """
    num_workers = psutil.cpu_count(logical=False)
    ray.init(num_cpus=num_workers)
    logging.info("Initializing {} ray actors...".format(num_workers))
    workers = [CosineSimilarityCalculator.remote(pathway2embedding) for _ in range(num_workers)]
    pmids = list(pmid2emebedding.keys())
    start = 0
    # For the final run
    # Use a little bit buffer for the total jobs
    # step = 1000 * num_workers
    step = 200
    end = start + step
    pmid2similarity = {}
    while start < len(pmid2emebedding):
        pmids_sub = pmids[start:end]
        counter = 0
        time1 = time.time()
        logging.info("Starting calculation...")
        for pmid in pmids_sub:
            embedding = pmid2emebedding[pmid]
            workers[counter % num_workers].calculate_similarity.remote(pmid,
                                                                       embedding)
            counter += 1
        for worker in workers:
            pmid2similarity.update(ray.get(worker.get_pmid2similarity.remote()))
            worker.clean.remote()
        time2 = time.time()
        logging.info("Done calculation for {} - {} : {} seconds.".format(start, end, (time2 - time1)))
        start = end
        end += step
        if end > len(pmid2emebedding):
            end = len(pmid2emebedding)
    return pmid2similarity


def plot_cor_batch_results(file_name: str,
                           col_index: int = 1,
                           color_col_name: str = 'Tdark',
                           need_violin_plot: bool = False):
    """
    This methos is used to analyze the results from batch correlation analysis between impact scores
    and pubmed abstract similarity scores.
    :param file_name:
    :param out_file_name:
    :param col_index:
    :return:
    """
    cor_df = pd.read_csv(file_name, sep='\t')
    cor_df.set_index('Gene', inplace=True, drop=False)
    cor_df['-Log10(p-value)'] = cor_df.iloc[:, col_index + 1].map(lambda x : -math.log10(x))
    # For gene dev levels
    gene_dev_df_file = '../../src/main/resources/UniProtGeneDevLevelsTypes_100721.txt'
    gene_dev_df = pd.read_csv(gene_dev_df_file, sep='\t')
    # There are some dupicated gene names for the same UniProt accessions
    gene_dev_df.drop_duplicates('GeneName', inplace=True)
    gene_dev_df.set_index('GeneName', inplace=True)
    cor_df['DevLevel'] = gene_dev_df.loc[cor_df.index]['DevLevel']
    cor_df['Tdark'] = cor_df['DevLevel'].map(lambda x : x == 'Tdark')
    # Check if genes are annotated in Reactome
    reactome_gene_file_name = '../../results/impact_analysis/ReactomeGenes_Ver_77.txt'
    reactome_gene_file = open(reactome_gene_file_name, 'r')
    reactome_genes = [line.strip() for line in reactome_gene_file.readlines()]
    reactome_gene_file.close()
    cor_df['Annotated'] = cor_df['Gene'].map(lambda g : g in reactome_genes)
    # scatter between correlation and p-values
    fig = px.scatter(cor_df,
                     x=cor_df.columns[col_index],
                     y='-Log10(p-value)',
                     color=color_col_name,
                     hover_data=['Gene', cor_df.columns[col_index], cor_df.columns[col_index + 1]])
    fig.write_html(file_name + "." + color_col_name.lower() + ".scatter.html")
    if need_violin_plot:
        fig = px.violin(cor_df,
                        y = 'Pearson',
                        x='Annotated' if color_col_name == 'Tdark' else 'Tdark',
                        color=color_col_name,
                        box = True,
                        points='all',
                        hover_data=cor_df.columns)
        fig.write_html(file_name + "." + color_col_name.lower() + ".violin.html")
    return cor_df


if __name__ == '__main__':
    dir_name = '../../results/impact_analysis/nlp_files/'
    file_name = dir_name + 'Average_Inhibition_impact_pubmed_score_cor.txt'
    plot_cor_batch_results(file_name, need_violin_plot=True, color_col_name='Annotated')
    plot_cor_batch_results(file_name, need_violin_plot=True, color_col_name='Tdark')


@ray.remote
class CosineSimilarityCalculator(object):
    def __init__(self, pathway2embedding):
        self.pathway2embedding = pathway2embedding
        self.pmid2similarity = dict()
        # Logging cannot work here
        print("Initialized CosineSimilarityCalculator: {}.".format(self))

    def calculate_similarity(self, pmid, embedding):
        pathway2similarity = calculate_cosine_similiarity(embedding,
                                                          self.pathway2embedding)
        similarity_mean = np.mean(list(pathway2similarity.values()))
        self.pmid2similarity[pmid] = similarity_mean

    def get_pmid2similarity(self):
        return self.pmid2similarity

    def clean(self):
        self.pmid2similarity.clear()