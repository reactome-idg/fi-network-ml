# This script is used to process the topic analysis results.

from sentence_transformers import util
import numpy as np
import pandas as pd
from typing import Tuple
import statsmodels as stats
import math
import scanpy as sc
import plotly.express as px

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
    print("The size of list for cor: {}".format(len(cos_list)))
    return stats.pearsonr(cos_list, impact_list), stats.spearmanr(cos_list, impact_list)


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