"""
This is script is used to find a set of parameters to run random forest via grid search.
"""
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
import time
import plotly.express as px
from random import random

def run_grid_search(data_matrix_file: str = '/Users/conleyp/projects/fi-run/test_matrix.csv',
                    out_file_name: str = 'gs_iter_results.csv',
                    percent_of_sampling: float = 0.01,
                    parameters_for_search: dict = {
                                                    'n_estimators': [100],
                                                    'criterion': ['gini'],
                                                    'max_depth' : [2, 4, 6, 8, 10, 14, 18],
                                                    'min_samples_split': [2, 4, 6, 8],
                                                    'min_samples_leaf': [1, 2, 4, 6],
                                                    'max_features': ['sqrt', 'log2'],
                                                    'class_weight': ['balanced', 'balanced_subsample']
                                                    },
                    n_jobs: int = 12,
                    random_state = None):
    print('Loading data...')
    time0 = time.time()
    df = pd.read_csv(data_matrix_file)
    df.GenePair = df.GenePair.str.replace('\t', '_')
    df.set_index('GenePair', inplace=True)
    print('Finished loading with the size: {}.'.format(df.shape))
    time1 = time.time()
    print('Total time used: {}'.format(time1 - time0))

    scoring = ["accuracy", "balanced_accuracy", "precision", "recall", "f1", "roc_auc"]
    clf = RandomForestClassifier(random_state=random_state)
    size_pos = int(df[df.FI == 1].shape[0]*percent_of_sampling)
    print('Positive sampling size: {}'.format(size_pos))
    size_neg = int(df[df.FI == 0].shape[0]*percent_of_sampling)
    print('Negative samling size: {}'.format(size_neg))
    dfs = []

    time00 = time.time()
    for i in range(10):
        print('Working on iteration: {}'.format(i))
        time0 = time.time()
        df2 = pd.concat([
            df[df.FI == 1].sample(size_pos),
            df[df.FI == 0].sample(size_neg)
        ])
        X = df2[df2.columns[~df2.columns.isin(['FI'])]]
        y = df2['FI']
        grid_search = GridSearchCV(estimator=clf,
                                   param_grid=parameters_for_search,
                                   scoring=scoring,
                                   cv=10,
                                   n_jobs=n_jobs,
                                   verbose=3,
                                   refit=False)
        grid_search.fit(X, y)
        dfs.append(pd.DataFrame(grid_search.cv_results_))
        time1 = time.time()
        # This output may be mingled with the output from GridSearchCV because of multiple processes used.
        print('Done of iteration {}. Total time: {}'.format(i, (time1 - time0)))
    time01 = time.time();
    print('Total time for grid search: {}'.format(time01 - time00))
    pd.concat(dfs, keys=['iter'+str(i) for i in list(range(10))]).to_csv(out_file_name)


def plot_grid_search_results(result_file_name: str = '/Volumes/ssd/results/reactome-idg/fi-network-ml/grid_search/'
                                                     'gs_iter_results_01_sampling_121722.csv',
                             score_name: str = 'test_f1',
                             score_threshold: float = 0.14,
                             score_std_threshold: float = 0.05,
                             out_file_name: str = '/Volumes/ssd/results/reactome-idg/fi-network-ml/'
                                                  'grid_search/test_f1.html'):
    df = pd.read_csv(result_file_name)
    df.rename(columns={'Unnamed: 0':'iteration',
                       'Unnamed: 1':'run'},
              inplace=True)
    # Drop columns starting with split
    cols_to_be_selected = df.columns.map(lambda c : not c.startswith('split'))
    df = df[df.columns[cols_to_be_selected]]
    # df = pd.melt(df, id_vars=['run', 'params', 'iteration'], value_vars=['mean_test_f1', 'std_test_f1'])
    # Make sure only std_test_f1 < 0.10 is selected
    print('Before filtering: {}'.format(df.shape))
    score_col_name = 'mean_' + score_name
    score_std_col_name = 'std_' + score_name
    df = df.loc[(df[score_std_col_name] < (df[score_col_name] * score_std_threshold)) &
                (df[score_col_name] > score_threshold)]
    # Jitter the dots a little bit so that we can view them easiler
    df['run_jitter'] = df['run'].map(lambda v: v + random() * 0.3)
    df_sum = df.groupby('run')[score_col_name].agg(['count', 'min', 'max', 'mean', 'std'])\
        .sort_values(by='mean', ascending=False)
    df_sum['std_percent'] = df_sum.apply(lambda row : row['std'] / row['mean'] * 100, axis=1)
    print(df_sum.head(20))
    print('After filtering: {}'.format(df.shape))
    fig = px.scatter(df,
                     x='run_jitter',
                     labels=['Run', score_name],
                     y=score_col_name,
                     error_y=score_std_col_name,
                     color='iteration',
                     hover_data=['param_class_weight',
                                 'param_max_depth',
                                 'param_max_features',
                                 'param_min_samples_leaf',
                                 'param_min_samples_split',
                                 score_std_col_name])
    fig.write_html(out_file_name)
    return df_sum


# run_grid_search(percent_of_sampling=0.1,
#                 parameters_for_search={
#                     'n_estimators': [100],
#                     'criterion': ['gini'],
#                     'max_depth': [8, 10],
#                     'min_samples_split': [2, 4, 6, 8],
#                     'min_samples_leaf': [2, 4, 6],
#                     'max_features': ['sqrt'],
#                     'class_weight': ['balanced', 'balanced_subsample']
#                 })

# The final setting to be used
run_grid_search(percent_of_sampling=0.1,
                parameters_for_search={
                    'n_estimators': [100],
                    'criterion': ['gini'],
                    'max_depth': [8],
                    'min_samples_split': [4],
                    'min_samples_leaf': [2],
                    'max_features': ['sqrt'],
                    'class_weight': ['balanced_subsample']
                })

