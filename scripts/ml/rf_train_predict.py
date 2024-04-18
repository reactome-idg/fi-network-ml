"""
This script is used to train a random forest using the training data set, measure its performance using the test 
data set. This script is also used to predict FIs using the matrix with gene/protein pairwise features.
Note: The matrices used for training, testing and prediction are generated by the Java code in this project. The parameters
used for the random forest model are found using grid-search by another Python script, grid_search_rf.py.
"""

import pandas as pd
import pickle as pkl
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import precision_recall_curve
from pathlib import Path
import numpy as np
from F1Analyzer import calculate_f1, analyze_rf_performance

import logging
logging_file_name = None
logging.basicConfig(level=logging.INFO,
                    filename=logging_file_name,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def train_rf(train_file: str,
             test_file: str,
             working_dir: str,
             repeats: int = 3,
             postfix_name: str = '',
             down_sampling: float = None):
    """This function is used to train a RF and then use the trained RFs to assign scores to the prediction file. The code was merged 
    from fit_predict.py, originally developed by Patrick Conley.

    Args:
        train_file (str): the data file used to train
        test_file (str): the data file used to test. No FIs in this file should be in the train_file.
        working_dir (str): the directory where outputs will be saved.
        repeats (int, optional): Since the selected parameters may result some stochastic results, it is recommended 
        to run the train multiple times. Defaults to 3.
        postfix_name (str, optional): a text that is used to postfix the training result files. Defaults to ''.
        down_sampling (float, optional): if downsampling is needed to checking. Default to None.
    """
    # Training data
    logger.info('Loading training file...')
    df = _load_data_file(train_file, down_sampling)
    X = df[df.columns[~df.columns.isin(['FI'])]]
    y = df['FI']
    logger.info('Done')

    # Test data
    logger.info('Loading test file...')
    df_test = _load_data_file(test_file, down_sampling=down_sampling)
    X_test = df_test[df_test.columns[~df_test.columns.isin(['FI'])]]
    logger.info('Done')

    # Keep the test results since we need to aggregate them together
    test_results = []

    logger.info('Training and testing...')
    for i in range(repeats):
        logger.info('Training {}...'.format(i))
        # The following parameters are selected by gridsearch using script grid_search_rf.py.
        rfc = RandomForestClassifier(
            n_estimators=100,
            class_weight='balanced_subsample',
            criterion='gini',
            max_depth=8,
            min_samples_split=4,
            min_samples_leaf=2,
            max_features='sqrt',
            n_jobs=16)

        rfc.fit(X, y)
        test_result_file = 'test_predict_scores_{}_{}.csv'.format(
            i, postfix_name)
        test_result = pd.DataFrame(rfc.predict_proba(X_test))
        test_result.index = df_test.index  # Use the original index
        test_result.to_csv(Path(working_dir, test_result_file))
        test_results.append(test_result)
        trained_rf_file = _get_trained_rfc_file(i, postfix_name)
        pkl.dump(rfc, open(Path(working_dir, trained_rf_file), "wb"))
        logger.info('Done')

    # Aggregate all test_results into one
    logger.info('Analyzing prediction results...')
    aggregated_test_result = aggreate_prediction_results(
        test_result_dfs=test_results)
    test_file_name = 'aggregated_test_predict_scores_{}.csv'.format(
        postfix_name)
    # Just want to make the same format
    aggregated_test_result.to_csv(
        Path(working_dir, test_file_name), index=None)
    logger.info('Done.')

    # Calculate precicison, recall and F1
    logger.info('Generate precision/recall data...')
    precision_recall_df, pr_file = generate_precision_recall(test_df=df_test,
                                                             aggregated_test_result=aggregated_test_result,
                                                             working_dir=working_dir,
                                                             postfix_name=postfix_name)
    logger.info('Done')

    logger.info('Plot the results...')
    # Analyze the results
    # In case pr_file is a Path
    if isinstance(pr_file, Path):
        pr_file = str(pr_file)
    html_file = pr_file.replace('.csv', '.html')
    analyze_rf_performance(precision_recall_df=precision_recall_df,
                           need_plot=False,
                           html_file=html_file)
    logger.info('Done')

    return aggregated_test_result, precision_recall_df


def _load_data_file(file_name: str,
                    down_sampling: float):
    df = pd.read_csv(file_name)
    if down_sampling is not None:
        df = df.sample(frac=down_sampling)
    df.GenePair = df.GenePair.str.replace('\t', '_')
    df.set_index('GenePair', inplace=True)
    # A bug in some old files
    if 'YeatPPI' in df.columns:
        df.rename(columns={'YeatPPI': 'YeastPPI'}, inplace=True)
    return df


def _get_trained_rfc_file(repeat: int,
                          postfix_name: str):
    trained_rf_file = 'trained_rfc_{}_{}.pkl'.format(repeat, postfix_name)
    return trained_rf_file


def aggreate_prediction_results(test_result_dfs: list) -> pd.DataFrame:
    # column names are 0 and 1. They are integer!!!
    score_dict = {'rfc_{}'.format(
        index): test_result_df[1] for index, test_result_df in enumerate(test_result_dfs)}

    aggregated_df = pd.DataFrame(score_dict)
    mean = aggregated_df.apply(np.mean, axis=1, raw=True)
    stdev = aggregated_df.apply(np.std, axis=1, raw=True)
    aggregated_df['mean'] = mean
    aggregated_df['stdev'] = stdev
    aggregated_df['%cv'] = (aggregated_df['stdev']/aggregated_df['mean'])*100
    aggregated_df.index = test_result_dfs[0].index

    return aggregated_df


def generate_precision_recall(test_df: pd.DataFrame,
                              aggregated_test_result: pd.DataFrame,
                              working_dir: str = '',
                              postfix_name: str = ''):
    precision, recall, thresholds = precision_recall_curve(test_df['FI'],
                                                           aggregated_test_result['mean'])  # Use mean as the score
    # Create a DataFrame for precision-recall pairs
    precision_recall_df = pd.DataFrame({
        'Threshold': thresholds,
        'Precision': precision[:-1],  # Last value is meaningless, removing it
        'Recall': recall[:-1]  # Last value is meaningless, removing it
    })
    calculate_f1(precision_recall_df)

    # Dump it into a file
    file_name = 'precision_recall_{}.csv'.format(postfix_name)
    file = Path(working_dir, file_name)
    precision_recall_df.to_csv(Path(working_dir, file_name))
    return precision_recall_df, file


def rf_predict(predict_file: str,
               working_dir: str,
               repeats: int = 3,
               postfix_name: str = '',
               down_sampling: float = None) -> pd.DataFrame:
    """Generate the scores using the trained RFC. There should be three.

    Args:
        predict_file (str): _description_
        working_dir (str): _description_
        postfix_name (str, optional): _description_. Defaults to ''.

    Returns:
        pd.DataFrame: _description_
    """
    # Check if the trained RFCs are there
    logger.info('Loading trained RFCs...')
    trained_rfcs = []
    for i in range(repeats):
        file_name = _get_trained_rfc_file(i, postfix_name)
        file = Path(working_dir, file_name)
        if not file.exists():
            raise ValueError('{} not in {}!'.format(file_name, working_dir))
        rfc = _load_trained_rfc(file)
        trained_rfcs.append(rfc)
    logger.info('Done')

    # Test data
    logger.info('Loading predictioin file...')
    predict_df = _load_data_file(predict_file, down_sampling)
    logger.info('Done')

    logger.info('Calculate score...')
    predict_result_dfs = []
    for rfc in trained_rfcs:
        predict_result = rfc.predict_proba(predict_df)
        predict_result_df = pd.DataFrame(predict_result)
        predict_result_df.index = predict_df.index
        predict_result_dfs.append(predict_result_df)
    logger.info('Done')

    logger.info('Aggregate results...')
    aggregated_results = aggreate_prediction_results(predict_result_dfs)
    logger.info('Done')

    file_name = 'prediction_results_{}.csv'.format(postfix_name)
    file = Path(working_dir, file_name)
    aggregated_results.to_csv(file)
    logger.info('Saved results into: {}'.format(file))

    return aggregated_results


def _load_trained_rfc(file_name: Path):
    with open(file_name, 'rb') as f:
        rfc = pkl.load(f)
    return rfc


def test_aggreate_prediction_results():
    dir = '/Users/wug/OneDrive - Oregon Health & Science University/FI_Network_Construction/new_pipeline/random_forest/predict_scores'
    test_results_dfs = []
    for i in range(3):
        file_name = Path(dir, 'predict_scores{}.csv'.format(i))
        test_result_df = pd.read_csv(file_name)
        test_results_dfs.append(test_result_df)
    aggreated_df = aggreate_prediction_results(test_results_dfs)
    print(aggreated_df.head())
    print('Shape: {}'.format(aggreated_df.shape))


def test_train_predict_df():
    training_file = '/Volumes/ssd/results/reactome-idg/fi-network-ml/feature_files/training/feature_matrix_041720.csv'
    test_file = '/Volumes/ssd/results/reactome-idg/fi-network-ml/feature_files/test/feature_test_matrix_051120.csv'
    working_dir = '/Users/wug/temp'
    postfix = 'test'
    down_sampling = 0.10
    test_result_df, precision_recall_df = train_rf(train_file=training_file,
                                                   test_file=test_file,
                                                   working_dir=working_dir,
                                                   postfix_name=postfix,
                                                   down_sampling=down_sampling)
    print('prediction results:\n{}'.format(test_result_df.head()))
    print('Shape: {}'.format(test_result_df.shape))

    print('\nperformance: \n{}'.format(precision_recall_df.head()))
    print('Shape: {}'.format(precision_recall_df.shape))

    # Test the prediction
    predict_file = '/Volumes/ssd/results/reactome-idg/fi-network-ml/feature_files/prediction/prediction_061820.csv'
    prediction_result = rf_predict(predict_file=predict_file,
                                  working_dir=working_dir,
                                  postfix_name=postfix,
                                  down_sampling=down_sampling)
    print('\nPredict results:\n{}'.format(prediction_result.head()))
    print('Shape: {}'.format(prediction_result))


if __name__ == '__main__':
    # test_aggreate_prediction_results()
    test_train_predict_df()
