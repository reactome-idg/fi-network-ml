"""
This script is used to analyze F1 score based on precision/recall file geneated originally.
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import numpy as np
import sys
from sklearn.metrics import roc_curve, auc, precision_score
import pandas as pd


def analyze_rf_f1(precision_recall_file: str):
    """
    The following code is used to plot the performance measure for the trained Random Forest
    """
    # dir_name = '../../results/features_check'
    # precision_recall_file = dir_name + '/precision_recall.csv'
    # Try to get it from sys
    if precision_recall_file is None:
        precision_recall_file = sys.argv[1]
    pr_df = pd.read_csv(precision_recall_file)
    html_file_name = precision_recall_file.replace('.csv', '.html')
    plot_rf_performance(pr_df, html_file=html_file_name)


def analyze_auc(test_feature_file: str,
                test_result_file: str):
    # Load the files first
    test_feature_df = pd.read_csv(test_feature_file)
    test_result_df = pd.read_csv(test_result_file)
    fpr, tpr, thresholds = roc_curve(test_feature_df['FI'],
                                     test_result_df['mean'])

    df = pd.DataFrame({'Threshold': thresholds, 'FPR': fpr, 'TPR': tpr})
    print('ROC data frame: \n{}'.format(df.head()))
    df_file = test_result_file.replace('.csv', '_roc.csv')
    df.to_csv(df_file)

    roc_auc = auc(fpr, tpr)

    plt.figure()
    plt.plot(fpr, tpr, color='darkorange', lw=2,
             label='ROC curve (area = %0.3f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic')
    plt.legend(loc="lower right")

    # Save the plot to a file
    fig_file = test_result_file.replace('.csv', '_roc.png')
    plt.savefig(fig_file, bbox_inches='tight')

    plt.show()


def plot_rf_performance(precision_recall_df: pd.DataFrame,
                           need_plot: bool = True,
                           html_file: str = None):
    # Calculate F1 score based on precision and recall
    if 'F1' not in precision_recall_df.columns:
        calculate_f1(precision_recall_df)
    precision_recall_df.rename(
        columns={'Thresholds': 'Threshold'}, inplace=True)
    print(precision_recall_df.head())
    pr_df_plot = create_plot_df(precision_recall_df)
    if need_plot:
        plot = plot_measures(pr_df_plot)
        # Plot the range we want to avoid any exploration
        plot.set(xlim=(0.23559, 1.0),
                 title='Performance of Trained Random Forest')
        plt.show()

    if html_file:
        # Plot plotly
        fig = px.line(pr_df_plot,
                      x='Threshold',
                      y='Value',
                      color='Measure',
                      title='Performance of Trained Random Forest')
        fig.show()
        # Get the output file
        print('HTML file of the F1 plot: {}'.format(html_file))
        fig.write_html(html_file)

    return pr_df_plot


def calculate_f1(pr_df):

    def _calculate_f1(row):
        precision = row['Precision']
        recall = row['Recall']
        # Check for division by zero
        if precision == 0 or recall == 0:
            return np.nan  # Return NaN if either precision or recall is zero
        return 2 * (precision * recall) / (precision + recall)

    pr_df['F1'] = pr_df.apply(_calculate_f1, axis=1)


def plot_measures(pr_df):
    # Plot the results having three curves
    plot = sns.lineplot(x='Threshold', y='Value', hue='Measure', data=pr_df)
    return plot


def create_plot_df(pr_df):
    # Need to melt the dataframe to plot
    pr_df_melt = pd.melt(pr_df,
                         id_vars=['Threshold'],
                         var_name='Measure',
                         value_vars=['Precision', 'Recall', 'F1'],
                         value_name='Value')
    print(pr_df_melt.head())
    return pr_df_melt


def analyze_nbc_f1(need_plot=True):
    """
    The following function is used to analyze the performance the old Naive Bayers classifier
    """
    dir_name = '../../../../FINetworkBuild/results/2021'
    file_name = dir_name + '/ROC_100_122921.txt'
    # This is recorded in the combined logging
    total_neg_pairs = 3671970
    total_pos_pairs = 36811
    df = pd.read_csv(file_name, sep='\t')
    print(df.head())
    df.rename(columns={'Cutoff': 'Threshold',
              'True_Positive_Rate': 'Recall'}, inplace=True)
    # Need to calculate precision
    df['Precision'] = df.apply(lambda row: row['Recall'] * total_pos_pairs / (row['Recall'] * total_pos_pairs + row['False_Positive_Rate'] * total_neg_pairs),
                               axis=1)
    calculate_f1(df)
    print(df.head())
    out_file_name = '../../../../FINetworkBuild_RF/results/2022/NBC_ROC_100_122921_020423.txt'
    df.to_csv(out_file_name, sep='\t')
    df_plot = create_plot_df(df)
    if need_plot:
        plot = plot_measures(df_plot)
        plot.set(title='Performance of Trained Naives Bayer Classifier')
    return df_plot


def plot_two_f1(measure='F1'):
    """
    This function is used to plot two F1 together
    :return:
    """
    rf_df = analyze_rf_f1(False)
    rf_df['Model'] = 'RF'
    nbc_df = analyze_nbc_f1(False)
    nbc_df['Model'] = 'NBC'
    merged_df = pd.concat([rf_df, nbc_df], axis=0)
    # Do a filtering
    if measure is not None:
        merged_df = merged_df[merged_df['Measure'] == measure]
        plot = sns.lineplot(x='Threshold', y='Value',
                            data=merged_df, hue='Model')
        plot.set(ylabel=measure, title=measure +
                 ' Plot of Random Forest and Naive Bayes Classifier')
    else:
        # Plot all measures
        merged_df['Model_Measure'] = merged_df.apply(
            lambda row: row['Model'] + '_' + row['Measure'], axis=1)
        plot = sns.lineplot(x='Threshold', y='Value',
                            hue='Model_Measure', data=merged_df)
        plot.set(title='Performance Plot of Trained RF and NBC')
    return merged_df


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print('Usage: python {precision_recall_file.csv} {test_feature_file.csv} {test_result_file.csv}')
        sys.exit(1)
    analyze_rf_f1(sys.argv[1])
    if len(sys.argv) >= 4:  # Expect to run auc
        analyze_auc(test_feature_file=sys.argv[2],
                    test_result_file=sys.argv[3])
