"""
This script is used to analyze F1 score based on precision/recall file geneated originally.
"""

import pandas as pd

DIR_NAME = '../../results/features_check'
PRECISION_RECALL_FILE = DIR_NAME + '/precision_recall.csv'

pr_df = pd.read_csv(PRECISION_RECALL_FILE)
print(pr_df)

