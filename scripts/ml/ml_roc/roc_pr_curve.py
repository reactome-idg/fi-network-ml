# This script generates a precision-recall curve using plotly
# Test matrix generated via FeatureFileGenerator.jar
# Triplicate random forest classifiers are used on the test matrix to generate prediction scores
# The mean of the prediction scores are used for precision/recall/F1 calculations

import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve
from sklearn.metrics import RocCurveDisplay
from sklearn.metrics import precision_recall_curve
import plotly.graph_objects as go 


test_matrix_file = '/Users/conleyp/Projects/Reactome/fi_network_build/results/2023/generated_test_matrix_17Jan2023.csv'

# pickled random forest classifiers
clf1_file = '/Users/conleyp/Projects/reactome-idg/fi-network-ml/run/results/predict_scores/trained_rfc_0.pkl'
clf2_file = '/Users/conleyp/Projects/reactome-idg/fi-network-ml/run/results/predict_scores/trained_rfc_1.pkl'
clf3_file = '/Users/conleyp/Projects/reactome-idg/fi-network-ml/run/results/predict_scores/trained_rfc_2.pkl'

# prediction scores files
prd1_file = '/Users/conleyp/Projects/reactome-idg/fi-network-ml/run/results/predict_scores/predict_scores0.csv'
prd2_file = '/Users/conleyp/Projects/reactome-idg/fi-network-ml/run/results/predict_scores/predict_scores1.csv'
prd3_file = '/Users/conleyp/Projects/reactome-idg/fi-network-ml/run/results/predict_scores/predict_scores2.csv'

# Load Test matrix
test_df = pd.read_csv(test_matrix_file)
test_df.set_index('GenePair', inplace=True)

X = test_df[test_df.columns[~test_df.columns.isin(['FI'])]]
y = test_df['FI']
#del(df)

clf1 = pkl.load(open(clf1_file, 'rb'))
clf2 = pkl.load(open(clf2_file,'rb'))
clf3 = pkl.load(open(clf3_file,'rb'))

prob1 = clf1.predict_proba(X)
prob2 = clf2.predict_proba(X)
prob3 = clf3.predict_proba(X)

pd.DataFrame(prob1).to_csv(prd1_file)
pd.DataFrame(prob2).to_csv(prd2_file)
pd.DataFrame(prob3).to_csv(prd3_file)

scores0 = pd.read_csv(prd1_file)
scores1 = pd.read_csv(prd2_file)
scores2 = pd.read_csv(prd3_file)

dfp = pd.DataFrame({'rep0':scores0['1'], 'rep1':scores1['1'], 'rep2':scores2['1']})
dfp['mean'] = dfp.apply(np.mean, axis=1, raw=True)
dfp['stdev'] = dfp.apply(np.std, axis=1, raw=True)
rsdp = (dfp['stdev']/dfp['mean'])*100

s0=scores0[['1']].rename({'1':'1_rep1'}, axis=1)
s1=scores1[['1']].rename({'1':'1_rep2'}, axis=1)
s2=scores2[['1']].rename({'1':'1_rep3'}, axis=1)

combined = s0.merge(s1.merge(s2, left_index=True, right_index=True), left_index=True, right_index=True)

combined['1_mean'] = combined[['1_rep1','1_rep2','1_rep3']].apply(np.mean, axis=1, raw=True)
combined['1_min'] = combined[['1_rep1','1_rep2','1_rep3']].apply(np.min, axis=1, raw=True)
combined['1_max'] = combined[['1_rep1','1_rep2','1_rep3']].apply(np.max, axis=1, raw=True)
combined['1_std'] = combined[['1_rep1','1_rep2','1_rep3']].apply(np.std, axis=1, raw=True)
combined['1_%CV'] = combined['1_std']/combined['1_mean']*100
combined['Gene pair'] = y.index


# Create and save ROC Curve
RocCurveDisplay.from_predictions(y, combined['1_mean'])
plt.title('ROC Curve using mean prediction score')
#plt.show()
plt.savefig('ROC_mean_scores_20220123_remove.pdf')


# Create and save precision-recall curve HTML
df_pr = precision_recall_curve(y, combined['1_mean'])
df_pr = pd.DataFrame(df_pr).T
df_pr.columns = ['Precision', 'Recall', 'Thresholds']
df_pr['F1 Score'] = 2 * (df_pr['Precision'] * df_pr['Recall']) / (df_pr['Precision'] + df_pr['Recall'])


fig = go.Figure()
fig.add_trace(go.Scatter(y=df_pr['Precision'], x=df_pr['Thresholds'],
                    #mode='lines',
                    name='Precision'))

fig.add_trace(go.Scatter(y=df_pr['Recall'], x=df_pr['Thresholds'],
                    #mode='lines',
                    name='Recall'))
fig.add_trace(go.Scatter(y=df_pr['F1 Score'], x=df_pr['Thresholds'],
                    #mode='lines',
                    name='F1 Score'))
fig.write_html('precision_recall_f1_20220123_remove.html')