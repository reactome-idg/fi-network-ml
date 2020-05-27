"""
Summary of modeling work that helped or evaluated final results.
most work was evaluated on internal server vs ExaCLoud. 
"""
import pandas 
import numpy 
import pickle
import sklearn
from sklearn.model_selection import train_test_split, cross_val_predict
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.metrics import roc_auc_score, roc_curve, auc, classification_report, confusion_matrix
from sklearn.model_selection import cross_val_score
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn import linear_model
from sklearn.naive_bayes import CategoricalNB, BernoulliNB
from sklearn.svm import SVC


""" 
------------------------------------------------------------------
Wrangle Data 
------------------------------------------------------------------
ExaCloud volume /home/exacloud/lustre1/WuLab/idgml/:/opt/
"""
df = pandas.read_csv("/opt/feature_matrix_041720.csv")
df.GenePair = df.GenePair.str.replace('\t', '_')
df.set_index('GenePair', inplace=True)

df.shape 
sum(df.FI == 0)
# 9612200
sum(df.FI == 1)
# 96122

X = df[df.columns[~df.columns.isin(['FI'])]]
y = df['FI']

# percentage of data set aside for testing is usually around 80/20 or 70/30. we choose 0.25
# Splitting with same proportions of y (0's and 1's) that is handled via train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=42)


""" 
------------------------------------------------------------------
Initial random forest classifier - EDA 
------------------------------------------------------------------
"""
# clf = RandomForestClassifier(n_estimators=1000) 
clf = RandomForestClassifier(n_estimators=100, random_state=42)

model = clf.fit(X_train, y_train)
score = clf.score(X_test, y_test)

# cross validation mean 
cv = cross_val_score(clf, X_train, y_train, cv=10, scoring="accuracy")
numpy.mean(cv)
# 0.991 accuracy score 

""" 
------------------------------------------------------------------
Dimention reduction using PCA - EDA 
------------------------------------------------------------------
"""
pca = PCA(n_components = 0.95) 
pca.fit(X)
Xp = pca.transform(X)
Xp.shape
# 85 out of 106 features (85 / 106 * 100 = 80.18867924528303 % of data is explained via features)


clf2 = RandomForestClassifier(n_estimators=1000)
Xp_train, Xp_test, yp_train, yp_test = train_test_split(Xp, y, test_size=0.25, random_state=42)

clf2.fit(Xp_train, yp_train)
clf2.score(Xp_test, yp_test)

# cross validation mean 
cv2 = cross_val_score(clf2, Xp_train, yp_train, cv=10, scoring="accuracy")
numpy.mean(cv2)
# 0.9920 # very small improvement 

"""
------------------------------------------------------------------
BernoulliNB or categoricalNB - EDA 
------------------------------------------------------------------
"""
clf_b = BernoulliNB()

clf_b.fit(X_train, y_train)
prd = clf_b.predict(X_test)

metrics.accuracy_score(y_test, prd)
# 0.9818
roc_auc_score(y_test, prd)
# 0.6832

# ---------------------------------
clf_c = CategoricalNB()

clf_c.fit(X_train, y_train)
prd = clf_c.predict(X_test)

metrics.accuracy_score(y_test, prd)
# 0.9827
roc_auc_score(y_test, prd)
# 0.6793

""" 
------------------------------------------------------------------
Hyper parameter tuning - gridSearch (on smaller subsets of data - memory & time)
May have to check params individually on best model fit 
------------------------------------------------------------------
"""
param_grid = {
 'class_weight': ['balanced', 'balanced_subsample', None],
 'max_depth' : [2, 4, 6, 10, None],
 'max_features': ['auto', 'sqrt', 'log2'],
 'min_samples_leaf': [1, 2, 4],
 'min_samples_split': [2, 4, 10],
 'n_estimators': [100, 200, 500, 1000], 
 }

l = [] 
for i in range(10): 
	df2 = pandas.concat([df[df.FI == 1].sample(n=100), df[df.FI == 0].sample(n=10000)])
	X2 = df2[df2.columns[~df2.columns.isin(['FI'])]]
	y2 = df2['FI']

	X2_train, X2_test, y2_train, y2_test = train_test_split(X2, y2, test_size=0.25, random_state=42)

	grid_search = GridSearchCV(estimator=clf2, param_grid=param_grid, cv=10, n_jobs=-1, verbose=2)
	grid_search.fit(X2_train, y2_train)
	tunned_model = grid_search.best_estimator_
	l.append(tunned_model)


""" 
------------------------------------------------------------------
Balanced RF with tuned hyper-params 
------------------------------------------------------------------
"""
clf = RandomForestClassifier(bootstrap=True, ccp_alpha=0.0, class_weight='balanced',
                       criterion='gini', max_depth=10, max_features='auto',
                       max_leaf_nodes=None, max_samples=None,
                       min_impurity_decrease=0.0, min_impurity_split=None,
                       min_samples_leaf=1, min_samples_split=2,
                       min_weight_fraction_leaf=0.0, n_estimators=200,
                       n_jobs=None, oob_score=False, random_state=42,
                       verbose=0, warm_start=False)

clf.fit(X_train, y_train)

prd = clf.predict(X_test)
clf.score(X_test, y_test)
roc_auc_score(y_test, prd)
cm = confusion_matrix(y_test, prd)

"""
roc: 0.8959475290303648
array([[2150218,  252891],
       [   2466,   21506]])

"""

cv = cross_val_score(clf, X_train, y_train, cv=10, scoring="f1")
numpy.mean(cv)
"""array([0.14502133, 0.14443901, 0.14483052, 0.14585848, 0.14420461,
       0.14510427, 0.14499293, 0.14478786, 0.145924  , 0.14537118])"""

cv_roc = cross_val_score(clf, X_train, y_train, cv=10, scoring="roc_auc")
numpy.mean(cv_roc)
"""array([0.96045308, 0.96034245, 0.96159635, 0.9623916 , 0.96091947,
       0.96072559, 0.96100835, 0.95936615, 0.96259128, 0.96048622])"""

pickle.dump(clf, open( "/opt/data/output/clf_RF_balanced_89ROC.pkl", "wb" ) )
"""
------------------------------------------------------------------
Try linear svc
------------------------------------------------------------------
"""
clf = LinearSVC(random_state=42,  dual=True, max_iter=1000, class_weight='balanced')
clf.fit(X_train, y_train)

prd = clf.predict(X_test)
clf.score(X_test, y_test)
# 0.8885183477601283
roc_auc_score(y_test, prd)
# 0.903312511832865

cv_roc = cross_val_score(clf, X_train, y_train, cv=10, scoring="roc_auc")
numpy.mean(cv_roc)

cv = cross_val_score(clf, X_train, y_train, cv=10, scoring="f1")
numpy.mean(cv)

pickle.dump(clf, open( "/opt/data/output/clf_LinearSVC_balanced_90ROC.pkl", "wb" ) )
"""
------------------------------------------------------------------
Run models on validation data 
------------------------------------------------------------------
"""
df = pandas.read_csv("/opt/feature_test_matrix_051120.csv")
df.GenePair = df.GenePair.str.replace('\t', '_')
df.set_index('GenePair', inplace=True)
df.shape 
sum(df.FI == 0)
# 2196398
sum(df.FI == 1)
# 26076

X_val = df[df.columns[~df.columns.isin(['FI'])]]
y_val = df['FI']

# -----------------------------------
clf = pickle.load( open( "/opt/output/clf_RF_balanced_89ROC.pkl", "rb" ) )
prd = clf.predict(X_val)
clf.score(X_val, y_val)
roc_auc_score(y_val, prd)
# 0.8865219292534692

feature_importance = clf.feature_importances_
feature_importance = 100 * (feature_importance / feature_importance.max())

df_importance = pandas.DataFrame({'features':df[df.columns[~df.columns.isin(['FI'])]].columns.values, 'RF_importance':feature_importance})
pandas.to_csv("/opt/output/RF_feature_importances_validation_data.csv", index=False)

# -----------------------------------
clf = pickle.load( open( "/opt/output/clf_LinearSVC_balanced_90ROC.pkl", "rb" ) )
prd = clf.predict(X_val)
clf.score(X_val, y_val)
roc_auc_score(y_val, prd)
# 0.8877285851659766
