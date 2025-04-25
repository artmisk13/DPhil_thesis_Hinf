# In-house Python codes used in Chapter 4: Identifying variants associated with NTHi through GWAS

import pandas as pd
import os
import numpy as np

## 1) pseudo-random number generator to select subset of sequence reads from NCBI SRA search

from random import sample

os.chdir('/working/dir')

infile = 'list_accession_all.txt' # list of all accession numbers

with open(infile, 'r') as f:
    content = f.readlines()

pool = []
for line in content:
    line = line.rstrip('\n')
    pool.append(line)

n = 135 # number of sequence reads to be chosen

samples = sample(pool, n)

outfile = 'list_accession_sampled.txt' # randomly chosen accession numbers

with open(outfile, 'w') as fo:
    for sample in samples:
        fo.write(f"{sample}\n")
fo.close()

### --- ###


## 2) Random allocation of samples to 200 GWAS subsets

os.chdir('/dir/to/dataset_file')

infile = 'genomes_phenotype.csv' # contains 2 columns: "genome ID" and "phenotype" (1: invasive, 0: non-invasive)
outdir = '/dir/for/subset'

df_nthi = pd.read_csv(infile)

for i in range(1,201): # generate 200 subset
    df_sample = df_nthi.groupby(by='phenotype', group_keys=False).apply(lambda x:x.sample(100)) # 100 samples of each phenotype group for every iteration
    df_sample = df_sample.reset_index(drop=True)
    df_sample.drop(df_sample.columns[[1]], axis=1, inplace=True)
    
    df_sample.to_csv(outdir+'phenotype_' + str(i) + '.txt', sep='\t', index=False)

### --- ###


## 3) Capsule region masking (2)

os.chdir('/dir/to/unitigs_output')

with open('nthi_unitigs_output/unitigs.txt', 'r') as f: # unitigs output from step (3) in gwas_code.sh
    content = f.readlines()
    
all_kmers = []
for line in content:
    line = line.split(' ')[0]
    all_kmers.append(line)

with open('caps_reg_ONLY_unitigs_output/unitigs.txt', 'r') as f: # unitigs output from step (9.2) in gwas_code.sh
    content = f.readlines()
    
capsreg_kmers = []
for line in content:
    line = line.split(' ')[0]
    capsreg_kmers.append(line)

line_num_ex = []
line_num_in = []
seq_in = []
for idx, kmer in enumerate(all_kmers):
    if kmer in capsreg_kmers:
        line_num_ex.append(idx)
    else:
        line_num_in.append(idx)
        seq_in.append(kmer)
        continue

with open('nthi_unitigs_output/unitigs_no_capsreg.txt', 'w') as fo:
    for idx in line_num_in:
        fo.write(f'{content[idx]}')
fo.close()

### --- ###


## 4) Chi-square test

from scipy.stats import chi2_contingency

os.chdir('/dir/to/kmer_mapping_output') # significant kmers presence/absence in all NTHi isolates

df_main = pd.read_csv('map_presence_absence.csv') # column 1: "genome ID", col 2: "phenotype", col 3-end: presence/absence of each significant kmer

kmers = df_main.columns.values.tolist()
del kmers[0:2]

summaries = []
for kmer in kmers:
    df_pivot = pd.pivot_table(data=df_main, values='genome ID', index=kmer, columns='phenotype', aggfunc='count')
    list1 = df_pivot.values.tolist()[0]
    list2 = df_pivot.values.tolist()[1]
    data = [list1, list2]
    phenotype_neg = list2[0]/(list2[0]+list1[0])*100
    phenotype_pos = list2[1]/(list2[1]+list1[1])*100
    stat, p, dof, expected = chi2_contingency(data)
    summary = [kmer, phenotype_pos, phenotype_neg, p, stat, dof]
    summaries.append(summary)

df_summaries = pd.DataFrame(summaries, columns=['kmer', '1', '0', 'pvalue', 'chisq_stat', 'dof'])
df_summaries.to_csv('kmers_chi_square.csv', index=False)

## Filter out kmers with p-value >= 0.1, kmers with p-value < 0.1:
	## kmer names stored in a text file: kmer_names.txt
	## gene hit names stored in a text file: kmer_hits.txt
		## e.g. kmer "cds-WP_044331820.1_4;init" and "cds-WP_044331820.1_7;init" were 2 different kmers, but both annotated to the same gene hit, that is "cds-WP_044331820.1"

### --- ###


## 5) Predictors/features selection 1: Check for multicollinearity

from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

os.chdir('/working/dir')

in_kmernames = '/dir/to/kmer_names.txt'
in_genehit = '/dir/to/kmer_hits.txt'
in_map = '/dir/to/kmer_mapping_output/map_presence_absence.csv'

with open(in_kmernames, 'r') as f:
    content = f.readlines()

kmers = []
for line in content:
    line = line.rstrip('\n')
    kmers.append(line)

with open(in_genehit, 'r') as f:
    content = f.readlines()

kmers_hits = []
for line in content:
    line = line.rstrip('\n')
    kmers_hits.append(line)

df_source = pd.read_csv(in_map)
df_main = df_source[['ID','phenotype'] + kmers]

features_tmp = []

kmer_clust = []
clust_result = []

for hit in kmer_hits:
    columns = [header for header in df_main.columns.tolist() if hit in header]
    df_tmp = df_main[columns]
    
    if len(df_tmp.columns) > 1: # there's a genehit that has only 1 kmer a/w invasiveness 
        corr_mtx_hit = df_tmp.corr()

        # hierarchal clustering to cluster predictor and choose representative
        dist_matrix = 1 - corr_mtx_hit

        linked = linkage(squareform(dist_matrix), method='average') # other methods: single, complete, weighted, centroid, median, ward

        threshold = 0.3 # equal to threshold of correlation coefficient of 0.7
        clusters = fcluster(linked, threshold, criterion='distance')

        cluster_dict = {}
        for col, cluster_id in zip(corr_mtx_hit.columns, clusters):
            cluster_dict.setdefault(cluster_id, []).append(col)

        for key in list(cluster_dict.keys()):
            tmp_list = cluster_dict[key]
            for kmer in tmp_list:
                kmer_clust.append(kmer)
                clust_result.append(key)
        
        representatives = [predictors[0] for predictors in cluster_dict.values()]

        df_in = df_main[representatives]
        corr_mtx_in = df_in.corr()
        corr_mtx_in.to_csv('prediction_multicollinearity_check_mtx-' + hit + '-included.csv')

        # save to main list
        features_tmp.append(representatives)
    
    else:
        features_tmp.append(columns)
        
# flatenned the main list
selected_features = [
    kmer for item in features_tmp for kmer in item
]

# to visually check if there is high corellation between features of different genehit

df_check = df_main[selected_features]
corr_mtx_check = df_check.corr()

corr_mtx_check.to_csv('prediction_multicollinearity_check_mtx_included_all.csv')

## Result: exclude more kmers if pair of kmer showed high correlation coefficient (> 0.7) --> list excluded kmers in a text file: excluded_kmer_names.txt

with open('excluded_kmer_names.txt', 'r') as f:
    content = f.readlines()

remove = []
for line in content:
    line = line.rstrip('\n')
    remove.append(line)

selected_features_new = [feature for feature in selected_features if feature not in remove]

tmp = set(selected_features_new)
selected_features_new = list(tmp)

with open('kmer_list_feature_selection_1.txt', 'w') as fo:
    for kmer in selected_features_new:
        fo.write(f'{kmer}\n')
fo.close()

### --- ###


## 6) Predictors/features selection 2: Importance metrics from Random Forest

from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.ensemble._forest import _generate_unsampled_indices
from rfpimp import *

df_sel1 = df_source[['ID','phenotype'] + selected_features_new]

# RandomForestClassifier to decide importances via column permutation

def mkdf(columns, importances):
    I = pd.DataFrame(data={'Feature':columns, 'Importance':importances})
    I = I.set_index('Feature')
    I = I.sort_values('Importance', ascending=False)
    return I

def oob_classifier_accuracy(rf, X_train, y_train):
    X = X_train.values
    y = y_train.values

    n_samples_bootstrap = len(X)
    n_samples = len(X)
    n_classes = len(np.unique(y))
    predictions = np.zeros((n_samples, n_classes))
    for tree in rf.estimators_:
        unsampled_indices = _generate_unsampled_indices(tree.random_state, n_samples, n_samples_bootstrap)
        tree_preds = tree.predict_proba(X[unsampled_indices, :])
        predictions[unsampled_indices] += tree_preds

    predicted_class_indexes = np.argmax(predictions, axis=1)
    predicted_classes = [rf.classes_[i] for i in predicted_class_indexes]

    oob_score = np.mean(y == predicted_classes)
    return oob_score

def permutation_importances(rf, X_train, y_train, metric):
    """
    Return importances from pre-fit rf; metric is function
    that measures accuracy or R^2 or similar. This function
    works for regressors and classifiers.
    """
    baseline = metric(rf, X_train, y_train)
    imp = []
    for col in X_train.columns:
        save = X_train[col].copy()
        X_train[col] = np.random.permutation(X_train[col])
        m = metric(rf, X_train, y_train)
        X_train[col] = save
        imp.append(baseline - m)
    return np.array(imp)

df_train, df_test = train_test_split(df_sel1, test_size=0.20)
features = selected_features_new + ['phenotype']

df_train = df_train[features]
df_test = df_test[features]

X_train, y_train = df_train.drop('phenotype',axis=1), df_train['phenotype']
X_test, y_test = df_test.drop('phenotype',axis=1), df_test['phenotype']

base_rf = RandomForestClassifier(n_estimators=100,
                            min_samples_leaf=5,
                            n_jobs=-1,
                            oob_score=True)

rf = clone(base_rf)
rf.fit(X_train, y_train)
oob = oob_classifier_accuracy(rf, X_train, y_train)

imp = permutation_importances(rf, X_train, y_train,
                              oob_classifier_accuracy)
I = mkdf(X_train.columns,imp)

# permutation importance is done 2x, the 2nd run with 'random' column to 'confuse' the approach - to test for consistence/stability

X_train2 = X_train.copy()
X_train2['random'] = np.random.random(size=len(X_train))
rf2 = clone(base_rf)
rf2.fit(X_train2, y_train)

imp2 = permutation_importances(rf2, X_train2, y_train,
                            oob_classifier_accuracy)
I2 = mkdf(X_train2.columns,imp2)

I.to_csv('feature_importance_random-forest_1_col-perm.csv')
I2.to_csv('feature_importance_random-forest_2_col-perm.csv')

## Result: compare two iterations and only feature with positive importance score in both were included in the logistic regression model
	## list of kmer names stored in a text file: selected_kmer_for_logres.txt

### --- ###


## 7) Logistic Regression

import statsmodels.api as sm
from sklearn.metrics import roc_auc_score, accuracy_score, roc_curve, confusion_matrix, ConfusionMatrixDisplay
import joblib

with open('selected_kmer_for_logres.txt', 'r') as f:
    content = f.readlines()
    
kmer_logres = []
for line in content:
    line = line.rstrip('\n')
    kmer_logres.append(line)

in_map = '/dir/to/kmer_mapping_output/map_presence_absence.csv'
df_source = pd.read_csv(in_map)
df_logres = df_source[['ID','phenotype'] + kmer_logres]

dependent_var = 'phenotype'
predictors = kmer_logres

def fit_logistic_regression(X_train, y_train):
    """
    Fits a logistic regression model using GLM and returns relevant metrics.
    """
    # Add constant for intercept
    X_train_const = sm.add_constant(X_train)
    
    # Fit the GLM model with Binomial family (logistic regression)
    model = sm.GLM(y_train, X_train_const, family=sm.families.Binomial()).fit()
    
    # Calculate pseudo R² (McFadden)
    ll_null = sm.GLM(
        y_train,
        sm.add_constant(pd.DataFrame({'const': np.ones(len(X_train))})),
        family=sm.families.Binomial()
    ).fit().llf
    ll_model = model.llf
    pseudo_r2 = 1 - (ll_model / ll_null)
    
    # Get AIC
    aic = model.aic
    
    # Get p-values
    pvalues = model.pvalues
    
    # Get model coefficients
    coefficients = model.params
    
    return model, pseudo_r2, aic, pvalues, coefficients

def evaluate_model(model, X_test, y_test):
    """
    Evaluates the fitted model on the test set and returns ROC AUC and Accuracy.
    """
    # Add constant to test data
    X_test_const = sm.add_constant(X_test)
    
    # Predict probabilities
    y_pred_prob = model.predict(X_test_const)
    
    # Predict class labels (using 0.5 threshold)
    y_pred = (y_pred_prob >= 0.5).astype(int)
    
    # Calculate ROC AUC
    roc_auc = roc_auc_score(y_test, y_pred_prob)
    
    # To make the ROC curve
    fpr, tpr, thresholds = roc_curve(y_test, y_pred_prob)
    
    # Calculate Accuracy
    accuracy = accuracy_score(y_test, y_pred)
    
    return roc_auc, fpr, tpr, accuracy, y_test, y_pred

# implement backward elimination to find the best set of predictors

def backward_elim(X, y, significance_level=0.05):
    """
    Perform backward elimination to select the best predictors based on AIC and accuracy.
    """
    # Split data into train and test sets
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )
    
    y_train = list(y_train)
    y_test = list(y_test)
    
    # Start with all predictors
    predictors = list(X.columns)
    best_model = None
    best_aic = float('inf')
    best_accuracy = 0

    while len(predictors) > 0:
        X_train_sub = X_train[predictors]
        X_test_sub = X_test[predictors]
        
        model, r2, aic, pvalues, coefficients = fit_logistic_regression(X_train_sub, y_train)
        roc_auc, fpr, tpr, accuracy, y_test, y_pred = evaluate_model(model, X_test_sub, y_test)
        
        # Check if this model is better
        if aic < best_aic or (aic == best_aic and accuracy > best_accuracy):
            best_model = model
            best_aic = aic
            best_accuracy = accuracy
            best_roc_score = roc_auc
            best_fpr = fpr
            best_tpr = tpr
            y_test_cm = y_test
            y_pred_cm = y_pred
        else:
            # If no improvement, stop elimination
            break
            
        # Find the predictor with the highest p-value: predictor elimination is based on the pvalue in each model
        pvalues = pvalues.drop('const')  # Drop intercept
        max_pvalue = pvalues.max()
        if max_pvalue > significance_level:
            # Remove the predictor with the highest p-value
            predictor_to_remove = pvalues.idxmax()
            predictors.remove(predictor_to_remove)
        else:
            # If all predictors are significant, stop elimination
            break
    
    return best_model, predictors, best_aic, best_accuracy, best_roc_score, best_fpr, best_tpr, y_test_cm, y_pred_cm

X = df_logres[predictors]
y = df_logres[dependent_var]

best_model, best_predictors, best_aic, best_accuracy, best_roc_score, best_fpr, best_tpr, \
y_test_cm, y_pred_cm = backward_elim(X, y)

# Visualise confusion matrix

cm = confusion_matrix(y_test_cm, y_pred_cm)
disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=None)
disp.plot(cmap=plt.cm.Blues)

plt.show()

# Plot the ROC curve
plt.figure(figsize=(8, 6))
plt.plot(best_fpr, best_tpr, label=f"ROC Curve (AUC = {best_roc_score:.4f})")
plt.plot([0, 1], [0, 1], linestyle="--", color="gray", label="Random Guess")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve")
plt.legend()
plt.grid()
plt.savefig('Model_Logistic_Regression_ROC.svg', format='svg')
plt.show()

# Save the best model to load and reuse in the future!

joblib.dump(best_model, 'Model_Logistic_Regression_NTHI_invasive_correct.pkl')

# Save other important information: model summary, list of predictors, the ROC curve

tab_sum_1 = best_model.summary().tables[0].as_html()
df_sum_1 = pd.read_html(tab_sum_1, header=None, index_col=0)[0]

tab_sum_2 = best_model.summary().tables[1].as_html()
df_sum_2 = pd.read_html(tab_sum_2, header=0, index_col=0)[0]

df_sum_1.to_csv('Model_Logistic_Regression_summary1.csv')
df_sum_2.to_csv('Model_Logistic_Regression_summary2.csv')

### --- ###