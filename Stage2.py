import pandas as pd
import numpy as np  
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import f1_score, average_precision_score, precision_recall_curve, matthews_corrcoef
from sklearn.neighbors import LocalOutlierFactor
from scipy.stats import norm
from imblearn.ensemble import BalancedRandomForestClassifier

df = pd.read_table("data/Halgreen 2009 - Table 6.txt", sep=" ") 

for threshold in [0, 1]:
    df["label"] = (df["categorya"] > threshold).astype(bool)
    
    features = ["Dscore","SScore","size","enclosure","philic","phobic"]
    X = df[features]
    y = df["label"]
    
    lof = LocalOutlierFactor()
    outlier_scores = lof.fit_predict(X.values)
    non_outlier_indices = [i for i, score in enumerate(outlier_scores) if score == 1]
    X_cleaned = X.iloc[non_outlier_indices]
    y_cleaned = y.iloc[non_outlier_indices]    
    
    kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    metrics = []
    
    for train_index, test_index in kf.split(X_cleaned, y_cleaned):
        X_train, X_test = X_cleaned.iloc[train_index], X_cleaned.iloc[test_index]
        y_train, y_test = y_cleaned.iloc[train_index], y_cleaned.iloc[test_index]
       
        rf =  BalancedRandomForestClassifier(random_state=42)
        rf.fit(X_train, y_train) 
               
        y_pred_rf = rf.predict(X_test)        
        y_pred_rf_proba = rf.predict_proba(X_test)[:, 1]
        
        f1_rf = f1_score(y_test, y_pred_rf)
        mcc_rf =  matthews_corrcoef(y_test, y_pred_rf)
        z_score = mcc_rf / (1 - mcc_rf**2) ** 0.5
        mcc_rf_pvalue = 2 * (1 - norm.cdf(abs(z_score)))
        app_rf = average_precision_score(y_test, y_pred_rf_proba)
        
        metrics.append([mcc_rf, mcc_rf_pvalue, f1_rf, app_rf, threshold, "RF"])
    
    metrics_df = pd.DataFrame(metrics, columns=["MCC", "MCC p-value", "F1 Score", "Average Precision Score", "Threshold", "Model"])
    metrics_df[["MCC", "Threshold"]].to_csv(f"outputs/cross_validation_data_threshold_{threshold}.csv")

confidence_cut_level = 0.5
df["y"] = df["categorya"]>0

lof = LocalOutlierFactor()
outlier_scores = lof.fit_predict(X.values)
non_outlier_indices = [i for i, score in enumerate(outlier_scores) if score == 1]
X_cleaned = X.iloc[non_outlier_indices]
y_cleaned = y.iloc[non_outlier_indices]  

rf = BalancedRandomForestClassifier(random_state=42) 
rf.fit(X_cleaned.values, y_cleaned)
   
lof2 = LocalOutlierFactor(novelty=True)
lof2.fit(X_cleaned.values)

df_targets = pd.read_csv("outputs/sitemap_results_with_quality_metrics.csv")
df_targets["UniProtKB Gene Name ID"] = df_targets["Entry"]

df_targets2 = df_targets[(df_targets["Entire Region Percent High Quality"]>=confidence_cut_level) & 
                         (df_targets["Percent High Quality"]>=confidence_cut_level)]

X_targets = df_targets2[["Dscore","SiteScore","size","enclosure","philic","phobic"]]
outlier_scores = lof2.predict(X_targets.values)
X_good_targets = X_targets[outlier_scores > -1]  

rf_probs = rf.predict_proba(X_good_targets.values)

df_targets2 = df_targets2[outlier_scores > -1] 
df_targets2["rf_probs"] = rf_probs[:,1]  
df_targets3 = df_targets2[df_targets2["rf_probs"]>=confidence_cut_level]

df_gene_map = pd.read_csv("data/uniprot_ensembl_map.txt")
df_merged = df_targets3.merge(df_gene_map, on="UniProtKB Gene Name ID")
df_hq = df_merged.drop_duplicates()

df_hq[["Gene name","residues", "Entire Region Percent High Quality", "Percent High Quality"]].drop_duplicates().to_csv("outputs/final_filtered_sites.csv", index=None)