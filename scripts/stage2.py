import pandas as pd
import numpy as np  
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import f1_score, average_precision_score, matthews_corrcoef, brier_score_loss
from scipy.stats import norm
from catboost import CatBoostClassifier
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.calibration import calibration_curve

df = pd.read_table("../data/Halgreen 2009 - Table 6.txt", sep=" ") 

for threshold in [0, 1]:
    df["label"] = (df["categorya"] > threshold).astype(bool)
    
    features = ["Dscore", "SScore", "size", "enclosure", "philic", "phobic"]
    X = df[features]
    y = df["label"]
    
    rkf = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=42)
    metrics = []
    
    for train_index, test_index in rkf.split(X, y):
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]

        pos_weight = (y_train == 0).sum() / (y_train == 1).sum()

        model = CatBoostClassifier(
            iterations=1000,
            learning_rate=0.05,
            depth=6,
            loss_function='Logloss',
            class_weights=[1.0, pos_weight],
            verbose=0,
            random_state=42
        )
        model.fit(X_train, y_train) 
               
        y_pred = model.predict(X_test)        
        y_proba = model.predict_proba(X_test)[:, 1]
        
        f1 = f1_score(y_test, y_pred)
        mcc = matthews_corrcoef(y_test, y_pred)
        z_score = mcc / (1 - mcc**2) ** 0.5
        mcc_pvalue = 2 * (1 - norm.cdf(abs(z_score)))
        ap = average_precision_score(y_test, y_proba)
        brier = brier_score_loss(y_test, y_proba)
        
        metrics.append([mcc, mcc_pvalue, f1, ap, brier, threshold, "CatBoost"])
    
    metrics_df = pd.DataFrame(metrics, columns=["MCC", "MCC p-value", "F1 Score", "Average Precision Score", "Brier Score", "Threshold", "Model"])
    metrics_df.to_csv(f"../outputs/cross_validation_data_threshold_{threshold}.csv", index=False)

# === Boxplot + Stripplot for MCC Distributions ===
df0 = pd.read_csv("../outputs/cross_validation_data_threshold_0.csv")
df1 = pd.read_csv("../outputs/cross_validation_data_threshold_1.csv")
df_all = pd.concat([df0, df1], ignore_index=True)
df_all["Threshold"] = df_all["Threshold"].astype(str)

plt.rcParams['font.family'] = 'Arial'
plt.figure(figsize=(6, 5))
sns.boxplot(data=df_all, x="Threshold", y="MCC", showfliers=False)
sns.stripplot(data=df_all, x="Threshold", y="MCC", color=".25", alpha=0.6, jitter=0.2)

plt.title("Cross validation of pocket druggability models", fontsize=16)
plt.xlabel("Target > threshold", fontsize=14)
plt.ylabel("Matthews correlation coefficient (MCC)", fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.savefig("../outputs/mcc_box_strip_plot.png", dpi=300)
plt.show()

# === Dynamically Select Best Threshold ===
mean_mcc_by_threshold = df_all.groupby("Threshold")["MCC"].mean()
best_threshold = float(mean_mcc_by_threshold.idxmax())
print(f"Best threshold based on average MCC: {best_threshold}")

# === Downstream logic ===
confidence_cut_level = 0.5
df["y"] = df["categorya"] > best_threshold

X_final = df[features]
y_final = df["y"]

final_pos_weight = (y_final == 0).sum() / (y_final == 1).sum()

final_model = CatBoostClassifier(
    iterations=1000,
    learning_rate=0.05,
    depth=6,
    loss_function='Logloss',
    class_weights=[1.0, final_pos_weight],
    verbose=0,
    random_state=42
)
final_model.fit(X_final, y_final)

#Note the file below is excluded from the online repository to comply with Schrodinger's terms
df_targets = pd.read_csv("../schrodinger/candidate_input_sitemap_data.csv")
df_targets["UniProtKB Gene Name ID"] = df_targets["Entry"]

df_targets2 = df_targets[
    (df_targets["Entire Region Percent High Quality"] >= confidence_cut_level) & 
    (df_targets["Percent High Quality"] >= confidence_cut_level)
]

print("AF quality filter - unique genes remaining:", df_targets2["UniProtKB Gene Name ID"].unique().shape[0])

X_targets = df_targets2[["r_sitemap_Dscore","r_sitemap_SiteScore","i_sitemap_size","r_sitemap_enclosure","r_sitemap_philic","r_sitemap_phobic"]]
catboost_probs = final_model.predict_proba(X_targets.values)

df_targets2["cat_probs"] = catboost_probs[:, 1]  
df_targets3 = df_targets2[df_targets2["cat_probs"] >= confidence_cut_level]

df_gene_map = pd.read_csv("../data/uniprot_ensembl_map.txt")
df_merged = df_targets3.merge(df_gene_map, on="UniProtKB Gene Name ID")
df_hq = df_merged.drop_duplicates()

df_hq[["Gene name_x", "s_sitemap_residues", "Entire Region Percent High Quality", "Percent High Quality", "cat_probs"]].drop_duplicates().to_csv("../outputs/final_filtered_sites.csv", index=None)