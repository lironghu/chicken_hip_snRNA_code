import scanpy as sc
import pandas as pd
import numpy as np
import xgboost as xgb
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import roc_curve, roc_auc_score

os.chdir("E:/KIZ/12.RJF_project/XGBoost")

# 
FILE_PATH = "E:/KIZ/12.RJF_project/chicken_hip_scRNA_harmony_annotation_8samples_final.h5ad"
TARGET_CELLTYPE = 'IPC1'  #IPC2
N_REPEATS = 1000

# 
adata = sc.read_h5ad(FILE_PATH)
adata_subset = adata[adata.obs['Celltype'] == TARGET_CELLTYPE].copy()
X = adata_subset.X
le = LabelEncoder()
y = le.fit_transform(adata_subset.obs['Group'])
feature_names = adata_subset.var_names.tolist()

#Iterative Evaluation
all_auc_scores = []
all_importances = pd.DataFrame(index=feature_names)
tprs = []
base_fpr = np.linspace(0, 1, 101)

for i in range(N_REPEATS):
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.25, random_state=1000 + i, stratify=y
    )

    model = xgb.XGBClassifier(
        n_estimators=100, max_depth=5, learning_rate=0.1,
        objective='binary:logistic', eval_metric='logloss',
        tree_method='gpu_hist', predictor='gpu_predictor', random_state=42
    )
    model.fit(X_train, y_train)

    probs = model.predict_proba(X_test)[:, 1]
    all_auc_scores.append(roc_auc_score(y_test, probs))
    
    fpr, tpr, _ = roc_curve(y_test, probs)
    tprs.append(np.interp(base_fpr, fpr, tpr))
    tprs[-1][0] = 0.0
    
    all_importances[f'run_{i}'] = model.feature_importances_

# 
auc_df = pd.DataFrame({'Run': range(1, N_REPEATS + 1), 'AUC': all_auc_scores})
auc_df.to_csv("./IPC1_1000_auc_scores.csv", index=False)

importance_summary = all_importances.sum(axis=1).sort_values(ascending=False)
top_20_genes = pd.DataFrame({
    'gene': importance_summary.index[:20],
    'total_importance': importance_summary.values[:20],
    'relative_importance_%': (importance_summary.values[:20] / importance_summary.sum()) * 100
})
top_20_genes.to_csv("./IPC1_top_20_feature_importance.csv", index=False)

#Visualization
mean_auc, std_auc = np.mean(all_auc_scores), np.std(all_auc_scores)
mean_tprs = np.array(tprs).mean(axis=0)
std_tprs = np.array(tprs).std(axis=0)

plt.figure(figsize=(18, 5))

# Plot 1: ROC
plt.subplot(1, 3, 1)
plt.plot(base_fpr, mean_tprs, color='blue', label=f'Mean ROC (AUC={mean_auc:.3f}±{std_auc:.3f})')
plt.fill_between(base_fpr, np.maximum(mean_tprs-std_tprs, 0), np.minimum(mean_tprs+std_tprs, 1), color='grey', alpha=0.3)
plt.plot([0, 1], [0, 1], 'r--')
plt.xlabel('FPR'); plt.ylabel('TPR'); plt.legend()

# Plot 2: Importance
plt.subplot(1, 3, 2)
sns.barplot(x='total_importance', y='gene', data=top_20_genes, palette='viridis')
plt.title('Top 20 Genes')

# Plot 3: AUC Hist
plt.subplot(1, 3, 3)
sns.histplot(all_auc_scores, kde=True)
plt.title('AUC Distribution')

plt.tight_layout()
plt.savefig("./Combined_Analysis_IPC1.pdf")
plt.show()