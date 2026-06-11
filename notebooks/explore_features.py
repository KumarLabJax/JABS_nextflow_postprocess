# %%
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt

# %%
nextflow_features_dir = "/projects/kumar-lab/USERS/nguyetu/SING-grant/postNextflow/final_nextflow_feature_data"
# %%
feature_dfs = {}
for file in list(Path(nextflow_features_dir).rglob("*.csv")):
    feature_df = pd.read_csv(file)
    feature_df = feature_df.drop(columns=['nextflow_version'])
    feature_dfs[file.stem] = feature_df

# %%
merged_df = feature_dfs['morphometrics'] \
    .merge(feature_dfs['gait_final'], on="NetworkFilename", how="inner") \
    .merge(feature_dfs['JABS_features_final'], on="NetworkFilename", how="inner")
# %%
# Remove very low variance col
col_var = merged_df.var(numeric_only=True) / (abs(merged_df.mean(numeric_only=True)) + 1e-6)
low_var_cols = col_var[col_var <= 1e-6].index
print(low_var_cols)
filtered_df = merged_df.drop(columns=low_var_cols)
# %%
# Remove columns where more than half the samples have NA
filtered_df = filtered_df.drop(columns = filtered_df.columns[filtered_df.isna().mean() > 0.5])
# %%
## The only non-numeric column should be NetworkFilename
filtered_df.select_dtypes(exclude="number").columns

# %%
filtered_df = filtered_df.set_index('NetworkFilename')
filtered_df

# %%
## Initial PCA
def run_pca_pipeline(df, n_components=None):
    # Impute numeric features
    num_df = df.select_dtypes(include="number")
    imputer = SimpleImputer(strategy="median")
    X_imputed = imputer.fit_transform(df)

    # Scale
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_imputed)

    # PCA
    pca = PCA(n_components=n_components)
    X_pca = pca.fit_transform(X_scaled)

    # Variance explained by each PC
    explained_var = pca.explained_variance_ratio_
    expl_var_dict = {f"PC{i+1}": explained_var[i] for i in range(len(explained_var))}

    # Build PCA dataframe
    pca_df = pd.DataFrame(X_pca, 
                          index = df.index,
                          columns=[f"PC{i+1}" for i in range(X_pca.shape[1])])

    return pca_df, expl_var_dict

# %%
pca_df, expl_var_dict = run_pca_pipeline(filtered_df)
pca_df['Sex'] = pca_df.index.map(lambda p: (Path(p).name).split("_")[1].lower())
pca_df['Strain'] = pca_df.index.map(lambda p: (Path(p).name).split("_")[2])
pca_df['Batch'] = pca_df.index.map(lambda p: (Path(p).parent.stem))
pca_df
# %%
hue_cols = ['Sex', 'Strain', 'Batch']

fig, axes = plt.subplots(1, 3, figsize=(18,6))

for ax, hue in zip(axes, hue_cols):    
    sns.scatterplot(x = 'PC1', y = 'PC2', hue = hue, data=pca_df, ax = ax)
    ax.set_xlabel(f"PC1 ({expl_var_dict['PC1']*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({expl_var_dict['PC2']*100:.1f}%)")

plt.show()

# %%
## Odd outlier, don't know why yet though
pca_df[pca_df['PC1'] > 200]

# %%
outlier = "videos/2025-08-20/100594_Female_Cdkl5_trimmed"
excluded_outlier_df = filtered_df.drop(index = outlier)
# %%
pca_df, expl_var_dict = run_pca_pipeline(excluded_outlier_df)
pca_df['Sex'] = pca_df.index.map(lambda p: (Path(p).name).split("_")[1].lower())
pca_df['Strain'] = pca_df.index.map(lambda p: (Path(p).name).split("_")[2])
pca_df['Batch'] = pca_df.index.map(lambda p: (Path(p).parent.stem))
pca_df
# %%
hue_cols = ['Sex', 'Strain', 'Batch']

fig, axes = plt.subplots(1, 3, figsize=(18,6))

for ax, hue in zip(axes, hue_cols):    
    sns.scatterplot(x = 'PC1', y = 'PC2', hue = hue, data=pca_df, ax = ax)
    ax.set_xlabel(f"PC1 ({expl_var_dict['PC1']*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({expl_var_dict['PC2']*100:.1f}%)")

plt.show()
# %%
half_done_strains = pca_df['Strain'].value_counts()[pca_df['Strain'].value_counts() >=9].index
half_done_strains = ['B6NJ', 'B6J', 'Shank3', 'Ube3a', 'Fmr1', 'Smarcc2']
half_done_strains_idx = pca_df[pca_df['Strain'].isin(half_done_strains)].index

within_half_done_df = filtered_df[filtered_df.index.isin(half_done_strains_idx)]

# %% PCA, just because
pca_df, expl_var_dict = run_pca_pipeline(within_half_done_df)
pca_df['Sex'] = pca_df.index.map(lambda p: (Path(p).name).split("_")[1].lower())
pca_df['Strain'] = pca_df.index.map(lambda p: (Path(p).name).split("_")[2])
pca_df['Batch'] = pca_df.index.map(lambda p: (Path(p).parent.stem))

hue_cols = ['Sex', 'Strain', 'Batch']

fig, axes = plt.subplots(1, 3, figsize=(18,6))

for ax, hue in zip(axes, hue_cols):    
    sns.scatterplot(x = 'PC1', y = 'PC2', hue = hue, data=pca_df, ax = ax)
    ax.set_xlabel(f"PC1 ({expl_var_dict['PC1']*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({expl_var_dict['PC2']*100:.1f}%)")

plt.show()
# %%
imputer = SimpleImputer(strategy="median")
df_imputed = imputer.fit_transform(within_half_done_df)
scaler = StandardScaler()
df_scaled = pd.DataFrame(
    scaler.fit_transform(df_imputed),
    index=within_half_done_df.index,
    columns=within_half_done_df.columns
)

# %%
df_scaled['Strain'] = df_scaled.index.map(lambda p: Path(p).name.split("_")[2])

strain_unique = df_scaled['Strain'].unique()
palette = sns.color_palette("tab20", n_colors=len(strain_unique))
strain_to_color = dict(zip(strain_unique, palette))

# Map row colors
row_colors = df_scaled['Strain'].map(strain_to_color)

# Drop strain from data before clustering
data_for_heatmap = df_scaled.drop(columns=['Strain'])

sns.clustermap(
    data_for_heatmap,
    cmap="vlag",
    center=0,
    figsize=(12,10),
    row_colors=row_colors,
    yticklabels=False,
)
import matplotlib.patches as mpatches
# Create legend handles
patches = [mpatches.Patch(color=color, label=strain) for strain, color in strain_to_color.items()]

# Place the legend (outside the heatmap)
plt.legend(handles=patches, bbox_to_anchor=(1, 1), title="Strain")
plt.show()

# %%
strain_to_color
# %%
features_corr_df = data_for_heatmap.corr()
plt.figure(figsize=(40, 35))
sns.heatmap(features_corr_df, center=0, cmap='vlag')
# %%
merged_df.shape
# %%
filtered_df.shape
# %%
na_cols = filtered_df.isna().sum().sort_values()
# %%
care_for = ['jumping_bout_behavior' in feat for feat in (na_cols.index)]
na_cols[care_for]