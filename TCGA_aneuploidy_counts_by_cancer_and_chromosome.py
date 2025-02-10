import pandas as pd
import numpy as np

df_exploded = pd.read_csv("../Final_Publication_Files/Supplemental Data 1 Full TCGA Dataset Expanded on MOA.csv", low_memory=False)

# Unexplode the 'moas' column by dropping and dropping duplicates
df_no_moas = df_exploded.drop(columns=['moas'])
df = df_no_moas.drop_duplicates(subset=['Sample ID', 'treatment',  "Study ID"])


arms = ['1p', '1q', '2p', '2q', '3p', '3q', '4p', '4q', '5p', '5q',
        '6p', '6q', '7p', '7q', '8p', '8q', '9p', '9q', '10p', '10q',
        '11p', '11q', '12p', '12q', '13q', '14q', '15q', '16p', '16q',
        '17p', '17q', '18p', '18q', '19p', '19q', '20p', '20q', '21q', '22q']

result_list = []
for ct in df['Cancer Type'].unique():
    df_cancer = df[df['Cancer Type'] == ct]
    df_cancer = df_cancer.drop_duplicates(subset=['Sample ID']) # drop exploded duplicated rows
    for arm in arms:
        gain_count = np.where(df_cancer[arm] == "Gain", 1,0).sum()
        loss_count = np.where(df_cancer[arm] == "Loss",1,0).sum()
        unc_count = np.where(df_cancer[arm] == "Unchanged",1,0).sum()
        na_count = np.where(df_cancer[arm] == pd.NA,1,0).sum()
        result = {'Cancer Type': ct, "Chromosome Arm":arm, "Gains":gain_count, "Losses":loss_count, "Unchangeds":unc_count, "NAs":na_count}
        result_list.append(result)

df_final = pd.DataFrame(result_list)
df_final.to_csv('TCGA_anueplopidy_counts_by_cancer_type.csv', index=False)
# This file is used to make "Supplemental Table 10 TCGA Cancer Ploidy States with Normalized Counts.csv" after some changes to header names

# Normalize counts for each cancer type by its total sample count
normalized_columns = ["Gains", "Losses", "Unchangeds", "NAs"]
for col in normalized_columns:
    df_final[f"{col} (Normalized)"] = df_final[col] / df_final["Number of samples"]

# Save the normalized dataset for downstream analysis and graphing
df_final.to_csv('TCGA_normalized_aneuploidy_counts_by_cancer_type.csv', index=False)
# This file is used to make "Supplemental Table 10 TCGA Cancer Ploidy States with Normalized Counts.csv" after some changes to header names
