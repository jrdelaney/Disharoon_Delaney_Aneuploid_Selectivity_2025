import pandas as pd

# Step 1: Create the mapping dictionary
mapping = {
    'Non-Small Cell Lung Cancer': ['LUNG'],
    'Endometrial Cancer': ['ENDOMETRIUM'],
    'Bladder Cancer': ['URINARY_TRACT'],
    'Seminoma': ['URINARY_TRACT'],
    'Non-Seminomatous Germ Cell Tumor': ['URINARY_TRACT'],
    'Esophagogastric Cancer': ['STOMACH', 'OESOPHAGUS'],
    'Pancreatic Cancer': ['PANCREAS'],
    'Renal Non-Clear Cell Carcinoma': ['KIDNEY'],
    'Sarcoma': ['SOFT_TISSUE'],
    'Breast Cancer': ['BREAST'],
    'Colorectal Cancer': ['LARGE_INTESTINE'],
    'Melanoma': ['SKIN'],
    'Renal Clear Cell Carcinoma': ['KIDNEY'],
    'Cervical Cancer': ['ENDOMETRIUM'],
    'Glioma': ['CENTRAL_NERVOUS_SYSTEM'],
    'Mature B-Cell Neoplasms': ['HAEMATOPOIETIC_AND_LYMPHOID_TISSUE']
}

# Step 2: Read both datasets into pandas DataFrames
TCGA_df = pd.read_csv('../../Final_Publication_Files/Supplemental Table 5 TCGA MOA Survival Prognosis Associations.csv')
Broad_df = pd.read_csv(
    '../../Final_Publication_Files/Supplemental Table 7 Broad MOA Viability Associations by Cancer Type.csv')

# Step 3: Map the cancer types in Dataset 1
TCGA_df['mapped_cancer_types'] = TCGA_df['cancer_type'].map(mapping)

# Step 4: Expand Dataset 1 to handle multiple mappings
TCGA_df_expanded = TCGA_df.explode('mapped_cancer_types').dropna(subset=['mapped_cancer_types'])
TCGA_df_expanded.rename(columns={'mapped_cancer_types': 'cancer_type'}, inplace=True)

# Ensure the 'treatment' and 'MOA' columns are in the same string format
TCGA_df_expanded['MOA'] = TCGA_df_expanded['MOA'].str.title()
Broad_df['MOA'] = Broad_df['MOA'].str.title()

# Step 5: Create sets of tuples for both datasets
# For Dataset 1
TCGA_set = set(TCGA_df_expanded[['MOA', 'arm', 'cancer_type']].apply(tuple, axis=1))

# For Dataset 2
Broad_set = set(Broad_df[['MOA', 'arm', 'cancer_type']].apply(tuple, axis=1))

# Step 6: Find common associations
common_associations = TCGA_set & Broad_set

# Step 7: Extract the common associations from both datasets
# Filter TCGA_df_expanded and Broad_df based on the common associations
TCGA_df_common = TCGA_df_expanded[TCGA_df_expanded[['MOA', 'arm', 'cancer_type']].apply(tuple, axis=1).isin(common_associations)]
Broad_df_common = Broad_df[Broad_df[['MOA', 'arm', 'cancer_type']].apply(tuple, axis=1).isin(common_associations)]

# Step 8: Merge the common associations dataframes on 'MOA', 'arm', and 'cancer_type'
merged_df = pd.merge(
    TCGA_df_common,
    Broad_df_common,
    on=['MOA', 'arm', 'cancer_type'],
    how='inner',
    suffixes=('_TCGA', '_Broad')
)

# Step 9: Select and rename the desired columns
# Columns from TCGA_df (Table 1) to include and prefix with 'TCGA_'
tcga_columns = ['log_rank_p_value', 'TCGA_median_loss_survival_months', 'TCGA_median_other_survival_months']
tcga_columns_prefixed = {col: f'TCGA_{col}' for col in tcga_columns}

# Columns from Broad_df (Table 2) to include and prefix with 'Broad_'
broad_columns = ['p_value', 'viability_shift']
broad_columns_prefixed = {col: f'Broad_{col}' for col in broad_columns}

# Rename the columns in the merged dataframe
merged_df.rename(columns={**tcga_columns_prefixed, **broad_columns_prefixed}, inplace=True)

# Step 10: Select the final columns to include in the output
final_columns = ['MOA', 'arm', 'cancer_type'] + list(tcga_columns_prefixed.values()) + list(broad_columns_prefixed.values())
final_df = merged_df[final_columns]
final_df['cancer_type'] = final_df['cancer_type'].str.title()

# Step 11: Save the final dataframe to a CSV file
final_df.to_csv('Final_Publication_Files/Supplemental Table 8 Common Effective Between TCGA and Broad MOAs.csv', index=False)