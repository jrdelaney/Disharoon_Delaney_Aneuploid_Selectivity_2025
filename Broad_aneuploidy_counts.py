import pandas as pd

# Read in files
df = pd.read_csv("Final_Publication_Files/Supplemental Data 2 Full Broad Dataset Expanded on MOA.csv", low_memory=False)

# Remove duplicate rows from data expansion onf MOA/Treatment
df.drop_duplicates(subset=['cell_line'], inplace=True)

# Melt the dataset to transform chromosome arms into rows for aggregation
chromosome_columns = [col for col in df.columns if col.startswith('chr')]
melted_df = df.melt(id_vars=['cancer_type'], value_vars=chromosome_columns, var_name='Chromosome Arm', value_name='State')

# Group by cancer type and chromosome arm, then count the occurrences of each state
result = melted_df.groupby(['cancer_type', 'Chromosome Arm', 'State']).size().reset_index(name='Count')

# Pivot to get counts for each state (-1, 0, 1) as columns
pivot_result = result.pivot(index=['cancer_type', 'Chromosome Arm'], columns='State', values='Count').reset_index()

# Rename the columns for clarity
pivot_result.columns.name = None  # Remove multi-level column name
pivot_result = pivot_result.rename(columns={-1: 'Losses', 0: 'Unchanged', 1: 'Gains'})

# Replace NaN values with 0 (for cases where no counts exist for a state)
pivot_result = pivot_result.fillna(0).astype({'Losses': int, 'Unchanged': int, 'Gains': int})

# Save the resulting table to a CSV file
pivot_result.to_csv('Supplemental Table 12 Broad Cancer Ploidy States.csv', index=False)
