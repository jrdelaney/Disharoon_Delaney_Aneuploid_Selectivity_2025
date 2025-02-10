import pandas as pd

df = pd.read_csv(
    '../../Final_Publication_Files/Supplemental Table 7 Broad MOA Viability Associations by Cancer Type.csv', low_memory=False)

# Count the occurrences of each MOA per 'arm' and 'viability_shift'
counts = df.groupby(['arm', 'viability_shift', 'MOA']).size().reset_index(name='count')

# Function to get top MOAs per group
def get_topk(group):
    max_count = group['count'].max()
    return group[group['count'] == max_count]

# Group by arm and viability, then get the shift
top_moas = counts.groupby(['arm', 'viability_shift']).apply(get_topk).reset_index(drop=True)

# Separate the results into sensitizing and desensitizing MOAs
top_sensitizing = top_moas[top_moas['viability_shift'] == 'Lower viability in loss lines']
top_desensitizing = top_moas[top_moas['viability_shift'] == 'Greater viability in loss lines']

# Combine the two tables and save the results
comb_sense = pd.concat([top_sensitizing,top_desensitizing])
comb_sense.to_csv('Supplemental Table 11 Broad Top Cytotoxic MOAs.csv')




