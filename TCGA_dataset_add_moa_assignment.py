import os
import json
import pandas as pd
import re
import ast
from collections import Counter

# Directory containing the JSONL files
directory = "gtp_out_tcga_broad"
results_list = []

# Helper function to fix apostrophes inside JSON-like structures
def fix_apostrophes(content):
    # Escape only apostrophes that are inside words like BRUTON'S and leave others intact
    return re.sub(r"(?<=\w)'(?=\w)", r"\\'", content)


# Iterate over all files in the directory
for filename in os.listdir(directory):
    if filename.endswith(".jsonl"):
        file_path = os.path.join(directory, filename)
        with open(file_path, 'r') as file:
            # Read each line (JSON object) from the file
            for line in file:
                data = json.loads(line.strip())
                # print(data)
                response_body = data['response']['body']
                choices = response_body.get('choices', [])
                if choices:
                    message_content = choices[0]['message']['content']
                    try:
                        # Fix apostrophe issues before parsing
                        fixed_content = fix_apostrophes(message_content)

                        # Safely parse the fixed message content into a Python list of dictionaries
                        response_list = ast.literal_eval(fixed_content.strip())
                        if isinstance(response_list, list):
                            for item in response_list:
                                # print(item)
                                results_list.append(item)
                    except (SyntaxError, ValueError) as e:
                        print(f"Skipping malformed line in file {filename}: {e}")

# Convert the results list to a DataFrame
df = pd.DataFrame(results_list)

# Now calculate the mode of 'match' values for each therapy/MOA pair
def calculate_mode(group):
    # Find the mode of the 'match' values
    match_mode = Counter(group['match']).most_common(1)[0][0]
    # Filter rows where the 'match' value equals the mode
    return group[group['match'] == match_mode]

# Group by Therapy and MOA, and apply the mode filter function
df_filtered = df.groupby(['therapy', 'MOA'], group_keys=False).apply(calculate_mode)
df_filtered = df_filtered.drop_duplicates()
pivot_df = df_filtered.pivot(index='therapy', columns='MOA', values='match')
# Optionally, fill any missing values with 'False' (if no match was found for a given MOA-therapy pair)
pivot_df = pivot_df.fillna(False)
pivot_df.to_csv('therapy_moa_match_results_tcga_broad.csv',index=True)
print(f"Results saved")

# Read the existing data
df = pd.read_csv('final_TCGA_datasets/7_TCGA_final_exploded_treatment_aneuploidy_by_treatment.csv')

# Ensure that 'treatment' and 'therapy' are in the same case for matching
df['treatment'] = df['treatment'].str.lower()
pivot_df.index = pivot_df.index.astype(str).str.lower()
therapy_moas = pivot_df.apply(lambda row: row[row == True].index.tolist(), axis=1).to_dict()

# Map the treatments in df to their MOAs
df['moas'] = df['treatment'].map(therapy_moas)

# Explode the DataFrame on the 'moas' column
df_exploded = df.explode('moas')
df_exploded.drop_duplicates(inplace=True)
# Save the final table
df_exploded.to_csv('final_TCGA_datasets/8_TCGA_final_exploded_MOA_exploded_treatment_aneuploidy_by_treatment.csv', index=False)

