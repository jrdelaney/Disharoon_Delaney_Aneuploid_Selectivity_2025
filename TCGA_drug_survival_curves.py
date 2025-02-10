import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from statsmodels.stats.multitest import multipletests


# Load TCGA data
df_exploded = pd.read_csv("../../Final_Publication_Files/Supplemental Data 1 Full TCGA Dataset Expanded on MOA.csv", low_memory=False)

# NOTE: Change this to run for either moa or drug
therapy_col = 'treatment'  # 'treatment' or 'moas' NOTE: 'treatment' is a typically a compound (drug) but not always

if therapy_col == 'treatment':
    # Unexplode the 'moas' column by dropping and dropping duplicates
    df_no_moas = df_exploded.drop(columns=['moas'])
    df = df_no_moas.drop_duplicates(subset=['Sample ID', 'treatment', "Study ID"])
else:
    df = df_exploded.copy()

drug_list = df[df[therapy_col] != 'radiation'][therapy_col].unique()
drug_list = df[df[therapy_col] != 'external'][therapy_col].unique()

# Define survival time metric
survival_time_metric = 'Progress Free Survival (Months)'

# Drop missing rows for survival
df_filtered = df[df[survival_time_metric].notna()].copy()

# Set chromosomes of interest
chrs_interest = [
    '1p', '1q', '2p', '2q', '3p', '3q', '4p', '4q', '5p', '5q',
    '6p', '6q', '7p', '7q', '8p', '8q', '9p', '9q', '10p', '10q',
    '11p', '11q', '12p', '12q', '13q', '14q', '15q', '16p', '16q',
    '17p', '17q', '18p', '18q', '19p', '19q', '20p', '20q', '21q', '22q'
]

# Create a binary group column for each chromosome
for chr in chrs_interest:
    df_filtered[f'{chr}_group'] = df_filtered[chr].apply(lambda x: 'Loss' if x == 'Loss' else 'Other')

# Refactor outcomes to simple terms
df_filtered = df_filtered.dropna(subset=['Disease-specific Survival status']).copy()
df_filtered['Survival_Status'] = df_filtered['Disease-specific Survival status'].apply(
    lambda x: 1 if '0:ALIVE OR DEAD TUMOR FREE' in x else 0
)

result_list = []
min_observations = 10 # Minimum number of samples per ploidy group in order to run stat test

# Iterate over each cancer type
for cancer in df_filtered['Cancer Type'].unique():
    df_tmp1 = df_filtered[df_filtered['Cancer Type'] == cancer].copy()
    print(f"\n[INFO] Analyzing Cancer Type: {cancer}, Count: {len(df_tmp1)}")

    # Iterate over each treatment
    for sig_therapy in drug_list:
        df_tmp2 = df_tmp1[df_tmp1[therapy_col] == sig_therapy].copy()
        if df_tmp2.empty:
            print(f"  [SKIP] No samples for treatment: {sig_therapy}")
            continue
        print(f"  [INFO]  Treatment: {sig_therapy}, Count: {len(df_tmp2)}")

        # Iterate over each chromosome of interest
        for chr in chrs_interest:
            df_tmp3 = df_tmp2.dropna(subset=[chr]).copy()
            if df_tmp3.empty:
                print(f"    [SKIP] {chr} has all missing data for this treatment.")
                continue

            group_loss = df_tmp3[df_tmp3[f'{chr}_group'] == 'Loss']
            group_other = df_tmp3[df_tmp3[f'{chr}_group'] == 'Other']

            if len(group_loss) < min_observations:
                print(f"    [SKIP] Not enough 'Loss' samples for {chr} (need >= {min_observations}, have {len(group_loss)})")
                continue
            if len(group_other) < min_observations:
                print(f"    [SKIP] Not enough 'Other' samples for {chr} (need >= {min_observations}, have {len(group_other)})")
                continue

            # Perform log-rank test
            log_rank = logrank_test(
                group_loss[survival_time_metric],
                group_other[survival_time_metric],
                event_observed_A=group_loss['Survival_Status'],
                event_observed_B=group_other['Survival_Status']
            )

            # Initialize Kaplan-Meier fitter instances
            kmf_loss = KaplanMeierFitter()
            kmf_other = KaplanMeierFitter()

            # Fit the models
            kmf_loss.fit(
                group_loss[survival_time_metric],
                event_observed=group_loss['Survival_Status'],
                label='Loss'
            )
            kmf_other.fit(
                group_other[survival_time_metric],
                event_observed=group_other['Survival_Status'],
                label='Other'
            )

            # Determine median survival times
            median_loss = kmf_loss.median_survival_time_
            median_other = kmf_other.median_survival_time_

            # Handle cases where median survival is na by skipping analysis
            if pd.isna(median_loss) or pd.isna(median_other):
                print(f"    [SKIP] Median survival time not reached for {chr}")
                continue

            # Append result to list
            result = {
                'treatment': sig_therapy.title(),
                'cancer': cancer.title(),
                'aneuploidy': chr,
                'log_rank_p_value': f'{log_rank.p_value:9.4f}',
                'median_loss': median_loss,
                'median_other': median_other,
                'n_loss': len(group_loss),
                'n_other': len(group_other)
            }
            result_list.append(result)

# Save results to CSV
result_df = pd.DataFrame(result_list)

if not result_df.empty:
    # Extract p-values as a float Series
    pvals = result_df['log_rank_p_value'].astype(float).values

    # Apply Benjamini–Hochberg (FDR) correction
    reject, pvals_corrected, _, _ = multipletests(pvals, method='fdr_bh')

    # Append corrected p-values to result_df
    result_df['log_rank_p_adj'] = pvals_corrected
# ---------------------------

# Choose file name depending on therapy_col
therapy_file_name = 'drug' if therapy_col == 'treatment' else 'MOA'

# Save results to CSV
csv_path = f"{therapy_file_name}_survival_curve_results.csv"
result_df.to_csv(csv_path, index=False)
print(f"\n[INFO] Results saved to {csv_path}")
