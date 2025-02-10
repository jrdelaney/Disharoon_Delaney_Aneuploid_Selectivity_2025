import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import os
import numpy as np

def analyze_survival_curves(
    df_filtered,
    cancer_type_filter,
    drug_name_filter,
    chromosome_arm,
    drug_name_clean,
    cancer_name_clean,
    therapy_col,
    survival_time_metric,
    save_dir,
    result_list
):
    min_observations = 10
    df_tmp1 = df_filtered[df_filtered['Cancer Type'] == cancer_type_filter].copy()
    # Iterate over each treatment
    df_tmp2 = df_tmp1[df_tmp1[therapy_col] == drug_name_filter].copy()
    # Check if there are enough samples for the treatment
    if df_tmp2.empty:
        print(f"No data found for {drug_name_filter} in {cancer_type_filter}")
        return None
    df_tmp3 = df_tmp2.dropna(subset=[chromosome_arm]).copy()
    if df_tmp3.empty:
        print(f"No data found for {chromosome_arm} in {drug_name_filter} in {cancer_type_filter}")
        return None

    group_loss = df_tmp3[df_tmp3[f'{chromosome_arm}_group'] == 'Loss']
    group_other = df_tmp3[df_tmp3[f'{chromosome_arm}_group'] == 'Other']

    if len(group_loss) < min_observations:
        print(f"Not enough loss observations.")
        return None
    if len(group_other) < min_observations:
        print(f"Not enough other observations.")
        return None

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

    # Handle cases where median survival time is not reached
    if pd.isna(median_loss) or pd.isna(median_other):
        print("Median can't be calculated")
        return None  # Skip if median survival is not reached in either group

    # Plot survival curves
    plt.figure(figsize=(5, 4))
    ax = plt.subplot(111)

    kmf_loss.plot_survival_function(ax=ax, ci_show=True)
    kmf_other.plot_survival_function(ax=ax, ci_show=True)

    plt.title(f'{drug_name_clean} Treatment of {chromosome_arm} Aneuploidy\n{cancer_name_clean} (p={log_rank.p_value:.4f})')
    plt.xlabel('Months')
    plt.ylabel('Survival Probability')
    plt.xlim(0, 60)
    plt.xticks(np.arange(0, 61, 12))  # Labeled ticks every 12 months
    plt.yticks(np.linspace(0, 1, 3))
    plt.legend(framealpha=1)
    plt.grid(True, linewidth=0.5)

    # Save plot
    fig_name = f'survival_curve_{chromosome_arm}_{cancer_name_clean}_{drug_name_clean}.png'
    fig_path = os.path.join(save_dir, fig_name)
    plt.savefig(fig_path,dpi=600, bbox_inches='tight')
    plt.close()

    # Append result to list
    result = {
        'treatment': drug_name_clean,
        'cancer': cancer_name_clean,
        'aneuploidy': chromosome_arm,
        'log_rank_p_value': f'{log_rank.p_value:9.4f}',
        'median_loss': median_loss,
        'median_other': median_other,
        'n_loss': len(group_loss),
        'n_other': len(group_other),
        'graph_name': fig_name
    }
    result_list.append(result)

analyses = [
    {
        'cancer_type_filter': 'Head and Neck Cancer',
        'drug_name_filter': 'carboplatin',
        'chromosome_arm': '9p',
        'drug_name_clean': 'Carboplatin',
        'cancer_name_clean': 'Head And Neck'
    },
    {
        'cancer_type_filter': 'Head and Neck Cancer',
        'drug_name_filter': 'paclitaxel',
        'chromosome_arm': '21q',
        'drug_name_clean': 'Paclitaxel',
        'cancer_name_clean': 'Head And Neck'
    },
    {
        "cancer_type_filter": "Non-Small Cell Lung Cancer",
        "drug_name_filter": "cisplatin",
        "chromosome_arm": "17p",
        "drug_name_clean": "Cisplatin",
        "cancer_name_clean": "Non-Small Cell Lung"
    },
    {
        "cancer_type_filter": "Non-Small Cell Lung Cancer",
        "drug_name_filter": "cisplatin",
        "chromosome_arm": "8p",
        "drug_name_clean": "Cisplatin",
        "cancer_name_clean": "Non-Small Cell Lung"
    },
    {
        "cancer_type_filter": "Non-Small Cell Lung Cancer",
        "drug_name_filter": "alimta",
        "chromosome_arm": "17p",
        "drug_name_clean": "Pemetrexed",
        "cancer_name_clean": "Non-Small Cell Lung"
    },
    {
        "cancer_type_filter": "Non-Small Cell Lung Cancer",
        "drug_name_filter": "docetaxel",
        "chromosome_arm": "5q",
        "drug_name_clean": "Docetaxel",
        "cancer_name_clean": "Non-Small Cell Lung"
    },
    {
        "cancer_type_filter": "Non-Small Cell Lung Cancer",
        "drug_name_filter": "vincristine",
        "chromosome_arm": "3p",
        "drug_name_clean": "Vincristine",
        "cancer_name_clean": "Non-Small Cell Lung"
    }
]

save_dir = 'FIGURES_all_TCGA_treatments_survival_curves'
os.makedirs(save_dir, exist_ok=True)

# Load TCGA data
df_exploded = pd.read_csv("../Final_Publication_Files/Supplemental Data 1 Full TCGA Dataset Expanded on MOA.csv", low_memory=False)

# Unexplode the 'moas' column by dropping and dropping duplicates
df_no_moas = df_exploded.drop(columns=['moas'])
df = df_no_moas.drop_duplicates(subset=['Sample ID', 'treatment',  "Study ID"])

therapy_col = 'treatment' # could use therapy if looking at specific drugs
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
min_observations = 10

for analysis in analyses:
    analyze_survival_curves(
        df_filtered = df_filtered,
        cancer_type_filter=analysis['cancer_type_filter'],
        drug_name_filter=analysis['drug_name_filter'],
        chromosome_arm=analysis['chromosome_arm'],
        drug_name_clean=analysis['drug_name_clean'],
        cancer_name_clean=analysis['cancer_name_clean'],
        therapy_col = therapy_col,
        survival_time_metric = survival_time_metric,
        save_dir = save_dir,
        result_list = result_list,
    )

# Save results to CSV
result_df = pd.DataFrame(result_list)
result_df.to_csv(os.path.join(save_dir, f'{save_dir}_survival_curve_results.csv'), index=False)
