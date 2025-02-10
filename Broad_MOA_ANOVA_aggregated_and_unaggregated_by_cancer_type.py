import os
import pandas as pd
import scipy.stats as stats
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from statsmodels.stats.multitest import multipletests  # <-- For BH correction

# ---------------------------
# Configuration and Setup
# ---------------------------

SAVE_DIR = '../Final_analyses'

# Create necessary directories
os.makedirs(SAVE_DIR, exist_ok=True)

# Define chromosome arms
CHROMOSOME_ARMS = [
    '1p', '1q', '2p', '2q', '3p', '3q', '4p', '4q', '5p', '5q',
    '6p', '6q', '7p', '7q', '8p', '8q', '9p', '9q', '10p', '10q',
    '11p', '11q', '12p', '12q', '13q', '14q', '15q', '16p', '16q',
    '17p', '17q', '18p', '18q', '19p', '19q', '20p', '20q', '21q', '22q'
]

cutoff = 10  # Minimum 10 observations per group


# ---------------------------
# Utility Functions
# ---------------------------
def perform_anova(filtered_data, chromosome_arm):
    """Perform a one-way ANOVA on the groups split by chromosome_arm."""
    groups = filtered_data.groupby(chromosome_arm)["cell_viability"].apply(list)
    F_statistic, p_value = stats.f_oneway(*groups)

    # Clamp very small p-values to avoid zero or underflow in BH correction
    if p_value < 1e-300:
        p_value = 1e-300

    return F_statistic, p_value


# ---------------------------
# Data Loading and Preprocessing
# ---------------------------
def load_and_preprocess_data(data_filepath, treatment_col):
    # Load file
    df = pd.read_csv(data_filepath)

    if treatment_col == 'drug_MOA':
        # Return the dataframe as-is if the treatment column is 'drug_MOA'
        return df

    # Unexplode the exploded data for MOA (remove duplicate rows for the same treatment)
    df = df.drop_duplicates(subset=['cell_line', 'drug_name'], keep='first')
    df = df.reset_index(drop=True)
    return df


# ---------------------------
# Analysis Functions
# ---------------------------
def run_stat_analysis(df, treatment, treatment_col, chromosome_arm, cancer_type=None):
    # Filter data for the selected drug and chromosome arm
    filtered_data = df[(df[treatment_col] == treatment) & df[chromosome_arm].notna()]
    filtered_data = filtered_data[[chromosome_arm, "cell_viability"]].dropna()

    # Replace numerical values with categorical labels
    filtered_data[chromosome_arm] = filtered_data[chromosome_arm].replace({-1: "Loss", 0: "Other", 1: "Other"})

    if filtered_data.empty:
        return None

    group_loss = filtered_data[filtered_data[chromosome_arm] == 'Loss']
    group_other = filtered_data[filtered_data[chromosome_arm] == 'Other']

    # Check cutoff
    if len(group_other) <= cutoff or len(group_loss) <= cutoff:
        return None

    # Perform ANOVA
    F_statistic, p_value = perform_anova(filtered_data, chromosome_arm)

    # Calculate means
    means = filtered_data.groupby(chromosome_arm)["cell_viability"].mean()
    loss_mean = means.get("Loss")
    other_mean = means.get("Other")

    # Determine shift direction
    if loss_mean is not None and other_mean is not None:
        if loss_mean > other_mean:
            shift = "Greater viability in loss lines"
        else:
            shift = "Lower viability in loss lines"
    else:
        shift = "Indeterminable viability difference"

    # Return result dictionary with group sizes included
    return {
        'Treatment': treatment,
        'chr_arm': chromosome_arm,
        'cancer_type': cancer_type if cancer_type else "All",
        'p_value': p_value,
        'sensitivity_change': shift,
        'loss_size': len(group_loss),
        'other_size': len(group_other)
    }


def analyze_by_cancer_type(cancer_type, df, chromosome_arms, treatment_col):
    # Analyze treatment effectiveness by a specific cancer type
    results = []
    print(f"Processing cancer type: {cancer_type}")

    # Subset data for the current cancer type
    filtered_data = df[df['cancer_type'] == cancer_type]

    # Identify the columns to check for duplicates
    columns_to_check = filtered_data.columns.difference([treatment_col])
    # Drop duplicates based on the selected columns
    cancer_data = filtered_data.drop_duplicates(subset=columns_to_check)

    # Loop through data by arm and check if there's data and run analysis
    for arm in chromosome_arms:
        chrom_arm_col = f"chr{arm}"
        if chrom_arm_col not in cancer_data.columns:
            print(f"Chromosome arm column {chrom_arm_col} not found in cancer type {cancer_type}. Skipping.")
            continue

        # Using threshold=0 here so that we include all results
        threshold = 0
        treatments_for_arm = cancer_data[cancer_data[chrom_arm_col] >= threshold][treatment_col].unique()

        for treatment in treatments_for_arm:
            result = run_stat_analysis(
                df=cancer_data,
                treatment=treatment,
                treatment_col=treatment_col,
                chromosome_arm=chrom_arm_col,
                cancer_type=cancer_type,
            )
            if result:
                print(
                    f"Generated results for Treatment: {treatment}, Chromosome Arm: {chrom_arm_col}, Cancer Type: {cancer_type}")
                results.append(result)

    return results


def analyze_aggregate(df, chromosome_arms, treatment_col):
    # Analyze all treatments across all cancer types (aggregate)
    results = []
    print("Processing aggregate dataset (All cancer types combined).")

    # Drop duplicates similarly, but no filter on cancer_type
    columns_to_check = df.columns.difference([treatment_col])
    all_data = df.drop_duplicates(subset=columns_to_check)

    for arm in chromosome_arms:
        chrom_arm_col = f"chr{arm}"
        if chrom_arm_col not in all_data.columns:
            print(f"Chromosome arm column {chrom_arm_col} not found. Skipping.")
            continue

        # Include all possible treatments for rows with >= 0
        threshold = 0
        treatments_for_arm = all_data[all_data[chrom_arm_col] >= threshold][treatment_col].unique()

        for treatment in treatments_for_arm:
            result = run_stat_analysis(
                df=all_data,
                treatment=treatment,
                treatment_col=treatment_col,
                chromosome_arm=chrom_arm_col,
                cancer_type=None  # We'll label these rows as "All"
            )
            if result:
                print(
                    f"Generated results for Treatment: {treatment}, Chromosome Arm: {chrom_arm_col}, Cancer Type: ALL")
                results.append(result)

    return results


# ---------------------------
# Main Execution
# ---------------------------
def main():
    # Load and preprocess data
    data_filepath = "../../Final_Publication_Files/Supplemental Data 2 Full Broad Dataset Expanded on MOA.csv"
    treatment_col = 'drug_MOA'
    dataset = load_and_preprocess_data(data_filepath, treatment_col)

    # Initialize results list
    all_results = []

    # Get unique cancer types
    cancer_types = dataset['cancer_type'].dropna().unique()
    print(f"Found {len(cancer_types)} unique cancer types.")

    # Use ProcessPoolExecutor for parallel processing
    with ProcessPoolExecutor() as executor:
        # Partial function with fixed dataframe and chromosome arms
        func = partial(analyze_by_cancer_type, df=dataset, chromosome_arms=CHROMOSOME_ARMS, treatment_col=treatment_col)

        # Submit all cancer type analyses to the executor
        futures = {executor.submit(func, cancer_type): cancer_type for cancer_type in cancer_types}

        for future in as_completed(futures):
            cancer_type = futures[future]
            try:
                result = future.result()
                all_results.extend(result)
            except Exception as exc:
                print(f"Cancer type {cancer_type} generated an exception: {exc}")

    # Convert results to DataFrame
    results_df = pd.DataFrame(all_results)

    # Perform BH multiple-testing correction for cancer-type-specific results
    if not results_df.empty:
        pvals = results_df['p_value'].values
        reject, pvals_corrected, _, _ = multipletests(pvals, method='fdr_bh')
        results_df['p_value_BH_adjusted'] = pvals_corrected

    # Save the cancer-type-specific results
    results_csv_path = os.path.join(SAVE_DIR, 'Final_Broad_MOA_V_Cancer_Type_V_Aneuploidy_Results.csv')
    results_df.to_csv(results_csv_path, index=False)
    print(f"Analysis with aggregation by cancer type completed. Results saved to '{results_csv_path}'.")

    # Now, run the aggregate analysis (all cancers combined)
    aggregate_results = analyze_aggregate(df=dataset, chromosome_arms=CHROMOSOME_ARMS, treatment_col=treatment_col)
    agg_results_df = pd.DataFrame(aggregate_results)

    # Perform BH multiple-testing correction on the aggregate results
    if not agg_results_df.empty:
        pvals = agg_results_df['p_value'].values
        reject, pvals_corrected, _, _ = multipletests(pvals, method='fdr_bh')
        agg_results_df['p_value_BH_adjusted'] = pvals_corrected

    # Save the aggregate results
    agg_results_csv_path = os.path.join(SAVE_DIR, 'Final_Broad_MOA_V_ALL_Cancers_Aneuploidy_Results.csv')
    agg_results_df.to_csv(agg_results_csv_path, index=False)
    print(f"Aggregate analysis (All cancers combined) completed. Results saved to '{agg_results_csv_path}'.")


if __name__ == "__main__":
    main()
