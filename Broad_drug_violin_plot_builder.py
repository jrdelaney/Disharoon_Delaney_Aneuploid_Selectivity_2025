import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
from scipy import stats


def strip_before_char(string, char):
    return string.split(char, 1)[-1]

def load_and_preprocess_data(data_filepath, treatment_col):
    df = pd.read_csv(data_filepath)
    if treatment_col == 'drug_MOA':
        # Return the dataframe as is if the treatment column is 'drug_MOA'
        return df
    # Unexplode the exploded data for MOA (remove duplicate rows for the same treatment)
    # Drop any duplicate rows that share the same (cell_line, drug_name)
    df = df.drop_duplicates(subset=['cell_line', 'drug_name'], keep='first')
    # Reset the index to return a flat DataFrame
    df = df.reset_index(drop=True)
    return df

def analyze_aneuploidy_sensitivity(
        data_filepath,
        treatment_col,
        cancer_type_filter,
        treatment_filter,
        chromosome_arm,
        cancer_name_clean
):
    dataset = load_and_preprocess_data(data_filepath, treatment_col)
    
    # Filter data for the specific cancer type and treatment (drug, MOA)
    print(f"Filtering data for cancer type: {cancer_type_filter} and treatment: {treatment_filter}...")
    filtered_data = dataset[
        (dataset['cancer_type'] == cancer_type_filter)
    ]
    print(f"Number of records after cancer type filter: {len(filtered_data)}")

    filtered_data = filtered_data[
        (filtered_data[treatment_col].str.title() == treatment_filter)
    ]
    print(f"Number of records after treatment filter: {len(filtered_data)}\n")

    # Select chromosome arm
    print(f"Selecting chromosome arm: {chromosome_arm}...")
    # Filter data for chromosome_arm
    filtered_data = filtered_data[filtered_data[chromosome_arm].notna()]
    print(f"Number of records after chromosome arm notna filter: {len(filtered_data)}")

    # Replace -1 with 'Loss', 0 and 1 with 'Other' for ploidy groups
    print(f"Replacing chromosome arm values for: {chromosome_arm}...")
    filtered_data[chromosome_arm] = filtered_data[chromosome_arm].replace({-1: 'Loss', 0: 'Other', 1: 'Other'})
    print(f"Value counts after replacement:\n{filtered_data[chromosome_arm].value_counts()}\n")

    # Ensure sufficient data is available
    print("Ensuring sufficient data for ANOVA...")
    group_loss = filtered_data[filtered_data[chromosome_arm] == 'Loss']
    group_other = filtered_data[filtered_data[chromosome_arm] == 'Other']
    print(f"Group 'Loss' count: {len(group_loss)}")
    print(f"Group 'Other' count: {len(group_other)}\n")

    # Check for NaN values in 'cell_viability'
    num_nan_z = filtered_data['cell_viability'].isna().sum()
    print(f"Number of NaN values in 'cell_viability': {num_nan_z}")

    print("Removing records with NaN 'cell_viability' values...\n")
    filtered_data = filtered_data.dropna(subset=["cell_viability"])
    print(f"Number of records after dropping NaNs: {len(filtered_data)}\n")

    if len(group_loss) > 10 and len(group_other) > 10:
        # Perform ANOVA
        print("Performing ANOVA...")
        groups = filtered_data.groupby(chromosome_arm)['cell_viability'].apply(list)
        F_statistic, p_value = stats.f_oneway(*groups)
        print(f"ANOVA results - F_statistic: {F_statistic}, p_value: {p_value}\n")

        if p_value < 0.05:
            # Prepare plot
            print("ANOVA significant. Preparing violin plot...")
            fig, ax = plt.subplots(figsize=(4, 3))
            graph_colors = {"Other": "#FBCEB1", "Loss": "#ADD8E6"}
            sns.violinplot(
                x=chromosome_arm,
                y="cell_viability",
                data=filtered_data,
                palette=graph_colors,
                order=['Loss', 'Other'],
                split=True,
                gap=-0.2,
                inner_kws=dict(box_width=8, whis_width=2, color="0.5"),
                ax=ax
            )

            # Clean drug name
            clean_drug = re.sub(r'[\\/*?:"<>|]', '_', treatment_filter)
            clean_drug = clean_drug.title()

            # Clean cancer type (remove underscores and format)
            clean_cancer = cancer_name_clean.title()

            # Modify title to add a new line after the MOA name
            ax.set_title(
                f"{clean_drug} Sensitivity in\n {clean_cancer} Cancer {chromosome_arm.split('chr')[1]} Loss (p={p_value:.4f})",
                fontsize=8
            )
            ax.set_ylabel("Viability", fontsize=8)
            ax.set_xlabel("", fontsize=0)
            ax.tick_params(axis='y', which='major', labelsize=6)
            ax.tick_params(axis='x', labelsize=8)

            # Save plot
            chart_name = f'violin_plot_{chromosome_arm}_{clean_drug}_{clean_cancer}_cancer.png'
            plt.savefig(chart_name, dpi=600)
            plt.close()
            print(f'Plot saved as {chart_name}\n')
        else:
            print('ANOVA not significant (p >= 0.05). No plot generated.\n')
    else:
        print('Not enough data to generate plot.\n')


# Parameters for each analysis
analyses = [
    {
        'cancer_type_filter': 'OVARY',
        'treatment_filter': 'Pf-06651600',
        'chromosome_arm': 'chr4p',
        'drug_name_clean': 'PF-06651600',
        'cancer_name_clean': 'Ovary'
    },
    {
        'cancer_type_filter': 'OVARY',
        'treatment_filter': 'Pf-06651600',
        'chromosome_arm': 'chr17p',
        'drug_name_clean': 'PF-06651600',
        'cancer_name_clean': 'Ovary'
    },
    {
        'cancer_type_filter': 'OVARY',
        'treatment_filter': 'Etomidate',
        'chromosome_arm': 'chr4p',
        'drug_name_clean': 'Etomidate',
        'cancer_name_clean': 'Ovary'
    },
    {
        'cancer_type_filter': 'OVARY',
        'treatment_filter': 'Etomidate',
        'chromosome_arm': 'chr17p',
        'drug_name_clean': 'Etomidate',
        'cancer_name_clean': 'Ovary'
    },
    {
        'cancer_type_filter': 'PANCREAS',
        'treatment_filter': 'Mk-8745',
        'chromosome_arm': 'chr3p',
        'drug_name_clean': 'MK-8745',
        'cancer_name_clean': 'Pancreas'
    },
    {
        'cancer_type_filter': 'PANCREAS',
        'treatment_filter': 'Losmapimod',
        'chromosome_arm': 'chr17p',
        'drug_name_clean': 'Losmapimod',
        'cancer_name_clean': 'Pancreas'
    },
]

# Execute analyses
for analysis in analyses:
    analyze_aneuploidy_sensitivity(
        data_filepath="../Final_Publication_Files/Supplemental Data 2 Full Broad Dataset Expanded on MOA.csv",
        treatment_col='drug_name',
        cancer_type_filter=analysis['cancer_type_filter'],
        treatment_filter=analysis['treatment_filter'],
        chromosome_arm=analysis['chromosome_arm'],
        # drug_name_clean=analysis['drug_name_clean'], # Depricated
        cancer_name_clean=analysis['cancer_name_clean']
    )
