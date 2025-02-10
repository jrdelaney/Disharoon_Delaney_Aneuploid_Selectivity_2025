Description of scripts and data utlized for Disharoon and Delaney, 2025.

The purpose of these scripts are to compute associations of aneuploidy data from the DepMap project and TCGA project with drug response.

Contributors:
Andrew Disharoon, all Python scripts
Joe Delaney, R script

Python scripts:
TCGA_dataset_create_full.py
Automates downloading, combining, and cleaning TCGA and cBioPortal datasets. Processes arm-level copy-number alteration data and integrates it with clinical, cancer type, and treatment data.
Tasks:
●	Create a directory for output files.
●	Combine arm-level copy-number data from multiple GISTIC files into a single dataset.
●	Merge the combined dataset with cancer type information.
●	Rename long column names for clarity.
●	Integrate clinical data from TCGA studies.
●	Retrieve treatment data for specified studies using the cBioPortal API.
●	Save treatment data as JSON and CSV.
●	Merge treatment data with the combined aneuploidy-clinical dataset.
●	Clean treatment data by normalizing names and handling multiple treatments.
●	Explode treatment lists into individual rows for analysis.
●	Save the final dataset for downstream analysis, including exploded treatment rows.
TCGA_dataset_add_moa_assignment.py
Processes JSONL files to extract therapy and mechanism of action (MOA) matches, calculates the mode for each therapy-MOA pair, and integrates the results with existing treatment and chromosome aneuploidy data.
Tasks:
●	Load JSONL files from a specified directory and extract therapy-MOA match data.
●	Fix issues with apostrophes in the JSON-like data.
●	Convert extracted data into a DataFrame and calculate the mode of "match" values for each therapy-MOA pair.
●	Create a pivot table with therapies as rows and MOAs as columns, indicating matches.
●	Save processed therapy-MOA matches to a CSV file.
●	Load an existing TCGA dataset and map therapies to their corresponding MOAs.
●	Explode the DataFrame to create individual rows for each therapy-MOA pair.
●	Save the final processed dataset for analysis, integrating therapy, MOA, and chromosome-level aneuploidy data.
TCGA_drug_survival_curves.py
Analyzes the relationship between survival outcomes, chromosome arm-level aneuploidy, and treatments or MOAs across different cancer types using Kaplan-Meier survival curves and log-rank tests.
Tasks:
●	Load a TCGA dataset containing survival data, chromosome arm-level aneuploidy statuses, and treatments or MOAs.
●	Specify the type of analysis (treatment or MOAs) via the therapy_col variable.
●	Filter the dataset to retain rows with non-missing survival times and chromosome arm data.
●	Group patients into "Loss" or "Other" categories for each chromosome arm.
●	For each cancer type, treatment or MOA, and chromosome arm:
○	Perform Kaplan-Meier survival analysis.
○	Calculate log-rank test p-values.
●	Collect analysis results, including median survival times, p-values, and sample sizes.
●	Apply the Benjamini-Hochberg procedure for FDR correction.
●	Save final results to a CSV file with a filename reflecting the type of analysis.
TCGA_MOA_survival_curves.py
Evaluates the relationship between survival outcomes and chromosome arm-level aneuploidy in cancer patients, stratified by therapy or MOA using Kaplan-Meier survival curves and log-rank tests.
Tasks:
●	Load a TCGA dataset containing survival, treatment, and chromosome arm-level aneuploidy data.
●	Specify analysis based on therapies or MOAs.
●	Filter the dataset to include rows with non-missing survival data.
●	Group patients by chromosome arm status ("Loss" or "Other") for all chromosomes of interest.
●	For each cancer type, therapy or MOA, and chromosome arm:
○	Perform Kaplan-Meier survival analysis.
○	Conduct log-rank tests.
●	Collect key metrics, including median survival times, p-values, and sample sizes.
●	Apply Benjamini-Hochberg correction for p-values.
●	Save results to a CSV file with adjusted p-values and survival statistics.
TCGA_drug_survival_curves_plot.py
Analyzes survival outcomes for cancer patients based on chromosome arm-level aneuploidy and specific drug treatments using Kaplan-Meier survival curves and log-rank tests.
Tasks:
●	Load a preprocessed TCGA dataset with expanded treatment and chromosome arm-level data.
●	Create binary groupings for each chromosome arm ('Loss' or 'Other') based on aneuploidy status.
●	Clean survival data by refactoring survival status into binary terms and removing rows with missing data.
●	Iterate through predefined analyses specifying cancer type, drug treatment, and chromosome arm:
○	Perform log-rank tests.
○	Plot Kaplan-Meier survival curves.
○	Save plots as high-resolution images.
○	Extract survival statistics, including median survival time and sample sizes.
●	Save analysis results to a CSV file.
TCGA_MOA_survival_curves_plot.py
Analyzes survival outcomes for cancer patients treated with specific MOAs and examines the impact of chromosome arm-level aneuploidy using Kaplan-Meier survival curves and log-rank tests.
Tasks:
●	Load a preprocessed TCGA dataset with treatment and chromosome arm-level aneuploidy data.
●	Filter the dataset based on survival time availability.
●	Create binary groups ('Loss' or 'Other') for each chromosome arm based on aneuploidy status.
●	Refactor survival outcomes into binary values and remove rows with missing data.
●	Iterate through predefined analyses specifying cancer type, MOA, and chromosome arm:
○	Perform Kaplan-Meier survival analysis.
○	Conduct log-rank tests.
○	Save survival curves as high-resolution images.
○	Collect results, including p-values, median survival times, and sample sizes.
●	Output a CSV file with analysis results, including MOA, chromosome arm, cancer type, survival statistics, and plot filenames.
TCGA_aneuploidy_counts_by_cancer_and_chromosome.py
Calculates and normalizes chromosome arm-level aneuploidy counts for each cancer type from a TCGA dataset.
Tasks:
●	Load a TCGA dataset containing chromosome arm-level aneuploidy data and treatment information.
●	Remove the "moas" column and duplicates based on unique combinations of sample ID, treatment, and study ID.
●	For each cancer type and chromosome arm:
○	Calculate the number of samples with gains, losses, unchanged states, and missing values.
●	Create a DataFrame with counts for each cancer type and chromosome arm.
●	Normalize counts by dividing by the total number of samples in each cancer type.
●	Save two CSV files:
○	Raw counts of aneuploidy states (TCGA_aneuploidy_counts_by_cancer_type.csv).
○	Normalized counts (TCGA_normalized_aneuploidy_counts_by_cancer_type.csv).
TCGA_circos_plot.py and TCGA_negative_circos_plot.py
Generates a circos plot to visualize relationships between chromosome arm-level aneuploidy, cancer types, and treatments based on survival prognosis data from TCGA.
Tasks:
●	Load a dataset containing survival prognosis associations and cancer type abbreviations.
●	Merge datasets to replace cancer type names with abbreviations.
●	Filter data to include chromosome arms of interest, cancer types, and treatments.
●	Group data into sectors: Aneuploidy (chromosome arms), Cancer Types, Treatments.
●	Define colors for connections and create a "from-to" table for relationships.
●	Calculate sector order, spacing, and assign colors.
●	Configure and initialize the circos plot:
○	Map sectors to groups and assign colors.
○	Configure connections between sectors.
○	Add labels and annotations.
●	Generate and save the plot as a high-resolution PNG file.
Match_common_moa_TCGA_Broad.py
Compares TCGA and Broad datasets to identify shared associations between MOAs, chromosome arms, and cancer types. Integrates survival prognosis data from TCGA with viability data from Broad.
Tasks:
●	Define a mapping dictionary to align cancer types in TCGA with Broad datasets.
●	Load TCGA and Broad datasets into pandas DataFrames.
●	Map cancer types in TCGA using the mapping dictionary and expand rows for multiple mappings.
●	Standardize the MOA column in both datasets.
●	Identify common associations based on MOA, arm, and cancer type.
●	Filter datasets to retain rows with common associations.
●	Merge filtered datasets on MOA, arm, and cancer type.
●	Select relevant columns from TCGA and Broad datasets.
●	Rename and organize selected columns.
●	Save the merged dataset to a CSV file.
Broad_aneuploidy_counts.py
Processes the Broad dataset to summarize chromosome arm-level aneuploidy states for each cancer type.
Tasks:
●	Read the Broad dataset containing cancer type and chromosome arm-level data.
●	Remove duplicate rows based on the cell_line column.
●	Transform the dataset by melting chromosome arm columns into rows.
●	Group data by cancer type, chromosome arm, and aneuploidy state, counting occurrences.
●	Pivot grouped data to organize counts of losses, gains, and unchanged states.
●	Rename columns for clarity and replace missing values with 0.
●	Convert counts to integers.
●	Save the summarized table to a CSV file.
Broad_drug_ANOVA_aggregated_and_unaggregated_by_cancer_type.py and Broad_MOA_ANOVA_aggregated_and_unaggregated_by_cancer_type.py
Performs statistical analysis to examine the relationship between drug/MOA treatments, chromosome arm-level aneuploidy, and cell viability across different cancer types using ANOVA tests.
Tasks:
●	Configure the analysis environment, including file paths, chromosome arms, and cutoff for minimum group size.
●	Define utility functions for loading data, preprocessing, and running statistical analyses.
●	Load the Broad dataset, preprocess data to remove duplicates, and filter by treatment column.
●	Analyze the dataset by cancer type using parallel processing:
○	For each cancer type:
■	Filter data by cancer type and treatments.
■	For each chromosome arm:
■	Compute group sizes.
■	Perform ANOVA tests.
■	Determine sensitivity changes.
●	Collect results for all treatments and chromosome arms.
●	Save cancer-type-specific results to a CSV file.
●	Apply Benjamini-Hochberg correction for p-values.
●	Perform aggregate analysis across all cancer types.
●	Save aggregate results to a separate CSV file.
Output:
●	Cancer-type-specific results showing relationships between drug treatments, chromosome arms, and cell viability.
●	Aggregate results across all cancer types.
Broad_drug_violin_plot_builder.py
Analyzes the relationship between chromosome arm-level aneuploidy and cell viability for specific treatments and cancer types. Generates violin plots for significant ANOVA test results.
Tasks:
●	Read the dataset and preprocess by removing duplicate rows based on cell_line and drug_name.
●	Filter data for a specific cancer type and treatment.
●	Select a specific chromosome arm and create categorical labels ("Loss" or "Other").
●	Verify that both groups have more than 10 observations.
●	Perform ANOVA to test cell viability differences between groups.
●	If p-value < 0.05, generate a violin plot:
○	Customize and save the plot with a formatted title.
●	Iterate through predefined analyses specifying cancer type, treatment, and chromosome arm.
Output: Saved violin plots for significant analyses and logs of data processing and statistical tests.
Broad_top_MOA_by_loss_event.py
Identifies the most common MOAs associated with lower or greater cell viability in cell lines exhibiting chromosome arm loss based on the Broad MOA viability dataset.
Tasks:
●	Read the dataset containing viability shift data for MOAs by chromosome arm.
●	Group data by arm, viability shift, and MOA, counting occurrences.
●	Define a function to identify the most common MOA(s) for each arm and viability shift combination.
●	Apply the function to extract top MOAs for each group.
●	Separate results into:
○	Sensitizing MOAs (lower viability in loss lines).
○	Desensitizing MOAs (greater viability in loss lines).
●	Combine results into a single table.
●	Save the final results to a CSV file.
Output: CSV file listing the top MOAs for each chromosome arm and viability shift category.
Broad_MOA_radial_network.py
Generates a polar plot to visualize top MOAs associated with lower or greater cell viability in cell lines with chromosome arm loss.
Tasks:
●	Read the dataset containing top MOAs for each chromosome arm and viability shift.
●	Process data to create ordered lists of MOAs, counts, and chromosome arms for both "positive" and "negative" groups.
●	Normalize radii for visual scaling and map color gradients to chromosome arms.
●	Set up a polar coordinate plot:
○	Define angles and base radii for categories.
○	Plot radiating lines with colors and lengths representing data.
○	Add text annotations for MOA, count, and chromosome arm.
●	Configure plot aesthetics.
●	Save the plot as a high-resolution TIFF file.
Broad_cell_line_sensitivity_cancer_dot_plot.py
Generates a broken-axis scatter plot to visualize significant drug viability shifts by cancer type based on p-values. Separates results into reduced sensitivity and greater sensitivity groups.
Tasks:
●	Read a CSV file containing drug viability association data.
●	Preprocess data:
○	Map full cancer type names to abbreviations.
○	Transform p-values into -log10(p-value).
○	Update sensitivity labels and map to colors.
●	Apply jitter to x-axis positions based on y-axis values.
●	Create two subplots with a broken y-axis:
○	One for positive -log10(p-value) values.
○	One for negative -log10(p-value) values.
●	Add color-coded data points and highlight specific drugs.
●	Add custom diagonal lines indicating the y-axis break.
●	Add legends for sensitivity shifts and highlighted drugs.
●	Save the plot as a high-resolution PNG file.
Broad_cell_line_sensitivity_comparison_plot.py
Generates a scatter plot to visualize significant associations between therapies and chromosome arms based on p-values.
Tasks:
●	Load and preprocess data:
○	Read a CSV file containing drug viability associations by chromosome arm.
○	Compute -log10(p-value).
○	Assign sensitivity types with specific colors.
○	Retain and order chromosome arms as categorical data.
●	Apply jitter and offsets to separate overlapping points.
●	Define markers and colors for specific drugs.
●	Create scatter plot with jittered x-axis positions and differentiated colors.
●	Set dynamic x-axis labels based on filtered chromosome arms.
●	Customize title, labels, and legend.
●	Save the plot as "dot_plot_significant_Broad_drugs_arm_jittered.png" or display it interactively.

R script:
CCLE_Aneuploidy
This folder contains input and output data from scripts designed to calculate aneuploidy from the .seg data provided by CCLE curators.
arm_edges.tsv lists the protein coding genes and locations that define beginning and ends of chromosome arms, for the purposes of the aneuploidy analysis here.
data_CCLEcna_2012_Nature_hg19.seg contains the copy-number variant calls per cell line downloaded from the CCLE project.
"Repurposing" files are the DepMap project files downloaded for this analysis. 
CCLE_chromosome_CNA_caller.R thresholds the .seg data to produce an aneuploidy call file, "CCLE_arm_signed_aneuploidy.tsv" 
