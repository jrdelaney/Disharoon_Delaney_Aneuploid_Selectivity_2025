import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D  # Import for custom legend handles
from matplotlib.ticker import MultipleLocator, FormatStrFormatter  # Import for tick formatting

# Define the cancer type names and mapping
cancer_type_mapper = {
    "Breast": "BC",
    "Central Nervous System": "CNSC",
    "Haematopoietic And Lymphoid Tissue": "HLTC",
    "Kidney": "KIDC",
    "Large Intestine": "LIC",
    "Liver": "LIVC",
    "Lung": "LC",
    "Oesophagus": "ESC",
    "Ovary": "OVC",
    "Pancreas": "PC",
    "Skin": "SC",
    "Stomach": "STC",
    "Upper Aerodigestive Tract": "UATC",
    "Urinary Tract": "UTC"
}

# Load the CSV file
file_path = "Final_Publication_Files/Supplemental Table 3 Broad Drug Viability Associations by Cancer Type.csv"
df = pd.read_csv(file_path)

# Ensure the required columns are present
required_columns = ['p_value', 'viability_shift', 'cancer_type', 'drug']
if not all(col in df.columns for col in required_columns):
    raise ValueError(f"The CSV file must contain these columns: {', '.join(required_columns)}")

# Transform p-values to -log10(p-value)
df['-log10(p-value)'] = -np.log10(df['p_value'])

# Update viability_shift labels
df['viability_shift'] = df['viability_shift'].map({
    'Greater viability in loss lines': 'Reduced Sensitivity',
    'Lower viability in loss lines': 'Greater Sensitivity'
})

# Ensure viability_shift is categorical
df['viability_shift'] = df['viability_shift'].astype('category')

# Define colors for sensitivity types
sensitivity_colors = {
    'Reduced Sensitivity': 'blue',
    'Greater Sensitivity': 'red'
}

# Map sensitivity types to colors
df['color'] = df['viability_shift'].map(sensitivity_colors)

# Multiply '-log10(p-value)' by -1 for 'Reduced Sensitivity' to make them negative
df.loc[df['viability_shift'] == 'Reduced Sensitivity', '-log10(p-value)'] *= -1

# Remove BONE and SOFT_TISSUE from the dataset due to few significant results
df = df[~df['cancer_type'].isin(['BONE', 'SOFT_TISSUE'])]

# Replace underscores with spaces and convert to title case for better labels
df['cancer_type'] = df['cancer_type'].str.replace('_', ' ').str.title()

# Map cancer types to their abbreviations using cancer_type_mapper
df['cancer_type_abbrev'] = df['cancer_type'].map(cancer_type_mapper)

# For any cancer types not in the mapper, keep the original name
df['cancer_type_abbrev'] = df['cancer_type_abbrev'].fillna(df['cancer_type'])

# Get unique abbreviated cancer types in desired order
abbreviated_cancer_types = sorted(df['cancer_type_abbrev'].unique())

# Convert 'cancer_type_abbrev' column to categorical with the specified order
df['cancer_type_abbrev'] = pd.Categorical(df['cancer_type_abbrev'], categories=abbreviated_cancer_types, ordered=True)

# Map categorical cancer types to numeric positions
df['cancer_type_numeric'] = df['cancer_type_abbrev'].cat.codes

# Apply variable jitter to the x-axis positions to create pyramid effect
np.random.seed(42)  # For reproducibility

max_jitter = 1    # Maximum jitter at y=0
min_jitter = 0      # Minimum jitter at maximum |y|
max_abs_y = df['-log10(p-value)'].abs().max()
relative_y = df['-log10(p-value)'].abs() / max_abs_y

# Define the decay rate for visual peaks
decay_rate = 5

# Calculate jitter_strength with exponential decay
df['jitter_strength'] = min_jitter + (max_jitter - min_jitter) * np.exp(-decay_rate * relative_y)

# Apply jitter with variable strength
df['cancer_type_jittered'] = df['cancer_type_numeric'] + np.random.uniform(-1, 1, len(df)) * df['jitter_strength']

# Get the minimum positive and maximum negative values for the y-axis
positive_values = df['-log10(p-value)'][df['-log10(p-value)'] > 0]
negative_values = df['-log10(p-value)'][df['-log10(p-value)'] < 0]

y_pos_min = positive_values.min()
y_pos_max = positive_values.max()
y_neg_max = negative_values.max()
y_neg_min = negative_values.min()

# Set the break points dynamically
upper_ylim = y_pos_max + 0.25
upper_ymin = y_pos_min - 0.05
lower_ylim = y_neg_max + 0.05
lower_ymin = y_neg_min - 0.25

# Create the figure and subplots (reduced figure width)
fig = plt.figure(figsize=(10, 6))  # Reduced width from 14 to 10
gs = gridspec.GridSpec(2, 1, height_ratios=[upper_ylim - upper_ymin, lower_ylim - lower_ymin], hspace=0.05)

# Upper plot (positive values)
ax_upper = plt.subplot(gs[0])
scatter_upper = ax_upper.scatter(
    df['cancer_type_jittered'][df['-log10(p-value)'] > 0],
    df['-log10(p-value)'][df['-log10(p-value)'] > 0],
    c=df['color'][df['-log10(p-value)'] > 0],
    alpha=0.5
)
# Plot specific drug values with custom markers and colors in the upper plot
for drug_data, color, marker, label in [
    (df[(df['drug'] == 'PITAVASTATIN') & (df['-log10(p-value)'] > 0)], 'yellowgreen', 'D', 'Pitavastatin'),
    (df[(df['drug'] == 'ETOMIDATE') & (df['-log10(p-value)'] > 0)], 'yellow', '^', 'Etomidate'),
    (df[(df['drug'] == 'PF-06651600') & (df['-log10(p-value)'] > 0)], 'aqua', 's', 'Pf-06651600'),

]:
    ax_upper.scatter(
        drug_data['cancer_type_jittered'],
        drug_data['-log10(p-value)'],
        color=color,
        marker=marker,
        s=50,
        label=label
    )
ax_upper.set_ylim(upper_ymin, upper_ylim)
ax_upper.spines['bottom'].set_visible(False)
ax_upper.tick_params(labeltop=False)
# Hide the x ticks on the top
ax_upper.set_xticks([])

# Set Y-axis major ticks at intervals of 1 and format as integers
ax_upper.yaxis.set_major_locator(MultipleLocator(1))
ax_upper.yaxis.set_major_formatter(FormatStrFormatter('%d'))

# Lower plot (negative values)
ax_lower = plt.subplot(gs[1], sharex=ax_upper)
scatter_lower = ax_lower.scatter(
    df['cancer_type_jittered'][df['-log10(p-value)'] < 0],
    df['-log10(p-value)'][df['-log10(p-value)'] < 0],
    c=df['color'][df['-log10(p-value)'] < 0],
    alpha=0.5
)
# Plot specific drug values with custom markers and colors in the lower plot
for drug_data, color, marker, label in [
    (df[(df['drug'] == 'PITAVASTATIN') & (df['-log10(p-value)'] < 0)], 'yellowgreen', 'D', 'Pitavastatin'),
    (df[(df['drug'] == 'ETOMIDATE') & (df['-log10(p-value)'] < 0)], 'yellow', '^', 'Etomidate'),
    (df[(df['drug'] == 'PF-06651600') & (df['-log10(p-value)'] < 0)], 'aqua', 's', 'PF-06651600'),
]:
    ax_lower.scatter(
        drug_data['cancer_type_jittered'],
        drug_data['-log10(p-value)'],
        color=color,
        marker=marker,
        s=50,
        label=label
    )
ax_lower.set_ylim(lower_ymin, lower_ylim)
ax_lower.spines['top'].set_visible(False)
ax_lower.xaxis.tick_bottom()

# Set Y-axis major ticks at intervals of 1 and format as integers
ax_lower.yaxis.set_major_locator(MultipleLocator(1))
ax_lower.yaxis.set_major_formatter(FormatStrFormatter('%d'))

# Adjust x-axis labels
ax_lower.set_xticks(np.arange(len(abbreviated_cancer_types)))
ax_lower.set_xticklabels(abbreviated_cancer_types, rotation=90, ha='right', fontsize=14)

# Hide x-axis labels on the upper plot
plt.setp(ax_upper.get_xticklabels(), visible=False)

# Draw diagonal lines to indicate the break
d = .015  # Size of diagonal lines in axes coordinates
kwargs = dict(transform=ax_upper.transAxes, color='k', clip_on=False)
ax_upper.plot((-d, +d), (-d, +d), **kwargs)        # Top-left diagonal
ax_upper.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # Top-right diagonal

kwargs.update(transform=ax_lower.transAxes)  # Switch to the lower axes
ax_lower.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # Bottom-left diagonal
ax_lower.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # Bottom-right diagonal

# Create custom legend handles for sensitivity types
custom_handles = [
    Line2D([], [], marker='o', color='w', markerfacecolor='blue', markersize=10, label='Greater Viability'),
    Line2D([], [], marker='o', color='w', markerfacecolor='red', markersize=10, label='Reduced Viability'),
    Line2D([], [], marker='^', color='yellow', markersize=10, label='Etomidate'),
    Line2D([], [], marker='s', color='aqua', markersize=10, label='Pf-06651600'),
    # Line2D([], [], marker='D', color='yellowgreen', markersize=10, label='Alfacalcidol'),
    Line2D([], [], marker='D', color='yellowgreen', markersize=10, label='Pitavastatin'),
]

# Add legends
legend = ax_lower.legend(handles=custom_handles, title='Viability', loc='lower right', frameon=False)
legend.set_title(None)

# Set labels
ax_upper.set_ylabel('-log10(p-value)')
ax_lower.set_ylabel('log10(p-value)')
ax_lower.set_xlabel('Cancer Type')

# Set titles
ax_upper.set_title('Significant Broad Therapies by Cancer Type')

# Tight layout
plt.tight_layout()

# Save or display the plot
plt.savefig("dot_plot_significant_Broad_drugs_cancer_type_broken_axis.png", dpi=600, bbox_inches='tight')
plt.show()
