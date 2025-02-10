import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the CSV file
file_path = "Final_Publication_Files/Supplemental Table 3 Broad Drug Viability Associations by Cancer Type.csv"
df = pd.read_csv(file_path)

# Ensure the required columns are present
required_columns = ['p_value', 'viability_shift', 'arm']
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

# Process chromosome arms
df['arm'] = df['arm'].str.replace('chr', '', regex=False)

# Original order of chromosome arms
CHROMOSOME_ARMS = [
    '1p', '1q', '2p', '2q', '3p', '3q', '4p', '4q', '5p', '5q',
    '6p', '6q', '7p', '7q', '8p', '8q', '9p', '9q', '10p', '10q',
    '11p', '11q', '12p', '12q', '13q', '14q', '15q', '16p', '16q',
    '17p', '17q', '18p', '18q', '19p', '19q', '20p', '20q', '21q', '22q'
]

# Filter CHROMOSOME_ARMS to only include those with data
filtered_chromosome_arms = [arm for arm in CHROMOSOME_ARMS if arm in df['arm'].unique()]

# Convert 'arm' column to categorical with the specified filtered order
df['arm'] = pd.Categorical(df['arm'], categories=filtered_chromosome_arms, ordered=True)

# Map categorical chromosome arms to numeric positions
df['arm_numeric'] = df['arm'].cat.codes

# Define offsets for sensitivity categories
offsets = {
    'Reduced Sensitivity': -0.2,  # Shift left
    'Greater Sensitivity': 0.2   # Shift right
}

# Apply jitter and offset
np.random.seed(42)  # For reproducibility
jitter_strength = 0.1
df['arm_jittered'] = df.apply(
    lambda row: row['arm_numeric'] + np.random.uniform(-jitter_strength, jitter_strength) + offsets[row['viability_shift']],
    axis=1
)

# Get the minimum value for the y-axis
y_min = df['-log10(p-value)'].min()

# Create the plot
plt.figure(figsize=(14, 8))
# Plot points with jitter for each sensitivity type, excluding specific drug values
for sensitivity in df['viability_shift'].unique():
    subset = df[(df['viability_shift'] == sensitivity) & ~df['drug'].isin(['ETOMIDATE', 'PF-06651600', 'PITAVASTATIN'])]
    plt.scatter(
        subset['arm_jittered'],
        subset['-log10(p-value)'],
        label=sensitivity,
        color=sensitivity_colors[sensitivity],
        alpha=0.5
    )

# Plot specific drug values with custom markers and colors
etomidate_data = df[df['drug'] == 'ETOMIDATE']
pf_06651600_data = df[df['drug'] == 'PF-06651600']
pita_data = df[df['drug'] == 'PITAVASTATIN']

# Yellow green diamonds
plt.scatter(
    pita_data['arm_jittered'],
    pita_data['-log10(p-value)'],
    color='yellowgreen',
    marker='D',
    s=50,
    label='Pitavastatin'
)

# Red triangles for ETOMIDATE
plt.scatter(
    etomidate_data['arm_jittered'],
    etomidate_data['-log10(p-value)'],
    color='yellow',
    marker='^',
    s=50,
    label='Etomidate'
)

# Purple squares for PF-06651600
plt.scatter(
    pf_06651600_data['arm_jittered'],
    pf_06651600_data['-log10(p-value)'],
    color='aqua',
    marker='s',
    s=50,
    label='Pf-06651600'
)


# Dynamically set x-axis labels based on filtered and ordered chromosome arms
plt.xticks(
    ticks=np.arange(len(filtered_chromosome_arms)),
    labels=filtered_chromosome_arms,
    rotation=90,
    ha='right',
    fontsize=14
)

# Customize the plot
plt.title('Significant Broad Therapies')
plt.xlabel('Chromosome Arm')
plt.ylabel('-log10(p-value)')
plt.ylim(y_min, None)  # Set y-axis lower limit to dataset minimum
plt.xlim(-0.5, len(filtered_chromosome_arms) - 0.5)
legend = plt.legend(title='Viability', loc='upper left', frameon=False)
legend.set_title(None)

plt.tight_layout()

# Save or display the plot
plt.savefig("dot_plot_significant_Broad_drugs_arm_jittered.png")
plt.show()
