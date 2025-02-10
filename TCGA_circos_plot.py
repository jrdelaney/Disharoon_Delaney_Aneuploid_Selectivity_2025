import pandas as pd
from pycirclize import Circos
from pycirclize.parser import Matrix
from pycirclize.utils import calc_group_spaces, ColorCycler

# Define the plot name for saving the figure
plot_name = 'positive_circos_all_survival'

# Load the survival curve results data
data = pd.read_csv('../Final_Publication_Files/Supplemental Table 1 TCGA Positive Survival Prognosis Associations.csv')

# Load a file containing cancer type abbreviations
# NOTE: This is the TCGA specific portion of Supplemental Table 9 Cancer Type Abbreviations.csv
cancer_abbreviations = pd.read_csv('cancer_type_abbreviations.csv')

# Merge the data with cancer abbreviations to add the abbreviation for cancer types
data = data.merge(cancer_abbreviations, left_on='cancer_type', right_on='Cancer', how='left')

# Replace the 'cancer' column with its abbreviation
data['cancer_type'] = data['Abbreviation']

# Drop the extra columns from the merge if not needed
data = data.drop(columns=['Cancer', 'Abbreviation'])

# Define chromosomes of interest for analysis in ascending order
chrs_interest = [
    '1p', '1q', '2p', '2q', '3p', '3q', '4p', '4q', '5p', '5q',
    '6p', '6q', '7p', '7q', '8p', '8q', '9p', '9q', '10p', '10q',
    '11p', '11q', '12p', '12q', '13q', '14q', '15q', '16p', '16q',
    '17p', '17q', '18p', '18q', '19p', '19q', '20p', '20q', '21q', '22q'
]

# Identify chromosomes present in the data
chromosomes = [chr for chr in chrs_interest if chr in data['arm'].unique()]

# Count occurrences of each cancer type, sorted in ascending order
cancer_counts = data.groupby('cancer_type').size().sort_values(ascending=True)
cancers = cancer_counts.index.tolist()

# Count occurrences of each treatment, sorted in descending order
treatment_counts = data.groupby('drug').size().sort_values(ascending=False)
treatments = treatment_counts.index.tolist()

# Define groups for sectors: Aneuploidy, Cancer Type, and Treatment
groups = {
    'Aneuploidy': chromosomes,
    'Cancer Type': cancers,
    'Treatment': treatments
}

# Define colors for connections between different groups
connection_colors = {
    ('arm', 'cancer_type'): 'lightsteelblue',
    ('arm', 'drug'): 'lightsteelblue',
    ('cancer_type', 'drug'): 'aquamarine',
    ('cancer_type', 'arm'): 'aquamarine',
    ('drug', 'cancer_type'): 'bisque',
    ('drug', 'arm'): 'bisque'
}

# Create a "from-to" table for connections, including values and colors
fromto_table_df = pd.concat([
    data.groupby(['arm', 'drug']).size()
        .reset_index(name='value')
        .rename(columns={'arm': 'from', 'drug': 'to'})
        .assign(color=connection_colors[('arm', 'drug')]),

    data.groupby(['cancer_type', 'arm']).size()
        .reset_index(name='value')
        .rename(columns={'cancer_type': 'from', 'arm': 'to'})
        .assign(color=connection_colors[('cancer_type', 'arm')]),

    data.groupby(['drug', 'cancer_type']).size()
        .reset_index(name='value')
        .rename(columns={'drug': 'from', 'cancer_type': 'to'})
        .assign(color=connection_colors[('drug', 'cancer_type')])
])

# Define the order of sectors in the plot
ordered_sectors = groups['Aneuploidy'] + groups['Cancer Type'] + groups['Treatment']

# Convert the "from-to" table into a matrix for plotting
matrix = Matrix.parse_fromto_table(fromto_table_df)

# Calculate spaces between groups for visual clarity
group_sizes = [len(groups['Aneuploidy']), len(groups['Cancer Type']), len(groups['Treatment'])]
spaces = calc_group_spaces(group_sizes, space_bw_group=15, space_in_group=2)

# Set up colors for each group
ColorCycler.set_cmap("tab20")  # Choose a colormap
group_color_list = ['blue', 'lightgreen', 'orange']  # Assign a color to each group

# Map each sector to a corresponding group color
group_colors = {}
for group_name, color in zip(groups.keys(), group_color_list):
    for item in groups[group_name]:
        group_colors[item] = color

# Assign colors for links (chords) based on connections
link_cmap = [
    (row['from'], row['to'], row['color']) for _, row in fromto_table_df.iterrows()
]

# Initialize the Circos plot with sectors, spaces, and color mappings
circos = Circos.initialize_from_matrix(
    matrix,
    space=spaces,
    order=ordered_sectors,  # Sector order
    cmap=group_colors,  # Sector colors
    link_cmap=link_cmap,  # Link colors
    label_kws=dict(size=8, orientation="vertical"),  # Label settings
    link_kws=dict(direction=1)  # Link arrow direction
)

# Add group annotations to the plot
for group_name, items in groups.items():
    group_deg_lim = circos.get_group_sectors_deg_lim(items)  # Get angle limits for the group
    circos.rect(r_lim=(97, 100), deg_lim=group_deg_lim, fc=group_colors[items[0]], ec="black", lw=1)  # Draw group indicator
    group_center_deg = sum(group_deg_lim) / 2  # Calculate group center angle
    circos.text(f"{group_name}", r=75, deg=group_center_deg, adjust_rotation=True, size=12)  # Add group label

# Add plot annotation for dataset details
circos.text(f"Improved Prognosis\nN={len(data)}", size=20)

# Generate the plot
fig = circos.plotfig()

# Save the plot to a file with high resolution
fig.savefig(f"{plot_name}_TCGA_chord_plot_filtered.png", dpi=1000)