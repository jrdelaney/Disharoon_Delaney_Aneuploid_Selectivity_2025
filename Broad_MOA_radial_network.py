import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.patheffects as patheffects

# Load data
df = pd.read_csv('../../Final_Publication_Files/Supplemental Table 11 Broad Top Cytotoxic MOAs.csv')
df['chr_arm'] = df['arm'].str.replace('chr', '', regex=False)

# Define the desired order of chromosome arms
CHROMOSOME_ARMS = [
    '1p', '1q', '2p', '2q', '3p', '3q', '4p', '4q', '5p', '5q',
    '6p', '6q', '7p', '7q', '8p', '8q', '9p', '9q', '10p', '10q',
    '11p', '11q', '12p', '12q', '13q', '14q', '15q', '16p', '16q',
    '17p', '17q', '18p', '18q', '19p', '19q', '20p', '20q', '21q', '22q'
]

# Convert 'chr_arm' column to categorical with the specified order
df['chr_arm'] = pd.Categorical(df['chr_arm'], categories=CHROMOSOME_ARMS, ordered=True)

# Separate labels and counts by association and sort by 'chr_arm'
positive_df = df[df['viability_shift'] == 'Lower viability in loss lines'].sort_values(by='chr_arm')
negative_df = df[df['viability_shift'] == 'Greater viability in loss lines'].sort_values(by='chr_arm')

# Create ordered lists with placeholders for positive values
positive_labels = []
positive_counts = []
positive_arms = []

for arm in CHROMOSOME_ARMS:
    if arm in positive_df['chr_arm'].values:
        positive_labels.append(positive_df[positive_df['chr_arm'] == arm]['MOA'].values[0])
        positive_counts.append(positive_df[positive_df['chr_arm'] == arm]['count_of_significant_associations'].values[0])
        positive_arms.append(arm)
    else:
        positive_labels.append('None')
        positive_counts.append(3)
        positive_arms.append(arm)

# Create ordered lists with placeholders for negative values
negative_labels = []
negative_counts = []
negative_arms = []

for arm in CHROMOSOME_ARMS:
    if arm in negative_df['chr_arm'].values:
        negative_labels.append(negative_df[negative_df['chr_arm'] == arm]['MOA'].values[0])
        negative_counts.append(negative_df[negative_df['chr_arm'] == arm]['count_of_significant_associations'].values[0])
        negative_arms.append(arm)
    else:
        negative_labels.append('None')
        negative_counts.append(3)
        negative_arms.append(arm)

# Normalize count values to set them within a specific radius range
min_radius = 0
positive_radii = np.array(positive_counts) # + min_radius
negative_radii = np.array(negative_counts) # + min_radius

# Define custom color maps for Reds and Blues with a darker starting color
reds_colormap = cm.get_cmap('Reds')
blues_colormap = cm.get_cmap('Blues')
custom_reds = mcolors.LinearSegmentedColormap.from_list('custom_reds', reds_colormap(np.linspace(0.3, 1, 256)))
custom_blues = mcolors.LinearSegmentedColormap.from_list('custom_blues', blues_colormap(np.linspace(0.3, 1, 256)))

# Assign colors to chromosome arms in sequence based on custom colormaps
positive_color_map = [custom_reds(i / (len(CHROMOSOME_ARMS) - 1)) for i in range(len(CHROMOSOME_ARMS))]
negative_color_map = [custom_blues(i / (len(CHROMOSOME_ARMS) - 1)) for i in range(len(CHROMOSOME_ARMS))]

# Calculate the angle for each label
angles_positive = np.linspace(-3 * np.pi / 8, 3 * np.pi / 8, len(positive_labels), endpoint=True)
angles_negative = np.linspace(5 * np.pi / 8, 11 * np.pi / 8, len(negative_labels), endpoint=True)

# Initialize the plot in polar coordinates
fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={'projection': 'polar', 'frame_on': False})
ax.set_theta_offset(0)
ax.set_theta_direction(-1)

# Set base radius to start lines farther out from the center
base_radius = 3

# Update overall maximum radius to include the base radius
overall_max = base_radius + max(np.max(positive_radii), np.max(negative_radii))

# Adjust radial limits
ax.set_rmin(base_radius - 1)  # Optional: subtract 1 for visual padding
ax.set_rmax(overall_max + 1)  # Optional: add 1 for visual padding

# Plot positive labels with colors based on sequence
for angle, label, radius, color, count, arm in zip(
    angles_positive, positive_labels, positive_radii, positive_color_map, positive_counts, positive_arms
):
    linestyle = (0, (0.25, 0.25)) if count == 0 else '-'
    ax.plot(
        [angle, angle],
        [base_radius, base_radius + radius],
        color=color,
        linewidth=8,
        linestyle=linestyle
    )
    full_label = f"   {label} | {arm}   "
    ax.text(
        angle,
        base_radius + radius,
        full_label,
        ha='left',
        va='center',
        fontsize=6,
        rotation=90,
        rotation_mode='anchor',
        transform_rotates_text=True,
        path_effects=[patheffects.withStroke(linewidth=0.75, foreground="white")]
    )

# Plot negative labels with colors based on sequence
for angle, label, radius, color, count, arm in zip(
    angles_negative, negative_labels, negative_radii, negative_color_map, negative_counts, negative_arms
):
    linestyle = (0, (0.25, 0.25)) if count == 0 else '-'
    ax.plot(
        [angle, angle],
        [base_radius, base_radius + radius],
        color=color,
        linewidth=8,
        linestyle=linestyle
    )
    full_label = f"   {arm} | {label}   "
    ax.text(
        angle,
        base_radius + radius,
        full_label,
        ha='right',
        va='center',
        fontsize=6,
        rotation=270,
        rotation_mode='anchor',
        transform_rotates_text=True,
        path_effects=[patheffects.withStroke(linewidth=0.75, foreground="white")]
    )

# Adjust radial grid
tick_positions = np.linspace(base_radius, overall_max, num=5)
ax.set_rgrids(
    tick_positions,
    labels=np.round(tick_positions).astype(int)-base_radius,
    angle=270,
    rotation=270,
    fontsize=8,
    color='gray'
)

# Remove polar grid and axis ticks for a cleaner look
ax.set_xticks([])

plt.savefig('best_MOA_by_aneuploidy_explosionplot.tif', dpi=600, bbox_inches='tight')
plt.show()
