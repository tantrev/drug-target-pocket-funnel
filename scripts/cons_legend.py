import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import rcParams

# Set font to Arial
rcParams['font.family'] = 'Arial'

# New dimensions
width = 450
height = 38  # 75% of the original 50px

# Create horizontal gradient
gradient = np.linspace(-2, 2, width).reshape(1, -1)
gradient = np.vstack([gradient] * height)

# Define custom blue-white-red colormap
cmap = LinearSegmentedColormap.from_list("custom_bwr", ["#0000ff", "#ffffff", "#ff0000"])

# Create figure and axis
fig, ax = plt.subplots(figsize=(9, 2))  # Adjusted for new width and larger text
ax.imshow(gradient, aspect='auto', cmap=cmap, extent=[-2, 2, 0, 1])

# Draw black border around gradient box
ax.plot([-2, 2, 2, -2, -2], [0, 0, 1, 1, 0], color='black', linewidth=1.5)

# Add large numerical labels with tighter vertical spacing
ax.text(-2, -0.15, '-2\nLess conserved', ha='center', va='top', fontsize=30)
ax.text(0, -0.15, '0', ha='center', va='top', fontsize=30)
ax.text(2, -0.15, '2\nMore conserved', ha='center', va='top', fontsize=30)

# Clean up axes
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')

# Save the figure
plt.savefig("../outputs/conservation_legend.png", bbox_inches='tight', dpi=300)
plt.close()
