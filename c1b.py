import math
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.colors as mcolors
import numpy as np

frequencyaxis = range(405, 790)
frequency = [f*(10**12) for f in frequencyaxis]
# n = [math.sqrt(1+math.sqrt(1/(1.731-0.261*(f/10**15)**2))) for f in frequency]
n = [math.sqrt(1+math.sqrt(1/(1.731-0.261*(f/10**15)**2))) for f in frequency]

def freq_to_colours(f):
    if f < 405: 
        return None, None, None, 'Infra Red'
    elif 405 <= f < 480:
        return 1, 0, 0, 'Red'
    elif 480 <= f < 510:
        return 1, 127/255, 0, 'Orange'
    elif 510 <= f < 530:
        return 1, 1, 0, 'Yellow'
    elif 530 <= f < 600:
        return 0, 1, 0, 'Green'
    elif 600 <= f < 620:
        return 0, 1, 1, 'Cyan'
    elif 620 <= f < 680:
        return 0, 0, 1, 'Blue'
    elif 680 <= f <= 790:
        return 127/255, 0, 1, 'Violet'
    else:
        return None, None, None, 'Ultra Violet'


fig, ax = plt.subplots()

x_data = np.array(frequencyaxis)
y_data = np.array(n)

points = np.array([x_data, y_data]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

spectrum_colors = [
        '#8B00FF',  # Violet (400nm ≈ 750THz)
        '#0000FF',  # Blue (450nm ≈ 667THz)
        '#00FFFF',  # Cyan (490nm ≈ 612THz)
        '#00FF00',  # Green (530nm ≈ 566THz)
        '#FFFF00',  # Yellow (570nm ≈ 526THz)
        '#FF7F00',  # Orange (600nm ≈ 500THz)
        '#FF0000'   # Red (700nm ≈ 429THz)
    ]


# Create custom colormap
spectrum_cmap = mcolors.LinearSegmentedColormap.from_list('spectrum', spectrum_colors)
    
# Create LineCollection
lc = LineCollection(segments, cmap=spectrum_cmap, linewidth=2)
lc.set_array(x_data[:-1])  # Color based on frequency
    
line = ax.add_collection(lc)
ax.set_xlim(x_data.min(), x_data.max())
ax.set_ylim(y_data.min() - 0.01, y_data.max() + 0.01)
    
# Add colorbar
cbar = plt.colorbar(line, ax=ax)
cbar.set_label('Frequency (THz)', rotation=270, labelpad=15)
    
ax.set_xlabel("Frequency (THz)")
ax.set_ylabel("Refractive Index")
ax.set_title("Refractive Index vs Frequency - Smooth Spectrum Gradient")
ax.grid(True, alpha=0.3)
    
plt.tight_layout()
plt.show()
# ax.plot(frequencyaxis, n)
# plt.xlabel("frequency THz")
# plt.ylabel("refractive index")
# plt.show()


# initla problems thinking icould model the whole range from 0 to 1500 thz, missing the key information that the formula was only valid for visible light spectrum
# after fixing this, slightly bamboozled by converting between hz and thz and back again, so i seperated out axis numbers and numbers for use 

# explore matplotlibs linecollections, linear segmented colormaps, colorbar integration