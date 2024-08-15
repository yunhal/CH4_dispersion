import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate


filename = "/Users/yunhalee/Documents/methanDart/Gaussian_Puff_CH4/output_data/gplume_conc/gplume_output_WS_1.5_stackH_30_stability_['B']_radiaton_weak_WD_20_spread_20.npy"
methane_data = np.load(filename)


def apply_wind_effect(methane_data, wind_direction):
    """
    Apply the wind direction effect to the methane image data.
    
    Parameters: 
    methane_data: 2D array of methane transmittance data
    wind_direction: Wind direction in degrees (0-360)

    Return: Rotated methane data to simulate wind direction change
    """
    # Rotate the methane data based on wind direction
    rotated_methane = rotate(methane_data, angle=-wind_direction, reshape=False, mode='nearest')
    return rotated_methane

# Simulate for a specific wind direction (e.g., 45 degrees)
wind_direction = 180
methane_with_wind = apply_wind_effect(methane_data, wind_direction)


# Display the results
plt.figure(figsize=(10, 10))

plt.subplot(1, 3, 1)
plt.title("Methane")
plt.imshow(np.mean(methane_data, axis = -1).T, cmap='jet')

plt.subplot(1, 3, 2)
plt.title("Methane (Wind Effect)")
plt.imshow(np.mean(methane_with_wind, axis = -1).T, cmap='jet')


plt.show()