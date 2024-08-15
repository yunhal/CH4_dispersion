__all__ = [
    'get_stab_class',
    'compute_sigma_vals',
    'gplume',
    'plot_2d_plume_concentration',
    'run_simulation',
    'apply_noise',
    'compute_mean_concentration',
    'generate_wind_directions',
    'create_grid'

]


import numpy as np
from datetime import datetime
import math
from pyproj import Proj, transform, CRS
from datetime import timedelta
from math import sqrt
from pyproj import Transformer
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap
from multiprocessing import Pool
import os
from scipy.special import erfcinv
from scipy.ndimage import gaussian_filter

def get_stab_class(U, day_time, solar_radiation=None, cloud_cover_fraction=None):
    """
    Determine stability class based on wind speed, time of day, incoming solar radiation during the daytime, and cloud cover fraction during nighttime.
    """
    if day_time:
        if solar_radiation == 'strong':
            if U < 2:
                return ["A"]
            elif 2 <= U < 3:
                return ["A", "B"]
            elif 3 <= U < 5:
                return ["B"]
            else: #5 <= U 
                return ["C"]

        elif solar_radiation == 'moderate':
            if U < 2:
                return ["A", "B"]
            elif 2 <= U < 3:
                return ["B"]
            elif 3 <= U < 5:
                return ["B", "C"]
            elif 5 <= U < 6:
                return ["C", "D"]
            else:
                return ["D"]
        elif solar_radiation == 'weak':
            if U < 2:
                return ["B"]
            elif 2 <= U < 5:
                return ["C"]
            else: # 5 <= U < 6:
                return ["D"]
    else:
        if cloud_cover_fraction == 'low':
            if U < 3:
                return ["F"]
            elif 3 <= U < 5:
                return ["E"]
            elif 5 <= U < 6:
                return ["D"]
            else:
                return ["D"]
        elif cloud_cover_fraction == 'high':
            if U < 2:
                return ["E"]
            elif 2 <= U < 3:
                return ["E"]
            else: # 3 <= U :
                return ["D"]


def compute_sigma_vals(stab_classes, total_dist):
    """
    Compute disperson coefficients (horizontal and vertical)  based on the Pasquill stability classes
    """
    # Initialize arrays for sigma values
    sigma_y_vals = np.zeros_like(total_dist)
    sigma_z_vals = np.zeros_like(total_dist)

    # Define parameters for different stability classes
    params = {
        'A': [
            (0.1, 122.8, 0.9447, 24.1670, 2.5334),
            (0.15, 158.08, 1.0542, 24.1670, 2.5334),
            (0.20, 170.22, 1.0932, 24.1670, 2.5334),
            (0.25, 179.52, 1.1262, 24.1670, 2.5334),
            (0.3, 217.41, 1.2644, 24.1670, 2.5334),
            (0.4, 258.89, 1.4094, 24.1670, 2.5334),
            (0.5, 346.75, 1.7283, 24.1670, 2.5334),
            (float('inf'), 453.85, 2.1166, 24.1670, 2.5334)
        ],
        'B': [
            (0.2, 90.673, 0.93198, 18.3330, 1.8096),
            (0.3, 98.483, 0.98332, 18.3330, 1.8096),
            (0.4, 109.3, 1.09710, 18.3330, 1.8096),
            (0.5, 112.88, 1.11750, 18.3330, 1.8096),
            (float('inf'), 135.5, 1.1831, 18.3330, 1.8096)
        ],
        'C': [
            (0.2, 61.141, 0.91465, 12.5, 1.0857),
            (0.3, 61.141, 0.91465, 12.5, 1.0857),
            (0.4, 68.55, 0.92627, 12.5, 1.0857),
            (0.5, 78.844, 0.94518, 12.5, 1.0857),
            (float('inf'), 109.3, 1.054, 12.5, 1.0857)
        ],
        'D': [
            (0.3, 34.459, 0.86974, 8.333, 0.7250),
            (0.4, 32.093, 0.81066, 8.333, 0.7250),
            (0.5, 32.093, 0.64403, 8.333, 0.7250),
            (float('inf'), 33.504, 0.60486, 8.333, 0.7250)
        ],
        'E': [
            (0.4, 24.26, 0.83660, 6.25, 0.54287),
            (0.5, 23.331, 0.81956, 6.25, 0.54287),
            (float('inf'), 24.26, 0.83660, 6.25, 0.54287)
        ],
        'F': [
            (0.5, 15.209, 0.81558, 4.167, 0.34006),
            (float('inf'), 15.209, 0.81558, 4.167, 0.34006)
        ]
    }

    # Initialize arrays to store accumulated sigma values for each distance
    sigma_y_accum = np.zeros_like(total_dist)
    sigma_z_accum = np.zeros_like(total_dist)

    for stab_class in np.unique(stab_classes):
        class_sigma_y = np.zeros_like(total_dist)
        class_sigma_z = np.zeros_like(total_dist)

        for i, (limit, a, b, c, d) in enumerate(params[stab_class]):
            if i == 0:  # First limit applies to all distances up to that limit
                mask = (total_dist <= limit) 
            elif limit == float('inf'):  # Last limit applies to all distances above the previous limit
                mask = (total_dist > params[stab_class][i - 1][0]) 
            else:
                mask = (total_dist > params[stab_class][i - 1][0]) & (total_dist <= limit) 
            #print(mask, limit, total_dist)

            # Apply the mask and ensure no invalid values are passed to log and power functions
            valid_mask = mask & (total_dist > 0)  # Ensure total_dist > 0 to avoid log and power issues

            if np.any(valid_mask):
                big_theta = 0.017453293 * (c - d * np.log(total_dist[valid_mask]))  # Log is safe here
                class_sigma_y[valid_mask] = 465.11628 * total_dist[valid_mask] * np.tan(big_theta)
                class_sigma_z[valid_mask] = np.minimum(a * np.power(total_dist[valid_mask], b), 5000)

        # Accumulate computed values
        sigma_y_accum += class_sigma_y
        sigma_z_accum += class_sigma_z
        #print("Sigma Y Values for Class", stab_class, ":", class_sigma_y)

    # Divide by the number of stability classes to average
    sigma_y_vals = sigma_y_accum / len(np.unique(stab_classes))
    sigma_z_vals = sigma_z_accum / len(np.unique(stab_classes))

    return sigma_y_vals, sigma_z_vals



def plot_2d_plume_concentration(big_C, plot_title, output_dir):
    """
    Plot and save the averaged CH4 concentrations 
    """
    # Convert big_C unit from [ug m^-3] to [ppm]
    big_C = big_C / 714.6

    colors = ["white", "red"]  # White for low concentrations, red for high
    cmap = LinearSegmentedColormap.from_list("custom_red", colors, N=256)
    
    fig, axs = plt.subplots(1, 2, figsize=(14, 6)) 
    
    #Sum over the Z-axis (height dimension)
    C_summed_z = np.mean(big_C, axis=2)
    c1 = axs[0].imshow(C_summed_z.T, cmap=cmap, origin='lower', aspect='auto', vmin=0, vmax=10)
    axs[0].set_title(f'{plot_title}')
    axs[0].set_xlabel('Longitude Index (x)')
    axs[0].set_ylabel('Latitude Index (y)')
    cbar1 = fig.colorbar(c1, ax=axs[0])
    cbar1.set_label('Plume Concentration (ppm)')
    
    # um over the X-axis (longitude dimension)
    C_summed_x = np.mean(big_C, axis=0)
    c2 = axs[1].imshow(C_summed_x.T, cmap=cmap, origin='lower', aspect='auto', vmin=0, vmax=10)
    axs[1].set_title(f'{plot_title}')
    axs[1].set_xlabel('Latitude Index (y)')
    axs[1].set_ylabel('Height Index (z)')
    cbar2 = fig.colorbar(c2, ax=axs[1])
    cbar2.set_label('Plume Concentration (ppm)')
    
    # Adjust layout to avoid overlap
    plt.tight_layout()

    # Save the plot
    plt.savefig(os.path.join(output_dir, plot_title + ".png"))
    
    # Show the plot (remove if running in a non-interactive environment)
    # plt.show()

    # Close the figure to free memory
    plt.close(fig)




def gplume(Q, stab_class, x_p, y_p, z_p, x_r_vec, y_r_vec, z_r_vec, WS, WA_x, WA_y):
    """
    Calculate the CH4 concentrations using the Gaussian plume model.
    """

    # Shift coordinates according to wind direction
    x1 = x_r_vec - x_p 
    y1 = y_r_vec - y_p

    # Scalar product for angle calculation
    dot_product = x1* WA_x +  y1* WA_y
    magnitudes = WS * np.sqrt(x1**2 + y1**2)  # product of magnitude of vectors

    # Angle between wind and point (x, y)
    subtended = np.arccos(dot_product / (magnitudes + 1e-15))
    hypotenuse = np.sqrt(x1**2 + y1**2)  # distance to point x, y from stack

    # Distance along the wind direction to perpendicular line that intersects x, y
    downwind = x1 * WA_x + y1 * WA_y

    # Calculate crosswind distance
    crosswind = x1 * WA_y - y1 * WA_x
    
    #print("Debug downwind distance", downwind[-3:-1, -1, -1], downwind[-3:-1, -1, -1], downwind[-3:-1, -1, -1],downwind[-3:-1, -1, -1])

    # Compute sigma_y and sigma_z based on stability class and downwind distance
    sigma_y, sigma_z = compute_sigma_vals(stab_class, downwind / 1000)  # Assuming distance is in meters, convert to kilometers for sigma calculation
    
    #print("Debug gplume", sigma_y[-3:-1, -1, -1], sigma_z[-3:-1, -1, -1], downwind[-3:-1, -1, -1], stab_class)

    # Gaussian plume model calculation
    C = (Q / (2 * np.pi * WS * sigma_y * sigma_z)) * \
        np.exp(-0.5 * (crosswind ** 2) / sigma_y ** 2) * \
        (np.exp(-0.5 * (z_r_vec - z_p) ** 2 / sigma_z ** 2) + np.exp(-0.5 * (z_r_vec + z_p) ** 2 / sigma_z ** 2))
    
    # Convert from kg/m^3 to ug/m^3
    conversion_factor = 1e9 #  * 1.524
    C *= conversion_factor
    
    # Replace NAs with zeros
    C = np.where(np.isnan(C), 0, C)
    
    return C

def create_grid(dxy, dz, num_xygrid, num_zgrid):
    """
    Create a 3D grid using the given resolution
    """
    grid_x, grid_y, grid_z = np.meshgrid(
        np.arange(-dxy * num_xygrid/3, dxy * num_xygrid/3*2, dxy),  # X-axis grid
        np.arange(-dxy * num_xygrid/3, dxy * num_xygrid/3*2, dxy),  # Y-axis grid
        np.arange(0, dz * num_zgrid, dz),  # Z-axis grid
        indexing='ij'
    )
    return grid_x, grid_y, grid_z


def generate_wind_directions(WD_mean, spread, days):
    """
    Generate random wind directions based on a mean direction and spread.
    """
    wind_dir = WD_mean + spread * np.sqrt(2.) * erfcinv(2. * np.random.rand(24 * days, 1))
    return np.mod(wind_dir, 360)  # Ensure within 0-360 degrees


def compute_mean_concentration(grid_x, grid_y, grid_z, initial_Q, stack_x, stack_y, H, wind_speed, wind_dir, stab_class):
    """
    Compute the CH4 mean concentration over all wind directions
    """
    big_C = np.zeros((grid_x.shape[0], grid_y.shape[1], grid_z.shape[2], len(wind_dir)))

    for i in range(len(wind_dir)):
        wind_direction = wind_dir[i, 0]
        WA = np.radians(wind_direction)
        WA_x, WA_y = np.cos(WA), np.sin(WA)

        concentrations = gplume(initial_Q, stab_class, stack_x, stack_y, 
                                H, grid_x, grid_y, grid_z, wind_speed, WA_x, WA_y)
        big_C[:, :, :, i] += concentrations  # ug/m3 

    mean_conc = np.mean(big_C, axis=3)
    return mean_conc


def apply_noise(mean_conc, noise_level_multiplier=2, sigma=[1, 1, 1]):
    """
    Apply Gaussian noise to the mean concentration to simulate atmospheric turbulence
    """
    noise_level = mean_conc.mean() * noise_level_multiplier
    white_noise = np.random.normal(0, noise_level, mean_conc.shape)
    noise = gaussian_filter(white_noise, sigma=sigma)
    return mean_conc + noise


# Run the Gaussian plume simulation
def run_simulation(output_dir, dxy, dz, num_xygrid, num_zgrid, days, mean_wind_directions,
                   direction_spreads, day_or_night, incoming_solar_radiations, cloud_covers,
                   stack_x, stack_y, initial_Q, H, plotting_on=False):

    """
    Main function to compute CH4 concentrations using Gaussian Plume Model
    """
    grid_x, grid_y, grid_z = create_grid(dxy, dz, num_xygrid, num_zgrid)
    print("Grid dimensions:", grid_x.shape, grid_y.shape, grid_z.shape)

    for s in range(len(H)):
        for wind_speed in [1.5, 4, 6.5]:
            rough_steady_state_time = num_xygrid * dxy / np.sqrt(wind_speed)
            print(f"Steady state in ~{rough_steady_state_time} seconds")

            for WD_mean in mean_wind_directions:
                for spread in direction_spreads:
                    wind_dir = generate_wind_directions(WD_mean, spread, days)

                    for day_time in day_or_night:
                        for incoming_solar_radiation in incoming_solar_radiations:
                            for cloud_cover in cloud_covers:
                                stab_class = get_stab_class(wind_speed, day_time, 
                                                            solar_radiation=incoming_solar_radiation, 
                                                            cloud_cover_fraction=cloud_cover)

                                print(f"H_{H[s]}, WS_{wind_speed}, stability_{stab_class}, WD_{WD_mean}, spread_{spread}, radiaton_{incoming_solar_radiation}")

                                mean_conc = compute_mean_concentration(grid_x, grid_y, grid_z, initial_Q, 
                                                                       stack_x[s], stack_y[s], H[s], 
                                                                       wind_speed, wind_dir, stab_class)
               
                                mean_conc = apply_noise(mean_conc) # Overwrite the same variable to reduce RAM usage 

                                file_name = f"gplume_output_WS_{wind_speed}_stackH_{H[s]}_stability_{stab_class}_radiaton_{incoming_solar_radiation}_WD_{WD_mean}_spread_{spread}.npy"
                                np.save(os.path.join(output_dir, file_name), mean_conc)

                                if plotting_on:
                                    plot_2d_plume_concentration(mean_conc, 
                                        f"Noisy_WS_{wind_speed}_stackH_{H[s]}_radiaton_{incoming_solar_radiation}_stability_{stab_class}_WD_spread_{spread}", output_dir)

    print("Simulation completed.")