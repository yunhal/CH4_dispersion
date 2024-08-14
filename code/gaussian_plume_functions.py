__all__ = [
    'latlon_to_utm',
    'get_stab_class',
    'compute_sigma_vals',
    'gplume',
    'convert_utm_to_latlon_df',
    'plot_2d_plume_concentration'
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

def convert_utm_to_latlon_df(df, easting_col, northing_col, zone, northern=True):
    # Create a Transformer object for UTM to WGS84 conversion
    hemisphere = 'north' if northern else 'south'
    transformer = Transformer.from_proj(
        proj_from=f"+proj=utm +zone={zone}, +{hemisphere} +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
        proj_to="epsg:4326",  # WGS84
        always_xy=True
    )

    df['lon'], df['lat'] = transformer.transform(df[easting_col].to_numpy(), df[northing_col].to_numpy())
    return df


def latlon_to_utm(latitude, longitude):
    # Define the geographic coordinate system WGS 84
    wgs84 = CRS('EPSG:4326')
    
    # Determine the UTM zone
    zone_number = int((longitude + 180) // 6) + 1
    hemisphere = 'north' if latitude >= 0 else 'south'
    
    # Define the UTM CRS based on the zone and hemisphere
    utm_crs = CRS(proj='utm', zone=zone_number, datum='WGS84', south=(hemisphere == 'south'))
    project = Proj(utm_crs)
    easting, northing = project(longitude, latitude)
    
    zone_letter = 'N' if latitude >= 0 else 'S'

    return easting, northing # , int(zone_number), zone_letter


def get_stab_class(U, day_time, solar_radiation=None, cloud_cover_fraction=None):
    """
    Determine stability class based on wind speed, time of day, 
    incoming solar radiation during the daytime, and cloud cover fraction during nighttime.
    
    Parameters:
        U (float): Wind speed in m/s.
        day_time (bool): True if it's daytime, False if it's nighttime.
        solar_radiation (str): Incoming solar radiation during daytime ('strong', 'moderate', 'weak').
        cloud_cover_fraction (str): Cloud cover fraction during nighttime ('low', 'moderate', 'high').
        
    Returns:
        str: Stability class.
    """
    if day_time:
        if solar_radiation == 'strong':
            if U < 2:
                return "A"
            elif 2 <= U < 3:
                return "B"
            elif 3 <= U < 5:
                return "B"
            elif 5 <= U < 6:
                return "C"
            else:
                return "D"
        elif solar_radiation == 'moderate':
            if U < 2:
                return "B"
            elif 2 <= U < 3:
                return "C"
            elif 3 <= U < 5:
                return "C"
            elif 5 <= U < 6:
                return "D"
            else:
                return "D"
        elif solar_radiation == 'weak':
            if U < 2:
                return "C"
            elif 2 <= U < 3:
                return "C"
            elif 3 <= U < 5:
                return "D"
            elif 5 <= U < 6:
                return "D"
            else:
                return "D"
    else:
        if cloud_cover_fraction == 'low':
            if U < 2:
                return "F"
            elif 2 <= U < 3:
                return "F"
            elif 3 <= U < 5:
                return "E"
            elif 5 <= U < 6:
                return "D"
            else:
                return "D"
        elif cloud_cover_fraction == 'high':
            if U < 2:
                return "E"
            elif 2 <= U < 3:
                return "E"
            elif 3 <= U < 5:
                return "D"
            elif 5 <= U < 6:
                return "D"
            else:
                return "D"



def compute_sigma_vals(stab_classes, total_dist):
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
            if np.any(mask):
                big_theta = 0.017453293 * (c - d * np.log(total_dist[mask]))
                class_sigma_y[mask] = 465.11628 * total_dist[mask] * np.tan(big_theta)
                class_sigma_z[mask] = np.minimum(a * (total_dist[mask] ** b), 5000)

        # Accumulate computed values
        sigma_y_accum += class_sigma_y
        sigma_z_accum += class_sigma_z
        #print("Sigma Y Values for Class", stab_class, ":", class_sigma_y)

    # Divide by the number of stability classes to average
    sigma_y_vals = sigma_y_accum / len(np.unique(stab_classes))
    sigma_z_vals = sigma_z_accum / len(np.unique(stab_classes))

    return sigma_y_vals, sigma_z_vals



def plot_2d_plume_concentration(big_C, plot_title):
    colors = ["white", "red"]  # white for low concentrations, red for high
    cmap = LinearSegmentedColormap.from_list("custom_red", colors, N=256)
    
    fig, ax = plt.subplots()

    # Sum over the Z-axis (height dimension)
    C_summed = np.mean(big_C, axis=2)

    # Create the heatmap
    c = ax.imshow(C_summed.T, cmap=cmap, origin='lower', aspect='auto', vmin=0, vmax=1000)
    ax.set_title(plot_title)
    ax.set_xlabel('Longitude Index')
    ax.set_ylabel('Latitude Index')

    # Add a color bar
    cbar = fig.colorbar(c, ax=ax)
    cbar.set_label('Plume Concentration (ug m-3)')

    # Sum over the X (longitude dimension)
    C_summed = np.mean(big_C, axis=0)

    # Create the plot
    fig, ax = plt.subplots()

    # Plotting with y on the x-axis and z on the y-axis
    c = ax.imshow(C_summed.T, cmap=cmap, origin='lower', aspect='auto', vmin=0, vmax=1000)

    ax.set_title(plot_title)
    ax.set_xlabel('Latitude Index (y)')
    ax.set_ylabel('Height Index (z)')

    # Add a color bar
    cbar = fig.colorbar(c, ax=ax)
    cbar.set_label('Plume Concentration (ug m-3)')

    # Display the plot
    plt.show()



def gplume(Q, stab_class, x_p, y_p, z_p, x_r_vec, y_r_vec, z_r_vec, WS, WA_x, WA_y):
    """Calculate the contaminant concentration using the Gaussian plume model."""

    
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

