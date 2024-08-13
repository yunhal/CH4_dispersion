__all__ = [
    'latlon_to_utm',
    'is_day',
    'get_stab_class',
    'compute_sigma_vals',
    'gpuff',
    'gplume',
    'process_chunk',
    'interpolate_wind_data',
    'plot_2d_concentration',
    'plot_2d_concentration_from_df',
    'plot_3d_concentration', 
    'plot_3d_concentration_with_slider',
    'calculate_end_time',
    'average_time_resolution',
    'convert_utm_to_latlon_df',
    'simulate_puff_concentration',
    'simulate_plume_concentration',
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
# This doesn't work somehow


    # Create a Transformer object for UTM to WGS84 conversion
    hemisphere = 'north' if northern else 'south'
    transformer = Transformer.from_proj(
        proj_from=f"+proj=utm +zone={zone}, +{hemisphere} +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
        proj_to="epsg:4326",  # WGS84
        always_xy=True
    )

    df['lon'], df['lat'] = transformer.transform(df[easting_col].to_numpy(), df[northing_col].to_numpy())
    return df

def calculate_end_time(start_time, chunk_size_minutes, n_chunks):

    # Calculate the total duration to be added to the start time based on n_chunks and chunk_size
    total_duration = timedelta(minutes=chunk_size_minutes * n_chunks)
    end_time = start_time + total_duration
    
    return end_time

def interpolate_wind_data(wind_data, dt):
    interpolated = [wind_data[0]]
    for i in range(1, len(wind_data)):
        step = np.linspace(wind_data[i-1], wind_data[i], num=int(60/dt))
        interpolated.extend(step[:-1])
    return np.array(interpolated)

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

def is_day(time):
    """Check if the time is during the day, considered as between 7 AM and 6 PM."""
    time_hour = time.hour
    return 7 <= time_hour <= 18

def get_stab_class(U, time):
    """Determine stability class based on wind speed and time of day."""
    if U < 2:
        return ["A", "B"] if is_day(time) else ["E", "F"]
    elif 2 <= U < 3:
        return ["B"] if is_day(time) else ["E", "F"]
    elif 3 <= U < 5:
        return ["B", "C"] if is_day(time) else ["D", "E"]
    elif 5 <= U < 6:
        return ["C", "D"] if is_day(time) else ["D"]
    else:
        return ["D"] if is_day(time) else ["D"]


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


def plot_3d_concentration(big_C, times):
    time_steps = big_C.shape[0]

    colors = ["white", "red"]  # White to red
    cmap = LinearSegmentedColormap.from_list("custom_red", colors, N=256)

    for t in range(time_steps):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Get the concentration data for this time step
        C_t = big_C[t]

        # Generate grid and flatten data
        X, Y, Z = np.meshgrid(np.arange(C_t.shape[0]), np.arange(C_t.shape[1]), np.arange(C_t.shape[2]), indexing='ij')
        X = X.flatten()
        Y = Y.flatten()
        Z = Z.flatten()
        C_t_flat = C_t.flatten()
        
        # Filter points where concentration is greater than zero
        mask = C_t_flat > 0
        X, Y, Z, C_t_flat = X[mask], Y[mask], Z[mask], C_t_flat[mask]

        # Sort points by Z for better depth visibility
        sorted_indices = np.argsort(Z)
        X, Y, Z, C_t_flat = X[sorted_indices], Y[sorted_indices], Z[sorted_indices], C_t_flat[sorted_indices]

        # Create the scatter plot with transparency
        sc = ax.scatter(X, Y, Z, c=C_t_flat, cmap=cmap, marker='o', alpha=0.5)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.title(f'Concentration at time {times[t]}')

        # Add color bar which maps values to colors.
        cbar = plt.colorbar(sc, ax=ax, format='%.1e')
        cbar.set_label('Concentration (ppm)')

        plt.savefig(f'3d_concentration_time_{t}.png')
        plt.show()
        plt.close()

def plot_2d_concentration(big_C, times):
    time_steps = big_C.shape[0]
    
    colors = ["white", "red"]  # white for low concentrations, red for high
    cmap = LinearSegmentedColormap.from_list("custom_red", colors, N=256)
    
    for t in range(time_steps):
        if t % 300 == 0:
            print(f"Plotting time step {t}")
            fig, ax = plt.subplots()

            # Sum over the Z-axis
            C_t_summed = np.sum(big_C[t], axis=2) 

            # Create the heatmap
            c = ax.imshow(C_t_summed.T, cmap=cmap, origin='lower', aspect='auto', vmin=0)
            ax.set_title(f'Concentration at time {times[t]}')
            ax.set_xlabel('Longitude Index')
            ax.set_ylabel('Latitude Index')

            # Add a color bar
            cbar = fig.colorbar(c, ax=ax)
            cbar.set_label('BigC Concentration (units)')

            # Display the plot
            plt.show()

            # # Extract the center slice along the first dimension 
            # center_pt = big_C[t].shape[0] // 2
            # center_slice = big_C[center_pt,:, :, t]

            # print(center_slice.shape)  # Should print (20, 10)

            # # Create the plot
            # fig, ax = plt.subplots()

            # # Plotting with y on the x-axis and z on the y-axis
            # c = ax.imshow(center_slice.T, cmap=cmap, origin='lower', aspect='auto', vmin=0)

            # ax.set_title('Concentration at Source Location')
            # ax.set_xlabel('Latitude Index (y)')
            # ax.set_ylabel('Height Index (z)')

            # # Add a color bar
            # cbar = fig.colorbar(c, ax=ax)
            # cbar.set_label('Plume Concentration (ug m-3)')

            # # Display the plot
            # plt.show()

def plot_2d_concentration_from_df(df):
    
    colors = ["white", "red"]  # white for low concentrations, red for high
    cmap = LinearSegmentedColormap.from_list("custom_red", colors, N=256)

    # subset key columns
    df = df [["times", "lat", "lon", "height", "Value"]]

    # Calculate mean values
    df_mean = df.groupby(['times', 'lat', 'lon'])['Value'].sum().reset_index()

    times = df_mean['times'].unique()
    
    for t in times:
        fig, ax = plt.subplots()

        # Filter data for the current time step
        current_data = df_mean[df_mean['times'] == t]

        # Pivot the data to create a 2D grid for plotting
        pivot_table = current_data.pivot_table(index='lat', columns='lon', values='Value', fill_value=0)

        # Create a meshgrid for pcolormesh
        grid_x, grid_y = np.meshgrid(pivot_table.columns, pivot_table.index)

        # Create the heatmap using pcolormesh
        heatmap = ax.pcolormesh(grid_x, grid_y, pivot_table, cmap=cmap, shading='auto', vmin=0) #, vmin=0, vmax=10)
        ax.set_title(f'Concentration at {pd.to_datetime(t)}')
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')

        cbar = fig.colorbar(heatmap, ax=ax)
        #cbar.set_ticks(np.linspace(0, 10, 11)) 
        cbar.set_label('Averaged Concentration (units)')

        plt.show()

import plotly.graph_objects as go

def plot_3d_concentration_with_slider(big_C, times):

    fig = go.Figure()

    grid_shape = big_C[0].shape
    X, Y, Z = np.meshgrid(np.arange(grid_shape[0]), np.arange(grid_shape[1]), np.arange(grid_shape[2]), indexing='ij')
    X = X.flatten()
    Y = Y.flatten()
    Z = Z.flatten()

    frames = []
    for t, C_t in enumerate(big_C):
        C_t_flat = C_t.flatten()

        # Filter out zero concentration values
        mask = C_t_flat > 0
        filtered_X = X[mask]
        filtered_Y = Y[mask]
        filtered_Z = Z[mask]
        filtered_C_t_flat = C_t_flat[mask]

        # Sort by Z for better visualization
        sorted_indices = np.argsort(filtered_Z)
        frames.append(go.Frame(
            data=[go.Scatter3d(
                x=filtered_X[sorted_indices],
                y=filtered_Y[sorted_indices],
                z=filtered_Z[sorted_indices],
                mode='markers',
                marker=dict(
                    size=2,
                    color=filtered_C_t_flat[sorted_indices],
                    colorscale='Reds',
                    opacity=0.6
                )
            )],
            name=str(times[t])
        ))

    # Set the initial state
    fig.add_trace(go.Scatter3d(
        x=frames[0].data[0]['x'],
        y=frames[0].data[0]['y'],
        z=frames[0].data[0]['z'],
        mode='markers',
        marker=dict(
            size=2,
            color=frames[0].data[0].marker.color,
            colorscale='Reds',
            opacity=0.6
        )
    ))

    # Add play and pause buttons
    fig.update_layout(
        updatemenus=[
            {
                "buttons": [
                    {
                        "args": [None, {"frame": {"duration": 500, "redraw": True}, "fromcurrent": True, "transition": {"duration": 300}}],
                        "label": "Play",
                        "method": "animate"
                    },
                    {
                        "args": [[None], {"frame": {"duration": 0, "redraw": False}, "transition": {"duration": 0}}],
                        "label": "Pause",
                        "method": "animate"
                    }
                ],
                "direction": "left",
                "pad": {"r": 10, "t": 87},
                "showactive": False,
                "type": "buttons",
                "x": 0.1,
                "xanchor": "right",
                "y": 0,
                "yanchor": "top"
            }
        ],
        sliders=[{
            "active": 0,
            "yanchor": "top",
            "xanchor": "left",
            "currentvalue": {
                "font": {"size": 20},
                "prefix": "Time:",
                "visible": True,
                "xanchor": "right"
            },
            "transition": {"duration": 300, "easing": "cubic-in-out"},
            "pad": {"b": 10, "t": 50},
            "len": 0.9,
            "x": 0.1,
            "y": 0,
            "steps": [
                {
                    "args": [[f.name], {"frame": {"duration": 500, "redraw": True}, "mode": "immediate", "transition": {"duration": 300}}],
                    "label": str(time),
                    "method": "animate"
                } for f, time in zip(frames, times)]
        }],
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z'
        )
    )

    fig.frames = frames

    return fig

def plot_2d_plume_concentration(big_C):
    colors = ["white", "red"]  # white for low concentrations, red for high
    cmap = LinearSegmentedColormap.from_list("custom_red", colors, N=256)
    
    fig, ax = plt.subplots()

    # Sum over the Z-axis (height dimension)
    C_summed = np.mean(big_C, axis=2)

    # Create the heatmap
    c = ax.imshow(C_summed.T, cmap=cmap, origin='lower', aspect='auto', vmin=0, vmax=1000)
    ax.set_title('Concentration')
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

    ax.set_title('Concentration')
    ax.set_xlabel('Latitude Index (y)')
    ax.set_ylabel('Height Index (z)')

    # Add a color bar
    cbar = fig.colorbar(c, ax=ax)
    cbar.set_label('Plume Concentration (ug m-3)')

    # Display the plot
    plt.show()

def average_time_resolution(big_C, times, output_dt_sec, output_dir,  grid_x, grid_y, grid_z):

    time_steps, grid_x_size, grid_y_size, grid_z_size = big_C.shape
    ("big_C shape", big_C.shape)
    
    # Flatten the array
    flat_big_C = big_C.ravel()

    # Get the multi-dimensional indices
    indices = np.unravel_index(np.arange(flat_big_C.size), big_C.shape)
    #print(indices[0])

    # Create a DataFrame based on known dimensions
    df = pd.DataFrame({
        'Value': flat_big_C,
        'time_index': indices[0],
        'lon_index': indices[1],
        'lat_index': indices[2],
        'height_index': indices[3]
    })

    # Fill the dataframe with actual grid values based on the 'grid_index'
    df['lon'] = df['lon_index'].apply(lambda x: grid_x[x, 0, 0])  
    df['lat'] = df['lat_index'].apply(lambda x: grid_y[0, x, 0])
    df['height'] = df['height_index'].apply(lambda x: grid_z[0, 0, x])
    df['times'] = df['time_index'].apply(lambda x: times[x])
    df['times'] = pd.to_datetime(df['times'])

    # Subset df 
    df = df [['lon', 'lat', 'height','Value','times']]

    df.set_index('times', inplace=True)

    print("df shape ", df.shape)

    # Resample and average only the 'Value' column while keeping the grid information
    resampled_df = df.groupby(
        [pd.Grouper(freq=f'{output_dt_sec}s'), 'lon', 'lat', 'height']
    ).mean().reset_index()

    # The 'time_index' no longer exists in the resampled_df since it's grouped by the new time intervals
    resampled_df['times'] = resampled_df['times'].dt.round(f'{output_dt_sec}s')

    print("final resampled_df shape", resampled_df.shape)

    output_filename = f'{output_dir}simulation_output_resampled_at_{output_dt_sec}sec.csv'
    resampled_df.to_csv(output_filename, index=False)

    return resampled_df

def simulate_puff_concentration(source_x, source_y, source_z, grid_x, grid_y, grid_z,  data, dt, emission_rate, times, chunk_size, n_chunks):

    # Wind data and interpolation
    WS_x = data['u_wind']
    WS_y = data['v_wind']
    WA = np.radians(data['wind_direction'])
    WA_x = np.cos(WA)
    WA_y = np.sin(WA)
    WS = np.sqrt(WS_x**2 + WS_y**2)

    # Emission rates
    Q_truth = np.full((len(WS)), emission_rate)

    # Setup chunk processing
    n_ints = len(WS) # ignore the last time 

    args = [(h, chunk_size, dt, n_ints, source_x, source_y, source_z, WS_x, WS_y, WS, Q_truth, grid_x, grid_y, grid_z, times) 
            for h in range(n_chunks)]

    print(f"chunk_size: {chunk_size}, n_chunks: {n_chunks}, dt: {dt}, n_ints: {n_ints}, "
      f"source_x: {source_x}, source_y: {source_y}, source_z: {source_z}, "
      f"WS_x length: {len(WS_x)}, WS length: {len(WS)}, Q_truth: {Q_truth[0]}, ")

    with Pool(processes=n_chunks) as pool:
        results = pool.map(process_chunk, args)
        #print("result", results)

    big_C = np.vstack(results)
    return big_C

def gpuff(Q, stab_class, x_p, y_p, x_r_vec, y_r_vec, z_r_vec, total_dist, H):

    """Calculate the contaminant concentration using the Gaussian puff model."""
    total_dist_km = total_dist / 1000
    sigma_y, sigma_z = compute_sigma_vals(stab_class, total_dist_km)

    print("Debug gplume", sigma_y[-3:-1, -1, -1], sigma_z[-3:-1, -1, -1], total_dist[-3:-1, -1, -1], stab_class)

    #print("debug sigma",stab_class, sigma_y[10,10,10])

    # Gaussian puff model calculation (Total reflection at z =0)
    C = (Q / ((2 * np.pi) ** 1.5 * sigma_y ** 2 * sigma_z)) * \
        np.exp(-0.5 * ((x_r_vec - x_p) ** 2 + (y_r_vec - y_p) ** 2) / sigma_y ** 2) * \
        (np.exp(-0.5 * (z_r_vec - H) ** 2 / sigma_z ** 2) + np.exp(-0.5 * (z_r_vec + H) ** 2 / sigma_z ** 2))

    # Convert from kg/m^3 to ug/m^3
    conversion_factor = 1e9 #  * 1.524
    C *= conversion_factor

    # Replace NAs (from zero total distance cases) with zeros
    C = np.where(np.isnan(C), 0, C)

    return C


def process_chunk(args):
    h, chunk_size, dt, n_ints, source_x, source_y, source_z, WS_x, WS_y, WS, Q_truth, grid_x, grid_y, grid_z, times = args

    chunk_start = int(h * chunk_size * 60 / dt)
    chunk_end = int(min((h + 1) * chunk_size * 60 / dt, n_ints))
    chunk_concentrations = np.zeros((chunk_end - chunk_start, grid_x.shape[0], grid_y.shape[1], grid_z.shape[2]))

    print("check chunk_concentration dimensions", chunk_concentrations.shape, chunk_start,chunk_end)

    for j in range(chunk_start, chunk_end):
        current_time = times[0] + timedelta(seconds=(chunk_start + j) * dt)
        stab_class = get_stab_class(WS.iloc[j], current_time)
        #print("debug process_chunk", WS[j], current_time, stab_class)

        for k in range(1, min(j - chunk_start, 300)+1):
            puff_x = source_x + np.sum(WS_x[j-k:j] * dt)
            puff_y = source_y + np.sum(WS_y[j-k:j] * dt)
            total_dist = np.sqrt((grid_x - puff_x)**2 + (grid_y - puff_y)**2)

            #print(f"in process_chunk, total_dist{total_dist[:2,:2, 0]},puff_x {puff_x}, source_x {source_x}, grid_x{grid_x[:1,:1,0]}, index{j, k, j}, wind{WS_x[j-k:j]},  Q{Q_truth[j]}" )

            # Compute concentrations using the vectorized gpuff function
            concentrations = gpuff(Q_truth[j], stab_class, puff_x, puff_y, grid_x, grid_y, grid_z, total_dist, source_z)
            chunk_concentrations[j - chunk_start] += concentrations

    return chunk_concentrations


def gplume_old(Q, stab_class, x_p, y_p, z_p, x_r_vec, y_r_vec, z_r_vec, WS, WA_x, WA_y):
    """Calculate the contaminant concentration using the Gaussian plume model."""

    
    # Calculate the total distance in the downwind direction
    total_dist = np.sqrt((x_r_vec - x_p)**2 + (y_r_vec - y_p)**2)
    
    # Compute sigma_y and sigma_z based on stability class and distance
    sigma_y, sigma_z = compute_sigma_vals(stab_class, total_dist/1000)
    
    print("Debug gplume", sigma_y[-3:-1, -1, -1], sigma_z[-3:-1, -1, -1], total_dist[-3:-1, -1, -1], stab_class)

    # Gaussian plume model calculation
    C = (Q / (2 * np.pi * WS * sigma_y * sigma_z)) * \
        np.exp(-0.5 * ((x_r_vec - x_p) ** 2 + (y_r_vec - y_p) ** 2) / sigma_y ** 2) * \
        (np.exp(-0.5 * (z_r_vec - z_p) ** 2 / sigma_z ** 2) + np.exp(-0.5 * (z_r_vec + z_p) ** 2 / sigma_z ** 2))
    
    # Convert from kg/m^3 to ug/m^3
    conversion_factor = 1e9 #  * 1.524
    C *= conversion_factor
    
    # Replace NAs with zeros
    C = np.where(np.isnan(C), 0, C)
    
    return C

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
    
    print("Debug downwind distance", downwind[-3:-1, -1, -1], downwind[-3:-1, -1, -1], downwind[-3:-1, -1, -1],downwind[-3:-1, -1, -1])

    # Compute sigma_y and sigma_z based on stability class and downwind distance
    sigma_y, sigma_z = compute_sigma_vals(stab_class, downwind / 1000)  # Assuming distance is in meters, convert to kilometers for sigma calculation
    
    print("Debug gplume", sigma_y[-3:-1, -1, -1], sigma_z[-3:-1, -1, -1], downwind[-3:-1, -1, -1], stab_class)

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

def simulate_plume_concentration(source_x, source_y, source_z, grid_x, grid_y, grid_z, data, emission_rate, times):
    #  Assuming single, average wind speed for simplicity
    WS_x = data['u_wind'].mean() 
    WS_y = data['v_wind'].mean()
    WS = np.sqrt(WS_x**2 + WS_y**2)
    WA_x = WS_x / WS
    WA_y = WS_y / WS



    # Emission rates
    Q_truth = emission_rate
    
    # Compute stability class for the entire simulation (Assume constant stability class for simplicity)
    stab_class = get_stab_class(WS, times[0]) 

    print("Debug plume", WS, stab_class)
    
    # Compute concentrations using the Gaussian Plume model
    concentrations = gplume(Q_truth, stab_class, source_x, source_y, source_z, grid_x, grid_y, grid_z, WS, WA_x, WA_y)
    
    print("Debug plume 2", concentrations.shape)

    return concentrations
