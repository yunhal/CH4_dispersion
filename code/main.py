import os
import numpy as np
import ssl
from hapi import *
from gaussian_plume_functions import *
from ch4_optics_functions import *

def main():

    # Configuration parameters for plume simulation
    output_dir = "../output_data/gplume_conc/"

    # grid parameters
    dxy = 20
    dz = 20
    num_xygrid = 250 # 2.5 km x 2.5 km
    num_zgrid = 100 # 2 km (to cover the day-time planetary boundary layer, which goes up to 1-2 km)

    # wind direction parameters
    days = 5 # only affect the wind directions data generation
    mean_wind_directions = [45]
    direction_spreads = [15, 30, 45]

    # stability parameters
    day_or_night = [True]
    incoming_solar_radiations = ["strong", "moderate", "weak"]
    cloud_covers = ["None"]

    # CH4 source parameters
    initial_Q = 1  # kg/s (later, it will be linearly scaled to Q_targets)
    stack_x = [0., 0]  # set all source facility at x = 0 (it doesn't matter where is located for our purpose)
    stack_y = [0., 0] # set all source facility at y = 0 (it doesn't matter where is located for our purpose)
    H = [30] # stack height (varying heights by each facility)

    # turning on and off for plotting
    plotting_on = True

    # Run the Gaussian plume simulation
    run_simulation(output_dir, dxy, dz, num_xygrid, num_zgrid, days, mean_wind_directions,
                    direction_spreads, day_or_night, incoming_solar_radiations, cloud_covers,
                    stack_x, stack_y, initial_Q, H, plotting_on)

    # Configuration for accessing the HITRAN database
    os.environ['GLOBAL_HOST'] = 'http://hitran.org'
    ssl._create_default_https_context = ssl._create_unverified_context

    # Constants for transmittance calculation
    satellites = ["sentinel2"]  # List of satellites
    hitran_data_dir = "/Users/yunhalee/Documents/methanDart/Gaussian_CH4_Modeling/input_data/hitran"
    output_dir = "../output_data/gplume_transmit/"
    conc_dir = "../output_data/gplume_conc/"

    # Fetch HITRAN data for each satellite band with CH4 absorption
    for satellite in satellites:
        fetch_hitran_data(satellite, hitran_data_dir)

    # Process concentration files and calculate/save transmittance
    Q_targets = [10000]  # kg/hr
    process_concentration_files(conc_dir, Q_targets, satellites, hitran_data_dir, output_dir, dz)

    print("Processing completed.")

if __name__ == "__main__":
    main()