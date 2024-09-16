__all__ = [
    'fetch_hitran_data',
    'calculate_transmittance_for_concentration',
    'process_concentration_files',
]

import os
import ssl
import numpy as np
import pandas as pd
import h5py
from hapi import *
import re


def fetch_hitran_data(satellite, hitran_data_dir):
    """
    Fetch HITRAN data for the given satellite.
    """
    obs_file = f"../input_data/{satellite}_band_wavelength.csv"
    bands = pd.read_csv(obs_file)

    # Subset the data based on CH4_AC > 0
    subset_bands = bands[bands["CH4_AC"] < 0]
    subset_bands = subset_bands[["Band", "Wavelength_start", "Wavelength_end", "CH4_AC"]]

    for index, row in subset_bands.iterrows():
        wavelength_start = row["Wavelength_start"]
        wavelength_end = row["Wavelength_end"]
        ch4_ac = row["CH4_AC"]
        band = row["Band"]

        # Compute wavelength range and wavenumber in cm^-1
        low_wavenumbers = 1.e7 / wavelength_end
        high_wavenumbers = 1.e7 / wavelength_start

        # Set the local folder to download the HITRAN data
        db_begin(hitran_data_dir)
        abs_name = f"CH4_{satellite}_{band}"
        fetch(abs_name, 6, 1, low_wavenumbers, high_wavenumbers)

def calculate_transmittance_for_concentration(mean_conc, Q_target, satellite, hitran_data_dir, filename, output_dir, dz):
    """
    Calculate the transmittance for a given concentration field
    """
    mean_conc_scaled = mean_conc / 1e9 * (Q_target / 3600)  # Convert to kg/m^3

    for hitran_file in os.listdir(hitran_data_dir):
        if hitran_file.endswith(".data"):
            pattern = re.compile(rf"CH4_{satellite}_(?P<band>.*?)\.data")
            match = pattern.match(hitran_file)

            if match:
                band = match.group("band")
                print(f"Loading HITRAN file for Satellite: {satellite}, Band: {band}")
            else:
                continue

            # Check if the transmittance data already exists
            file_name = f"transmittance_Q_{Q_target}_satellite_{satellite}_band_{band}_{filename}.h5"
            if os.path.isfile(os.path.join(output_dir, file_name)):
                print(f"File exists: {file_name}")
                continue

            # Load HITRAN data and compute absorption coefficients
            file_base_name, _ = os.path.splitext(hitran_file)
            nu, abs_coef = absorptionCoefficient_Lorentz(SourceTables=file_base_name,
                                                         Diluent={'air': 1.0},
                                                         HITRAN_units=True, Environment={'T': 296., 'p': 1.},
                                                         WavenumberStep=0.1)
            # Convert cm^2/molecule to m^2/kg
            ch4_molar_weight = 0.01604  # kg/mole
            Avogadro = 6.022e23  # molecules per mole
            convert_unit = Avogadro / ch4_molar_weight / 1e4
            abs_coef = abs_coef.astype(np.float32) * convert_unit

            # Mean of the absorption coefficient
            absorption_mean = abs_coef.mean()

            # Calculate extinction and transmittance
            total_extinction = np.sum(mean_conc_scaled[:, :, :] * dz * absorption_mean, axis=2, dtype=np.float32)

            # Calculate total transmittance based on Beer-Lambert Law
            # final_transmittance = np.exp(-np.mean(total_extinction, axis=-1, dtype=np.float32), dtype=np.float32)
            final_transmittance = np.exp(-total_extinction, dtype=np.float32)

            # Save output
            with h5py.File(os.path.join(output_dir, file_name), 'w') as hf:
                hf.create_dataset('CH4 transmittance', data=final_transmittance)

def scale_transmittance(final_transmittance, Q_target_old, Q_target_new):
    scaling_factor = Q_target_new / Q_target_old
    
    # Scale the transmittance
    final_transmittance_new = final_transmittance ** scaling_factor
    
    return final_transmittance_new

def process_concentration_files(conc_dir, Q_targets, satellites, hitran_data_dir, output_dir, dz):
    """
    Process all concentration files and calculate transmittance for each file
    """
    for filename in os.listdir(conc_dir):
        if filename.endswith('.npy'):  
            # Read concentration fields
            mean_conc = np.load(os.path.join(conc_dir, filename), allow_pickle=True)
            print(f"Reading concentrations from {filename}")

            # Convert mean_conc to float32 to reduce memory usage
            mean_conc = mean_conc.astype(np.float32)

            for Q_target in Q_targets:
                for satellite in satellites:
                    calculate_transmittance_for_concentration(mean_conc, Q_target, satellite, hitran_data_dir, filename, output_dir, dz)
        else:
            print(f"Skipping non-npy file: {filename}")