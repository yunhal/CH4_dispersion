# Gaussian CH4 Modeling

**A Python-based modeling tool for simulating methane plume dispersion using Gaussian Puff or Plume models and calculating the transmittance of atmospheric methane using HITRAN data (only for plume model).**


## Introduction

Gaussian CH4 Modeling is a Python-based tool designed to simulate methane concentrations using Gaussian Puff or Plume models and calculate the transmittance of methane through the atmosphere. The tool uses HITRAN database to compute the light absorption by CH4 at a given wavelength band uses for satellite measurements. I intended to generate a synthetic CH4 plume, but the tool can  be used with geo coordinate data with minor modifications. 

## Features

- **Gaussian Plume Simulation:** Simulate methane plumes based on various atmospheric conditions and source parameters.
- **HITRAN Data Integration:** Fetch and utilize HITRAN absorption data to calculate methane transmittance.
- **Plume Visualization:** Optionally visualize methane plume dispersion.
- **Customizable Parameters:** Easily adjust grid parameters, wind direction, stability classes, source parameters, and more.
- **Multiprocessor Processing:** Utilize multiprocessors for the Gaussian Puff model to speed up the calculations.


## Project Structure

```plaintext
Gaussian_CH4_Modeling/
├── input_data/                    # Directory for storing input data (e.g., HITRAN data)

├── output_data/                   # Output directory for simulation results (not included in the gitrepo)
│   ├── gplume_conc/               # Output concentration fields
│   └── gplume_transmit/           # Calculated transmittance files

├── code/                          # Source code for the project
│   ├── main.py                    # Main script to run the Gaussian Plume model and compute CH4 transmittance using HITRAN database
│   ├── gaussian_plume_functions.py  # Functions related to Gaussian plume modeling
│   ├── ch4_optics_functions.py    # Functions for CH4 transmittance calculations
│   ├── hapi.py                    # HITRAN python script
│   ├── compare_ch4_absorption_data.ipynb     # Compare HITRAN-based absorption coef and the values obtained from satellite retrieval data
│   └── gaussian_puff_and_plume_models.ipynb  # Main script that runs Gaussian Puff and Plume models

├── requirements.txt               # Python dependencies  (not included yet)
└── README.md                      # Project README file

```

## License and Contact Information

This project is licensed under the MIT License - see the LICENSE file for details

Created by Yunha Lee. For any questions or suggestions, feel free to reach out or submit an issue on GitHub.