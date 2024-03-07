# CryoCube

![multiicequake](https://github.com/JosepinaU/CryoCube/assets/36039541/10a313c1-5df4-4403-aca2-e61d8a55e365)



**Lexcube Visualization of DAS Records**

This repository is dedicated to building and visualizing interactive 3D Data Cubes from DAS (Distributed Acoustic Sensing) records.

To showcase the workflow, we use a 10 min long DAS record section from Rhonegletscher (Switzerland) from 2020 including 100 channels. We first fourier transform the time series records to compute spectrograms (1-100 Hz), which we organize as a 3D data cube with the axes time, frequency and DAS channel. The data cube is stored as cryo_subcube.zarr and visualized using lexcube.


### Features

- **Cube Creation**: Automated processing of DAS data including 1) Fast Fourier Transformation of the DAS time series data and 2) Data Cube generation of the DAS data in spectral domain (zarr-formatted). 
- **Cube Visualization**: Visualization of the spectral Data Cube using Lexcube (https://www.lexcube.org/) in a Python based Jupyter notebook.

## Getting Started

### Installation and Environment

The workflow requires a set of libraries stored in `requirements.txt`, which need to be installed in your preferred environment (e.g. conda) 

### Cube Creation

To create the zarr-formatted spectral Data Cube from cryoseismological DAS records, execute:
`python create_cube.py 24`

This script generates a Data Cube for a couple of minutes of a specific day (24/07/2020) of DAS data recorded at the Rhonegletscher (Switzerland). The resulting Data Cube features a spectrogram for each DAS channel and stored as cryo_subcube.zarr.

### Cube Visualization

To visualize the zarr-fomatted Data Cube (cryo_subcube.zarr), execute the notebook: 
`lexcube_visualization.ipynb`

