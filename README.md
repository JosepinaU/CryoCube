# CryoCube
**Lexcube Visualization of DAS Records**

This repository is dedicated to processing and visualizing data cubes. The project involves the creation as well as visualization of a 3D data cube.

We use a 10 min long DAS record section from Rhonegletscher (Switzerland) from 2020 including 100 channels. We fourier transform the time series records to compute spectrograms (1-100 Hz), which we organize as a 3D data cube with the axes time, frequency and DAS channel. The data cube is stored as cryo_subcube.zarr and visualized using lexcube.


### Features

- **Cube Creation**: Automated processing of DAS data to generate a zarr-format data cube.
- **Cube Visualization**: Visualize the data cube using Lexcube in a Python notebook.

## Getting Started

### Cube Creation

To create a zarr data cube that processes 24/7 of the data, execute:

`python create_cube.py 24`

This script generates a data cube for a specific day of DAS data recorded at the Rhonegletscher. The resulting cube features a spectrogram for each DAS channel.

### Cube Visualization

To visualize the zarr data cube, execute the notebook: `lexcube_visualization.ipynb`

