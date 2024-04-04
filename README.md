# CryoCube

## **Lexcube Visualization of DAS Records**

This repository is dedicated to building and visualizing interactive 3D Data Cubes from DAS (Distributed Acoustic Sensing) records.

To showcase the workflow, we use a 10 min long DAS record section from Rhonegletscher (Switzerland) from 2020 including 100 channels. We first fourier transform the time series records to compute spectrograms (1-100 Hz), which we organize as a 3D data cube with the axes time, frequency and DAS channel. The data cube is stored as cryo_subcube.zarr and visualized using lexcube.

<img align="center" src="https://github.com/JosepinaU/CryoCube/assets/36039541/e0b6f663-30ea-4a3f-b0c9-a98b7bcabafc">
<br/><br/><br/>
<img align="center" src="https://github.com/JosepinaU/CryoCube/assets/36039541/10a313c1-5df4-4403-aca2-e61d8a55e365">

## Features

- **Cube Creation**: Automated processing of DAS data including 1) Fast Fourier Transformation of the DAS time series data and 2) Data Cube generation of the DAS data in spectral domain (zarr-formatted). 
- **Cube Visualization**: Visualization of the spectral Data Cube using Lexcube (https://www.lexcube.org/) in a Python based Jupyter notebook.

## Getting Started

### Environment, Installation and Data

The workflow requires a set of libraries stored in `requirements.txt`, which need to be installed in your preferred environment (e.g. conda or python) 

```console
conda create --name lexcube
conda activate lexcube
conda install pip
pip install -r requirements.txt
```

OR

```console
python -m venv lexcube
source lexcube/bin/activate
pip install -r requirements.txt
```


Our test data set can be downloaded via: [testdata_rhone](https://cloud.scadsai.uni-leipzig.de/index.php/s/QMZSeXCfkqPFna7)

Unzip `./testdata_rhone.zip` to the project folder. The testdata features five *.hdf5 formatted files representing 30 s seconds of strain rate data for the same 2496 DAS channels each. The DAS cable covered the entire glacier, from accumuluation zone to ablation zone (~ 9 km). The channel spacing was set to 4 m.

### Cube Creation

To create the zarr-formatted spectral Data Cube from cryoseismological DAS records, execute:

```console
python3 create_cube.py
```

This script generates a spectral Data Cube from DAS data recorded on Rhonegletscher (Switzerland) in summer 2020.  

For all DAS channels the data is fourier transformed in the frequency band of 1-100 Hz with a resolution of d_f = 1 Hz and using a time window length = 1/d_f s (overlap = 0.5). The resulting Data Cube is a spectrogram for each DAS channel and it is stored as cryo_cube.zarr. To individually configure the cube creation adjust the key parameters documented in the script.

### Cube Visualization

To visualize the zarr-fomatted Data Cube (cryo_cube.zarr), execute the jupyter notebook: 

`lexcube_visualization.ipynb`

... and play. :)

## Get in touch with us

ScaDS.AI (Center for Scalable Data Analytics and Artificial Intelligence) Dresden/Leipzig is a center for Data Science, Artificial Intelligence and Big Data with locations in Dresden and Leipzig. It is one of the five new AI centers in Germany funded under the federal government’s AI strategy by the Federal Ministry of Education and Research and the Free State of Saxony. It is established as a permanent research facility at both locations with strong connections to the local universities: the TUD Dresden University of Technology and the Leipzig University. 
In the ScaDS.AI Earh & Environmental Sciene group (https://scads.ai/research/applied-ai-and-big-data/environment-and-earth-sciences/), we develop and apply multimodal and multiscale machine learning techniques and visualization tools that respect the characteristics and the heterogeneity of Earth System Science data. 

## Acknowledgements

The developers acknowledge the financial support by the Federal Ministry of Education and Research of Germany and by the Sächsische Staatsministerium für Wissenschaft Kultur und Tourismus in the program Center of Excellence for AI-research “Center for Scalable Data Analytics and Artificial Intelligence Dresden/Leipzig”, project identification number: ScaDS.AI.
