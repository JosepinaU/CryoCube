# CryoCube
Lexcube Visualization of DAS records

We use a 10 min long DAS record section from Rhonegletscher (Switzerland) from 2020 including 100 channels. We fourier transform the time series records to compute spectrograms (1-100 Hz), which we organize as a 3D data cube with the axes time, frequency and DAS channel. The data cube is stored as *.zarr and visualized using lexcube.
