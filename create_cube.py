import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import scipy
from datetime import datetime, time, timedelta
import zarr
import numcodecs
import time
import os
from os import listdir
from os.path import isfile, join
import re
from operator import itemgetter
from multiprocessing import Pool
import multiprocessing
from functools import partial

"""
This script is used to create a data cube (*.zarr) containing 
spectrograms for measurements of strain rate recorded by 
DAS (distributed acoustic sensing). The data cube features 
the dimensions DAS channel, time, and frequency.

Please note that the script is designed for the structure 
of the example DAS data set (testdata_rhone)[https://cloud.scadsai.
uni-leipzig.de/index.php/apps/files/files/5700666?dir=/
environment-earth/Shares/CryoCube_GIT/testdata_rhone], which is 
stored as *.hdf5 files, each holding 30 seconds of 
continuous strain rate data (rows are DAS channels, 
columns are data samples). Depending on the input data, 
you may need to adapt the workflow. 

It is recommended to run this script with as many CPU cores as 
possible (and a significant amount of memory).

How it's used in the terminal:

create_cube.py 

If you have any questions, feel free to reach out to
to us!
"""

# Key adjustable parameters:

# d_f: Frequency resolution in Hz - the desired frequency granularity of the spectrogram.
d_f = 1

# loc_a and loc_e: Specify the cable section to be processed (in meters).
# 0 signifies the start of the cable.
loc_a, loc_e = 0, 9200

# DATA_PATH: Path to the directory where the h5 data files are stored.
DATA_PATH = "testdata_rhone/"

# nFiles: Determines how many h5 files are processed. For example, setting nFiles=10
# processes the first ten chronologically sorted h5 files in the DATA_PATH directory.
nFiles = 5

# nCores: Specifies how many CPU cores shall be used for computation. Adjust according
# to your system's capabilities and the workload size. The maximum number can be set to
# multiprocessing.cpu_count() to use all available cores.
nCores = 8

# Additional parameters:

# file_length: Length of a single h5 file in seconds.
file_length = 30

# nu: Sampling frequency in Hz of the recorded data.
nu = 1000

# freq_max: Maximum frequency to be considered in the analysis. All values above
# this frequency will be cut off.
freq_max = 100

#path and name of resulting zarr-formatted data cube.
ZARR_NAME = "cryo_cube.zarr"


def getTimeFromFilename(filename, day, month):
    """
    Extracts the time from a DAS-h5-file's filename and calculates 
    the time elapsed since the start of the day.
    
    Args:
        filename (str): The filename to extract time information from.
        day (int): The day to calculate the time difference for.
        month (int): The month to calculate the time difference for.
    
    Returns:
        float: The time difference in seconds between the file's timestamp and the start of the specified day.
    """
    # Extracting hours, minutes, and seconds from the filename
    hours_minutes_seconds = re.findall("_......\\.", filename)[0][1:]
    #print(hours_minutes_seconds)
    hour = int(hours_minutes_seconds[0:2])
    minutes= int(hours_minutes_seconds[2:4])
    seconds= int(hours_minutes_seconds[4:6])
    return (datetime(2020, month, day, hour, minutes, seconds) - datetime(2020, month, day, 0, 0, 0)).total_seconds()


def getFilenames(day, month, mypath):
    """
    Collects the filenames in the data folder and sorts them by time.
    
    Args:
        day (int): The day of interest.
        month (int): The month of interest.
        mypath (str): Path to the folder containing DAS data.
    
    Returns:
        dict: A dictionary where keys are integers and values are filenames, sorted chronologically.
    """
    # Filtering and sorting files based on their timestamps
    two_digit_month = str(month) if month >9 else "0"+str(month)
    #mypath = "/work/users/fu591lmui/RhoneData/2020"+two_digit_month+str(day)+"/"

    files = list([f for f in listdir(mypath) if isfile(join(mypath, f))])
    q = list([(mypath+file, getTimeFromFilename(file, day, month)) for file in files])
    q = sorted(q, key = itemgetter(1))
    return  {i:qq[0] for i,qq in enumerate(q)}


def channel_fourier(data, nseg, seg_len, taper, ind_a, ind_e, ind_f):
    """
    Applies Fourier Transformation to segments of DAS records to compute spectrograms.
    
    Args:
        data (ndarray): The raw data from DAS channels.
        nseg (int): Number of segments.
        seg_len (int): Length of each segment.
        taper (ndarray): The taper function to apply before the Fourier transform.
        ind_a (int): Start index of the cable section.
        ind_e (int): End index of the cable section.
        ind_f (int): Number of frequency bins.
    
    Returns:
        ndarray: A 3D array containing the Fourier transform for each segment and channel.
    """
    # Preparing data segments and applying Fourier transform
    
    #dividing the data into segments each consisting of desired amount of data points
    segs = ([data[i*(seg_len//2):(i+2)*(seg_len//2)] for i in range(nseg)])
    
    #transposing the segments individually to gain time series for each channel
    segs = [segs[i].T[ind_a:ind_e] for i in range(nseg)]
    

    # the first loop iterates over all segments (each corresponding to a time point)
    # in the second loop, the fourier transform gets applied on each channel
    
    Fsegs=np.zeros((nseg, ind_e-ind_a, ind_f))
    
    for j in range(nseg):
        for channel_number, channel  in enumerate(segs[j]):
            
            # note that modified_log(x)=10*log(x) (conversion to 
            
            fourier_transformed = np.fft.rfft(taper*channel, n=seg_len)
            fourier_transformed = ((10*np.log(np.abs(fourier_transformed)))**2)[0:ind_f]
            fourier_transformed[0]=0
            Fsegs[j][channel_number]=fourier_transformed
    
    return Fsegs


def create_spectro_segment(file_index, d_f, nseg, seg_len, ind_a, ind_e, ind_f, q):
    """
    Creates a spectrogram segment for a given file index.
    
    Args:
        file_index (int): The index of the file to process.
        d_f, nseg, seg_len, ind_a, ind_e, ind_f: Parameters for spectrogram computation (see above).
        q (dict): Dictionary of filenames to process.
    
    Returns:
        ndarray: The spectrogram segment for the specified file.
    """
    # Processing a single file to generate its spectrogram segment
    global nu, file_length, nFiles
    #opening the DAS-data
    f = h5py.File(q[file_index],'r')
    dset=f['Acoustic']        

    #extracting metadata and the DAS-data (hereby defined "data")
    attr=dset.attrs 
    data = np.array(dset)
    #choosing a taper function (important for the fourier transform)
    taper = signal.windows.tukey(seg_len, 0.25)
    
    if file_index!=nFiles-1:

        g = h5py.File(q[file_index+1],'r')
        dset2=g['Acoustic']
        data2= np.array(dset2)
        #Concatenating data with data from adjacent file
        #because one of the windows overlapping with
        #the next file

        data = np.concatenate((data, data2[0:(nu*file_length)//(nseg)]), axis=0)
        Fsegs = channel_fourier(data, nseg, seg_len, taper, ind_a, ind_e, ind_f)

    else:
        Fsegs = channel_fourier(data, nseg-1, seg_len, taper, ind_a, ind_e, ind_f)

    return Fsegs

 

if __name__=='__main__':
    print("Max number of CPUs: ", multiprocessing.cpu_count())
    

    # only relevant if you want a cube for a specific 
    # time interval
    day= 24
    month=7
    sec=0
    minut=0
    hours=0

    print(f"Processed day: {day}.{month}.2020")

    startT = time.time() 


    seg_length=1/d_f #calculate window length corresponding to d_f
    
    ind_f = int(seg_length*freq_max+1)

    seg_len=int(seg_length*nu) #how many time points should be in one processing window

    nseg=int(2*(file_length/seg_length)) #amount of segments for the desired window length

    location_coords = np.arange(loc_a, loc_e, 4)
    freq_coords=scipy.fft.rfftfreq(int(nu/d_f), 1/nu)[:ind_f]



    # in next line the corresponding channel indices get calculated
    # since the distance between each channel is 4m
    ind_a, ind_e= loc_a//4, loc_e//4 

    # z-shape is in following order: (number time intervals, location, frequency)

    two_digit_month = str(month) if month >9 else "0"+str(month)
    two_digit_day = str(day) if day >9 else "0"+str(day)
        
    path_of_zarr = ZARR_NAME
    root = zarr.open(path_of_zarr, mode="w")
    z = root.zeros('data', shape=(nFiles*nseg-1,ind_e-ind_a,ind_f), chunks=(nseg,ind_e-ind_a,ind_f), dtype='float64')

    
    # getFilenames returns a dictionary full of the filenames
    # which are in the folder /RhoneData/YYYYMMDD

    filenames = getFilenames(day,month,DATA_PATH)

    dummy_file = h5py.File(filenames[0], 'r')
    attr = dummy_file['Acoustic'].attrs
    zeit = datetime.fromisoformat(attr["ISO8601 Timestamp"])
    time_coords = [(zeit+timedelta(i*file_length/nseg)).replace(tzinfo=None) for i in range (nFiles*nseg-1)]
    #time_coords = time_coords

    metadata=root.create_group('metadata')
    metadata.create_dataset('loc_coords', data=location_coords)
    metadata.create_dataset('freq_coords', data=freq_coords)
    time_coordinates=metadata.empty('time_coords', shape=(len(time_coords)),dtype='M8[D]')
    for i,times in enumerate(time_coords):
        time_coordinates[i]=times


    # In the following lines, multiple cpu-cores calculate 
    # a fft for each file simultanously.
    
    # Before that we split the whole data to be processed in 
    # digestible parts (of course arbitrary) to not overload the memory!
    # Here, we divide a list of indices (1,2,3,...,nFiles) into 13 parts.  
    # Since create_spectro_segment needs an index k specifying which 
    # k-th file needs to be processed, we now pass on those index
    # lists to all the cores of the processors so each file can 
    # be processed independently. Make the amount of diveded parts
    # bigger (in this case 13) if the memory can't hold that much data.

    n_div = 13
    index_list = np.arange(nFiles)
    
    for liste in np.array_split(index_list, n_div):
        with Pool(nCores) as p:
            res = p.map(partial(create_spectro_segment, d_f=d_f, nseg=nseg, seg_len=seg_len, ind_a=ind_a, ind_e= ind_e, ind_f = ind_f, q= filenames)  , liste)
        res = list(res)
        for i in liste:
            if i!=nFiles-1:
                z[i*nseg:(i+1)*nseg]=res[i-int(liste[0])]
            else:
                z[i*nseg:(i+1)*nseg-1]=res[i-int(liste[0])]
 
    print("Time for processing all this:", (-startT + time.time()))
    print("Number of processed files:", nFiles)
    print("Number of used cores:", nCores)
    print("Time per File: ", ((-startT + time.time())/nFiles))

    





        

