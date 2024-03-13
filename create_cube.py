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
This script is used to create a datacube containing
a spectrogram for each DAS-channel. However, since it is
programmed for the structure of our DAS-data, 
some things need to be changed to make it suitable for
other structures. It is recommended to run this with 
as many CPU-cores as possible (and very much memory).

How it's used in terminal:

create_cube.py 


If you have any questions, feel free to reach out to
me, my mail is felix.roth@studserv.uni-leipzig.de

"""

# Key adjustable parameters:

# d_f: Frequency resolution in Hz. Adjusts the frequency granularity of the spectrogram.
d_f = 1

# loc_a and loc_e: Specify the cable section to be processed (in meters).
# 0 signifies the start of the cable. Adjust these to focus on a specific segment.
loc_a, loc_e = 0, 9200

# DATA_PATH: Path to the directory where the h5 data files are stored.
# Adjust this to point to your data location.
DATA_PATH = "rhonedata/"

# nFiles: Determines how many h5 files are processed. For example, setting nFiles=10
# processes the first ten chronologically sorted h5 files in the DATA_PATH directory.
nFiles = 2

# nCores: Specifies how many CPU cores should be used for computation. Adjust according
# to your system's capabilities and the workload size. The maximum number can be set to
# multiprocessing.cpu_count() to use all available cores.
nCores = 8

# Additional parameters:

# file_length: Length of a single h5 file in seconds.
file_length = 30

# nu: Sampling frequency in Hz. Defines the rate at which data points are recorded.
nu = 1000

# freq_max: Maximum frequency to be considered in the analysis. All values above
# this frequency will be cut off.
freq_max = 100

 #path and name of resulting cube (which is a zarr-file)
ZARR_NAME = "cryo_cube.zarr"


def getTimeFromFilename(filename, day, month):
    """
    purpose of this function:
    
    takes in filename of a DAS-h5-file
    and returns the associated time of the data
    (each h5 file contains 30 s of data)
    
    structure of h5-files:
    rhone1khz_UTC_YYYYMMDD_HHMMSS.***.h5

    the function returns the time difference between the timepoint 
    of the file and the start of the day
    """
    hours_minutes_seconds = re.findall("_......\\.", filename)[0][1:]
    #print(hours_minutes_seconds)
    hour = int(hours_minutes_seconds[0:2])
    minutes= int(hours_minutes_seconds[2:4])
    seconds= int(hours_minutes_seconds[4:6])
    return (datetime(2020, month, day, hour, minutes, seconds) - datetime(2020, month, day, 0, 0, 0)).total_seconds()



def getFilenames(day, month, mypath):
    """
    This function is supposed to collect all
    filenames in the folder with the DAS-data
    and returns a dictionary where the keys are
    integers and the corresponding values
    are the filenames. The dictionary is then sorted
    chronologically.
    """
    
    two_digit_month = str(month) if month >9 else "0"+str(month)
    #mypath = "/work/users/fu591lmui/RhoneData/2020"+two_digit_month+str(day)+"/"

    files = list([f for f in listdir(mypath) if isfile(join(mypath, f))])
    q = list([(mypath+file, getTimeFromFilename(file, day, month)) for file in files])
    q = sorted(q, key = itemgetter(1))
    return  {i:qq[0] for i,qq in enumerate(q)}




def channel_fourier(data, nseg, seg_len, taper, ind_a, ind_e, ind_f):
    
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
    if not os.path.exists("data"):
        os.makedirs("data")
        
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
    # 13 parts (of course arbitrary) to not overload the memory!
    # How it's done in this case in this case is to divide
    # a list of indices (1,2,3,...,nFiles) into 13 parts. Since 
    # create_spectro_segment needs an index k which tells which 
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

    





        

