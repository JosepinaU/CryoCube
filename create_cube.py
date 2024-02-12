import h5py
import numpy as np
from scipy import signal
import scipy
from datetime import datetime, time, timedelta
import zarr
import time
import os
from os import listdir
from os.path import isfile, join
import re
from operator import itemgetter

import sys
from multiprocessing import Pool
import multiprocessing
from functools import partial

"""
This script is used to create a datacube containing
a spectrogram for each DAS-channel. However, since it is
programmed for the structure of our DAS-data, 
some things need to be changed to make it suitable for
other structures. It is recommended to run this with 
as many CPU-cores as possible. (And very much memory)

How it's used in terminal:

create_cube.py [day of data] 

where as [day of data] = 19 would be the data of the 19th of August
If you have any questions, feel free to reach out to
me, my mail is felix.roth@studserv.uni-leipzig.de

"""


def getTimeFromFilename(filename, day, month):
    """
    purpose of this function:
    
    takes in filename of a DAS-hd5-file
    and returns the associated time of the data
    (each hd5 file contains 30s of data)
    
    structure of hd5-files:
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



def getFilenames(day, month):
    """
    This function is supposed to collect all
    filenames in the folder with the DAS-data
    and returns a dictionary where the keys are
    integers and the corresponding values
    are the filenames. The dictionary is then sorted
    chronologically.
    """
    
    two_digit_month = str(month) if month >9 else "0"+str(month)
    mypath = "/work/users/fu591lmui/RhoneData/2020"+two_digit_month+str(day)+"/"

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
    
    Fsegs=np.zeros((nseg, ind_e- ind_a, ind_f))
    
    for j in range(nseg):
        for channel_number, channel  in enumerate(segs[j]):
            
            # note that modified_log(x)=10*log(x) (conversion to 
            
            fourier_transformed = np.fft.rfft(taper*channel, n=seg_len)
            fourier_transformed = (modified_log(np.abs(fourier_transformed))**2)[0:ind_f]
            fourier_transformed[0]=0
            Fsegs[j][channel_number]=fourier_transformed
    
    return Fsegs
   

def modified_log(x):
    return 10*np.log(x)

def create_spectro_segment(file_index, d_f, nseg, seg_len, ind_a, ind_e, ind_f, q):

    #opening the DAS-data
    f = h5py.File(q[file_index],'r')
    dset=f['Acoustic']        
    
    
    #extracting metadata and the DAS-data (hereby defined "data")
    attr=dset.attrs 
    data = np.array(dset)
    
    
    
    if file_index!=2880:
        g = h5py.File(q[file_index+1],'r')
        dset2=g['Acoustic']
        data2= np.array(dset2)
    

    #Concatenating data with data from adjacent file
    #because one of the windows is in the next file

    data = np.concatenate((data, data2[0:30000//(nseg)]), axis=0)


    #choosing a taper function (important for the fourier transform)
    taper = signal.windows.tukey(seg_len, 0.25)
    Fsegs = channel_fourier(data, nseg, seg_len, taper, ind_a, ind_e, ind_f)
    

    return Fsegs

 

if __name__=='__main__':
    print("max number of CPUs: ", multiprocessing.cpu_count())
    

    # get the day of interest out of command line
    day= int(sys.argv[1])
    month=7
    sec=0
    minut=0
    hours=0

    print("Processed day: ", day, ".", month, ".2020")

    startT = time.time() 
    # nFiles determines how many DAS-files get opened. One day holds 2880 files in 
    # our case. It is important to note that nFiles should always be smaller than 
    # the amount of files in the folder

    nFiles = 2879
    nCores = 100 # how many cores should be used for computation

    file_length=30 # length of file in seconds
    nu = 1000 # sampling frequence
    loc_a,loc_e = 0, 9200 #what cable length will be processed (in meter)
    freq_max = 100 #maximum frequency (all values above that get cut off)
    d_f = 1 # d_f is the desired frequency resolution


    number_of_parts=13 # in how many parts the data of the day 
    # should be divided into. If you have low memory, you can crank it up as much
    # as you like.
    
    
    
    # in next line the corresponding channel indices get calculated
    # since the distance between each channel is 4m
    ind_a, ind_e= loc_a//4, loc_e//4 


    
    

    seg_length=int(1/d_f) #calculate window length corresponding to d_f
    
    ind_f = seg_length*freq_max+1

    seg_len=seg_length*nu #how many time points should be in one processing window

    nseg=2*(file_length//seg_length) #amount of segments for the desired window length
    location_coords = np.arange(loc_a, loc_e, 4)
    freq_coords=scipy.fft.rfftfreq(int(nu/d_f), 1/nu)[:ind_f]
    time_coords= [] 


    # die Reihenfolge des z-shape ist (Anzahl Zeitsegmente, Location, Frequenz) 

    two_digit_month = str(month) if month >9 else "0"+str(month)
    two_digit_day = str(day) if day >9 else "0"+str(day)
    if not os.path.exists("data"):
        os.makedirs("data")
        
    path_of_zarr = '/work/users/fu591lmui/rhone/try1_create_spectro_local/19/data/'+two_digit_day+two_digit_month+'_1Hz_.zarr'
    z = zarr.open(path_of_zarr, mode='w', shape=(nFiles*nseg,ind_e-ind_a,ind_f), chunks=(nseg,ind_e-ind_a,ind_f), dtype='float64')
    print("segment length: ",seg_len)
    print("number of segments: ",nseg)

    
    # getFilenames returns a dictionary full of the filenames
    # which are in the folder /RhoneData/YYYYMMDD

    filenames = getFilenames(day,month)

    # In the following lines, multiple cpu-cores calculate 
    # a fft for each file simultanously.
    
    # Before that we split the whole data to be processed in 
    # 13 parts (of course arbitrary) to not overload the memory!
    # (here denoted as number_of_parts)
    # How it's done in this case in this case is to divide
    # a list of indices (1,2,3,...,nFiles) into 13 parts. Since 
    # create_spectro_segment needs an index k which tells which 
    # k-th file needs to be processed, we now pass on those index
    # lists to all the cores of the processors so each file can 
    # be processed independently. Make the amount of diveded parts
    # bigger (in this case 13) if the memory can't hold that much data.

    for liste in np.array_split(range(nFiles), number_of_parts):
        with Pool(nCores) as p:
            res = p.map(partial(create_spectro_segment, d_f=d_f, nseg=nseg, seg_len=seg_len, ind_a=ind_a, ind_e= ind_e, ind_f = ind_f, q= filenames)  , liste)
        res = list(res)
        for i in liste:
            z[i*nseg:(i+1)*nseg]=res[i-liste[0]]
 
    print("Time for processing all this:", (-startT + time.time()))
    print("number of processed files: ", nFiles)
    print("number of used cores: ", nCores)
    print("Time per File: ", ((-startT + time.time())/nFiles))

    





        

