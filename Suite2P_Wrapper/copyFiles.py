import Thor_Experiment
#import os
#import subprocess
from suite2p import run_s2p
import scipy.io
import shutil
import os
import numpy as np
import subprocess


# tif data folder where tif files will be extracted from raw data
tif_data_folder = 'E:/Users/samsoon.inayat/S_Drive/Tif_Data'
# processed_data_folder where the results of suite 2p will be stored
#processed_data_folder = 'E:/Users/samsoon.inayat/OneDrive - University of Lethbridge/pySuite2P/processed_data'
processed_data_folder = 'E:/Users/samsoon.inayat/S_Drive/Processed_Data'
nas_processed_data_folder = 'Z:/homes/brendan.mcallister/2P/Processed_Data'

#E:\Users\samsoon.inayat\S_Drive\Processed_Data
# Directory name where raw data is
#dir_name = '//mohajerani-nas.uleth.ca/storage/homes/samsoon.inayat/Data/183628/2019-07-01/1_003'

filename = 'recording_list.txt'
f = open(filename,'r')
dir_names = []
for drs in f:
    temp = drs[:-1]
    tempp = temp.replace(os.sep,os.altsep)
    dir_names.append(tempp)
f.close()
