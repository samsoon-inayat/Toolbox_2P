#import matplotlib.pyplot as plt
#import numpy as np

import Thor_Experiment
#import os
import subprocess

# tif data folder where tif files will be extracted from raw data
tif_data_folder = 'E:/Users/samsoon.inayat/S_Drive/Tif_Data'
# processed_data_folder where the results of suite 2p will be stored
processed_data_folder = 'E:/Users/samsoon.inayat/OneDrive - University of Lethbridge/pySuite2P/processed_data'

# Directory name where raw data is
dir_name = '//mohajerani-nas.uleth.ca/storage/homes/samsoon.inayat/Data/183628/2019-07-01/1_003'

filename = 'recording_list.txt'
f = open(filename,'r')
dir_names = []
for drs in f:
    dir_names.append(drs)
f.close()

for x in dir_names:
    print(x)
    dir_name = x[:-1]
    te = Thor_Experiment.Thor_Exp(dir_name,processed_data_folder,tif_data_folder)
    del(te)



#dir_path = 'explorer "%s"'%processed_data_folder
#subprocess.Popen(dir_path)

#for ii in range(10,20):
#    print(ii)