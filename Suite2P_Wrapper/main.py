#import matplotlib.pyplot as plt
#import numpy as np

import Thor_Experiment
#import os
import subprocess

# tif data folder where tif files will be extracted from raw data
tif_data_folder = 'T:/py_suite2p/tif_data_folder'
# processed_data_folder where the results of suite 2p will be stored
processed_data_folder = 'T:/py_suite2p/processed_data'

# Directory name where raw data is
dir_name = r'//mohajerani-nas.uleth.ca/storage/homes/samsoon.inayat/Data/183628/2019-07-01/1_003'

te = Thor_Experiment.Thor_Exp(dir_name,processed_data_folder,tif_data_folder)



#dir_path = 'explorer "%s"'%processed_data_folder
#subprocess.Popen(dir_path)

#for ii in range(10,20):
#    print(ii)