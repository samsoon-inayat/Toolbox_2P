#import matplotlib.pyplot as plt
#import numpy as np

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
nas_processed_data_folder = 'Y:/homes/samsoon.inayat/S_Drive/pySuite2p_Processed_Data'
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


for x in dir_names:
    print(x)
    dir_name = x
    te = Thor_Experiment.Thor_Exp(dir_name,processed_data_folder,tif_data_folder,nas_processed_data_folder)
    data_processing_status = 0
    for root, dirs, files in os.walk(te.nas_pd_dir):
        for name in files:
            if name.endswith(("Fall.mat")):
                data_processing_status = 1
                break
    if data_processing_status == 1:
        continue
#    if te.exp_params.get('zFastEnable') == '1':
#        continue
    ops = run_s2p.default_ops()
#        print(ops)
    ops.update({'nplanes':te.exp_params.get('nplanes')})
    ops.update({'save_path0':te.pd_dir})
    ops.update({'fs':float(te.exp_params.get('frameRate'))})
    ops.update({'save_mat':True})
    mat = scipy.io.loadmat(te.pd_dir + '/bidishift.mat')
    bidishift = mat['bidishift'];
#    print(bidishift[0][0])
#    break
    ops.update({'do_bidiphase':True})
    ops.update({'bidiphase':bidishift[0][0]})
    ops.update({'roidetect':True})
    ops.update({'do_registration':1})
    print(ops)
    db = {
            'h5py': [], # a single h5 file path
            'h5py_key': 'data',
            'data_path': [te.tif_dir_name], # a list of folders with tiffs
                                         # (or folder of folders with tiffs if look_one_level_down is True, or subfolders is not empty)
            }
    np.save('ops.npy', ops)
    np.save('db.npy', db)
    subprocess.call(['python','-u','-W','ignore','-m','suite2p','--ops', 'ops.npy','--db', 'db.npy'],shell=True)
#    break
    # opsEnd = run_s2p.run_s2p(ops = ops,db = db
    if os.path.exists(te.tif_dir_name):
        shutil.rmtree(te.tif_dir_name)
    else:
        print('tif folder already removed')
    del(te)
    



#dir_path = 'explorer "%s"'%processed_data_folder
#subprocess.Popen(dir_path)

#for ii in range(10,20):
#    print(ii)
