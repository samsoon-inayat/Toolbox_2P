import Thor_Experiment
import abf_processor
import os
import numpy, scipy.io

m_pd_dir = '//mohajerani-nas.uleth.ca/storage2/homes/samsoon.inayat/S_Drive/pySuite2p_Processed_Data'
#d_dir = '//mohajerani-nas.uleth.ca/storage/homes/samsoon.inayat/Data/183745/2019-07-01/1'
pd_dir = '//mohajerani-nas.uleth.ca/storage2/homes/samsoon.inayat/S_Drive/pySuite2p_Processed_Data'



filename = 'recording_list_r.txt'
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
    te = Thor_Experiment.Thor_Exp(dir_name,pd_dir,'')
    if hasattr(te,'abf'):
        abf_data = te.abf.channel_data
        abf_data.update({'si':te.abf.si})
        abf_data.update({'number_of_samples':te.abf.number_of_samples})
        fileName = os.path.join(te.pd_dir.replace(os.sep,os.altsep),'abf_data.mat').replace(os.sep,os.altsep)
        scipy.io.savemat(fileName, mdict = abf_data)
    del(te)
