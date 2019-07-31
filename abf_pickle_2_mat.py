import Thor_Experiment
import abf_processor
import os


m_pd_dir = '//mohajerani-nas.uleth.ca/storage2/homes/samsoon.inayat/S_Drive/pySuite2p_Processed_Data'
d_dir = '//mohajerani-nas.uleth.ca/storage/homes/samsoon.inayat/Data/183745/2019-07-01/1'
pd_dir = '//mohajerani-nas.uleth.ca/storage2/homes/samsoon.inayat/S_Drive/pySuite2p_Processed_Data'

te = Thor_Experiment.Thor_Exp(d_dir,pd_dir,'')

abf_data = te.abf.channel_data
abf_data.update({'si':te.abf.si})
abf_data.update({'number_of_samples':te.abf.number_of_samples})


import numpy, scipy.io
scipy.io.savemat(os.path.join(te.pd_dir,'abf_data.mat'), mdict = abf_data)
