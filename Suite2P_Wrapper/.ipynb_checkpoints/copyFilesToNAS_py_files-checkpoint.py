import shutil
import os


# tif data folder where tif files will be extracted from raw data
tif_data_folder = 'E:/Users/samsoon.inayat/S_Drive/Tif_Data'
# processed_data_folder where the results of suite 2p will be stored
#processed_data_folder = 'E:/Users/samsoon.inayat/OneDrive - University of Lethbridge/pySuite2P/processed_data'
processed_data_folder = 'E:/Users/samsoon.inayat/S_Drive/Processed_Data'
nas_processed_data_folder = '//mohajerani-nas.uleth.ca/storage/homes/brendan.mcallister/2P/Processed_Data'

#E:\Users\samsoon.inayat\S_Drive\Processed_Data
# Directory name where raw data is
#dir_name = '//mohajerani-nas.uleth.ca/storage/homes/samsoon.inayat/Data/183628/2019-07-01/1_003'

filename = 'recording_list_10.txt'
f = open(filename,'r')
dir_names = []
for drs in f:
    temp = drs[:-1]
    tempp = temp.replace(os.sep,os.altsep)
    dir_names.append(tempp)
f.close()
# print(dir_names)
ii = 1
for x in dir_names:
    dir_name = x
    str_pos= dir_name.find('Data')
    slash_pos = dir_name[str_pos:].find("/")
    rec_info = dir_name[str_pos+slash_pos:]
    spdir_name = processed_data_folder + rec_info
    dpdir_name = nas_processed_data_folder + rec_info


    ssuite2p_folder = spdir_name + '/suite2p'
    dsuite2p_folder = dpdir_name + '/suite2p'
    if not os.path.isdir(ssuite2p_folder):
        continue

    dir_list = os.listdir(ssuite2p_folder)
    # print(dir_list)
    ii = ii + 1
    for nn in dir_list:
        str_present = nn.find('plane')
        if str_present < 0:
            continue
        splane_dir = ssuite2p_folder + '/' + nn
        #print(splane_dir)
        dplane_dir = dsuite2p_folder + '/' + nn
        #print(dplane_dir)
        dir_list_plane = os.listdir(splane_dir)
        for pp in dir_list_plane:
            str_present = pp.find('.npy')
            if str_present < 0:
                continue
            sfilename = splane_dir + '/' + pp
            dfilename = dplane_dir + '/' + pp
            # print(sfilename)
            if not os.path.isdir(dplane_dir):
                os.makedirs(dplane_dir)
            if os.path.isfile(sfilename) and not os.path.isfile(dfilename):
                print(dfilename)
                shutil.copyfile(sfilename,dfilename)

    #te = Thor_Experiment.Thor_Exp(dir_name,processed_data_folder,tif_data_folder,nas_processed_data_folder)
