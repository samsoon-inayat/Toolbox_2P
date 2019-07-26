# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 09:13:11 2019

@author: samsoon.inayat
"""


import os
from os import path
import xml.etree.ElementTree as ET
import abf_processor
import numpy as np
#import matplotlib.pyplot as plt
#import time
from PIL import Image
from progressbar import ProgressBar
import pickle


#mat = scipy.io.loadmat('file.mat')


class Thor_Exp:

    def __init__(self,dir_name,p_dir_name_main,tif_data_folder):
#        print(tif_data_folder)
        self.d_dir = dir_name # data directory name
        self.pd_dir_main = p_dir_name_main #
        self.recording_info = self.__get_recording_info()
        self.pd_dir = p_dir_name_main + self.recording_info
        if not os.path.exists(self.pd_dir):
            os.makedirs(self.pd_dir)
        self.processing_status_filename = self.pd_dir + '/processing_status.pyp'
        self.processing_status = self.__load_processing_status()
        self.exp_params = self.__get_experiment_parameters()
        if self.exp_params.get('zFastEnable') == '1':
            return
        self.raw_filename = self.__get_raw_file_name()
        self.stim_params = self.__get_stim_parameters()
        self.__save_processing_status('abf',0)
        self.abf = abf_processor.abf_data(self.__get_abf_file_name(),self.pd_dir,self.stim_params)
        self.__save_processing_status('abf',1)

        self.__save_processing_status('tif',0)
        if not tif_data_folder.strip():
            self.tif_dir_name = ''
            self.processing_status.update({'tif':-1})
            if self.exp_params.get('zFastEnable') == '1': # check to see if multiple plane data exists
                self.exp_params.update({'nplanes':10})
                raise Exception('You need to write code for this')
            else:
                self.exp_params.update({'nplanes':1})
                self.exp_params['timepoints'] = len(self.abf.channel_data.get('frames_f'))-1
                self.exp_params['frames'] = len(self.abf.channel_data.get('frames_f'))-1
		# temp = self.exp_params.get('nplanes')
		# print(temp)
        else:
            self.tif_dir_name = tif_data_folder + self.recording_info
            self.__raw_to_tif()
            self.__save_processing_status('tif',1)
        # self.__run_suite2p()


#     def __run_suite2p(self):
#         ops = run_s2p.default_ops()
# #        print(ops)
#         ops.update({'nplanes':self.exp_params.get('nplanes')})
#         ops.update({'save_path0':self.pd_dir})
#         ops.update({'fs':float(self.exp_params.get('frameRate'))})
#         ops.update({'save_mat':True})
#         mat = scipy.io.loadmat(self.pd_dir + '/bidishift.mat')
#         bidishift = mat['bidishift'];
#         ops.update({'do_bidiphase':True})
#         ops.update({'bidiphase':bidishift[0][0]})
#         ops.update({'roidetect':True})
#         ops.update({'do_registration':1})
#         print(ops)
#         self.__save_processing_status('s2p',0)
# #        print(self.tif_dir_name)
#         db = {
#                 'h5py': [], # a single h5 file path
#                 'h5py_key': 'data',
#                 'data_path': [self.tif_dir_name], # a list of folders with tiffs
#                                              # (or folder of folders with tiffs if look_one_level_down is True, or subfolders is not empty)
#                 }
#         opsEnd = run_s2p.run_s2p(ops = ops,db = db)
#         self.__save_processing_status('s2p',1)


    def __save_processing_status (self,param_name,value):
        self.processing_status.update({param_name:value})
        with open(self.processing_status_filename, 'wb+') as f:
            pickle.dump(self.processing_status,f)

    def __load_processing_status (self):
        if path.exists(self.processing_status_filename):
            with open(self.processing_status_filename, 'rb') as f:
                return pickle.load(f)
        else:
            return {}


    def __raw_to_tif(self):

        databin_filename  = self.pd_dir + '/suite2p/plane0/data.bin'
        print(self.tif_dir_name)
        if not os.path.exists(self.tif_dir_name): # if tif directory doesn't exist create it
            os.makedirs(self.tif_dir_name)
        if self.exp_params.get('zFastEnable') == '1': # check to see if multiple plane data exists
            self.exp_params.update({'nplanes':10})
            raise Exception('You need to write code for this')
        else:
            self.exp_params.update({'nplanes':1})
            self.exp_params['timepoints'] = len(self.abf.channel_data.get('frames_f'))-1
            self.exp_params['frames'] = len(self.abf.channel_data.get('frames_f'))-1

            if path.exists(databin_filename):
                print('data.bin present skipping converting to tifs')
                return

            file_list = os.listdir(r"{}".format(self.tif_dir_name))
            if len(file_list) == self.exp_params.get('timepoints'):
                print('\n Raw to tif conversion already complete \n')
                return
            start_file_number = len(file_list)-10
            print(start_file_number)
            if start_file_number < 0:
                start_file_number = 0
            nFrames = int(self.exp_params.get('timepoints'))
            print(nFrames)
            rows = int(self.exp_params.get('pixelY'))
            cols = int(self.exp_params.get('pixelX'))
            print('\n Starting conversion raw to tif \n')
            f = open(self.raw_filename, 'rb')
            pbar = ProgressBar()
            for ii in pbar(range(nFrames)):
                if ii < start_file_number:
                    continue
                f.seek(ii*rows*cols*2)
                dbytes = np.fromfile(f,'uint16',rows*cols)
                frame = np.reshape(dbytes,(rows,cols))
                tif_filename = '/time%d_plane0_channel0.tif'%ii
                tif_filename = self.tif_dir_name + tif_filename
                im = Image.fromarray(frame)
                im.save(tif_filename)
    #            print(frame.shape)
    #            plt.figure(10)
    #            plt.imshow(frame)
    #            plt.title(ii)
    #            plt.show()
    #            time.sleep(0.300)
            f.close()
            print('\n Conversion of raw to tif complete \n')

    def __get_experiment_parameters(self):
        print('Fetching the following parameters from Experiment.xml file \n')
        params = [['LSM','frameRate','pixelX','pixelY','widthUM','heightUM'],
                  ['Timelapse','timepoints'],
                  ['ExperimentNotes','text'],
                  ['Streaming','zFastEnable','frames'],
                  ['ZStage','steps']]
        root = self.__get_xml_root()
        param_listd = {}
        len_params = len(params)
        for ii in range(len_params):
            temp = params[ii]
            print(temp)
            sectionName = temp[0]
            these_params = temp[1:]
            for param_name in these_params:
                fr = self.__get_parameter(root,sectionName,param_name)
                param_listd.update({param_name:fr})
        print('\n Done fetching parameters \n')
        return param_listd


    def __get_stim_parameters(self):
        filename = self.__get_stim_file_name()
        tree = ET.parse(filename)
        root = tree.getroot()
        param_list = []
        for child in root.findall('channels'):
            for child1 in child:
                value = child1.get('name')
                param_list.append(value)
        return param_list


    def __get_parameter(self,root,child_name,param_name):
        for child in root.findall(child_name):
            value = child.get(param_name)
        return value


    def __get_xml_root(self):
        dir_name = self.d_dir
        xml_filename = 'Experiment.xml'
        file_with_path = os.path.join(dir_name, xml_filename)
        print(file_with_path)
        tree = ET.parse(file_with_path)
        root = tree.getroot()
        return root

    def __get_raw_file_name(self):
        dir_name = self.d_dir
    #    'Image_0001_0001.raw'
        dlist = os.listdir(dir_name)
        matching = [s for s in dlist if 'raw' in s]
        matching = [s for s in matching if 'Image' in s]
        return dir_name + '/' + matching[0]


    def __get_abf_file_name(self):
        dir_name = self.d_dir
    #    'Image_0001_0001.raw'
        dlist = os.listdir(dir_name)
        matching = [s for s in dlist if 'abf' in s]
        return dir_name + '/' + matching[0]


    def __get_recording_info (self):
        dir_name = self.d_dir
        str_pos= dir_name.find('Data')
        slash_pos = dir_name[str_pos:].find("/")
        return dir_name[str_pos+slash_pos:]


    def __get_stim_file_name(self):
        dir_name = self.d_dir
        dlist = os.listdir(dir_name)
        matching = [s for s in dlist if 'xml' in s]
        matching = [s for s in matching if 'meta' in s]
        return dir_name + '/' + matching[0]
