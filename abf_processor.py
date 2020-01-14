# -*- coding: utf-8 -*-
import numpy as np
import pyabf
import pickle
from os import path

class abf_data:
    
    def __init__(self,abf_filename,pd_dir,stim_params):
        self.filename = abf_filename
        self.pd_filename = pd_dir + '/abf_data.pyp'
        if path.exists(self.pd_filename):
            print('Loading abf data '+self.pd_filename)
            self.channel_data = self.__load_abf_data()
            abf = pyabf.ABF(abf_filename)
            self.stim_params = stim_params
            self.number_of_samples = len(abf.data[0])
            self.si = 1e6/abf.dataRate
            self.dataRate = abf.dataRate
            self.abf = abf
        else:
            print('Loading abf file')
            abf = pyabf.ABF(abf_filename)
            print('Processing abf file')
            self.channel_data = self.__get_and_save_abf_data(abf,stim_params)
            self.stim_params = stim_params
            self.number_of_samples = len(abf.data[0])
            self.si = 1e6/abf.dataRate
            self.dataRate = abf.dataRate
           
        print('Successfuly loaded data \n')
        
    def load_pyabf(self): # if abf file object through pyabf needs to be generated later for any reason
        self.abf= pyabf.ABF(self.filename)
        
    def __load_abf_data(self):
        with open(self.pd_filename, 'rb') as f:
            channel_data = pickle.load(f)
        return channel_data
    
    def __find_rising_edge(self,signal,threshold,minimum_diff):
        dSignal = np.diff(signal)
        edge = (dSignal>=threshold).nonzero()
        edge = edge[0]
        temp = np.diff(edge)
        temp1 = (temp<=minimum_diff).nonzero()
        temp1 = temp1[0]+1
        return np.delete(edge,temp1)
         
    def __find_falling_edge(self,signal,threshold,minimum_diff):
        dSignal = np.diff(signal)
        edge = (dSignal<=threshold).nonzero()
        edge = edge[0]
        temp = np.diff(edge)
        temp1 = (temp<=minimum_diff).nonzero()
        temp1 = temp1[0]+1
        return np.delete(edge,temp1)
    
    def __get_rid_of_values(self,signal,difference):
        dSignal = np.diff(signal)
        inds = (dSignal>difference).nonzero()
        print(inds)
        inds = inds[0]+1
        return np.delete(signal,inds)
    
    def __get_rid_of_close_repetitions(self,edge,minimum_diff):
        dEdge = np.diff(edge)
        temp = (dEdge<=minimum_diff).nonzero()
        temp = temp[0]+1
        np.delete(edge,temp)
        return edge,temp
    
    def __get_and_save_abf_data(self,abf,stim_params):
        o = {}
        if abf.channelCount != len(stim_params):
            raise Exception('Channel count in abf is not equal to number of channels in stim param xml file')
        ecd = {}
        for ii in range(abf.channelCount):
            abf.adcNames[ii]= stim_params[ii]
            f_edges = self.__find_falling_edge(abf.data[ii],-0.5,2)
            r_edges = self.__find_rising_edge(abf.data[ii],0.5,2)
            
            if 'air_puff' in stim_params[ii] or 'stim' in stim_params[ii] or 'photo_sensor' in stim_params[ii]:
                cmdTxt = "o.update({'%s_raw':abf.data[ii]})"%stim_params[ii]
                eval(cmdTxt)
            
            if 'ch_' in stim_params[ii]:
                cmdTxt = "ecd.update({'%s':abf.data[ii]})"%stim_params[ii]
                eval(cmdTxt)
            
            if 'air_puff' in stim_params[ii]:
                [f_edges,inds] = self.__get_rid_of_close_repetitions(f_edges,5000)
                if len(inds)>0:
                    np.delete(r_edges,inds)
                [r_edges,inds] = self.__get_rid_of_close_repetitions(r_edges,5000)
                if len(inds)>0:
                    np.delete(f_edges,inds)
                if f_edges[0] < r_edges[0]: # for air puff the rising has to be first
                    np.delete(f_edges,0)
            
            if len(f_edges) != len(r_edges):
                diff = len(f_edges) - len(r_edges)
                if diff > 0:
                    f_edges = f_edges[:len(r_edges)]
                else:
                    r_edges = r_edges[:len(f_edges)]
    #            raise Exception('rising_falling different')
            cmdTxt = "o.update({'%s_f':f_edges})"%stim_params[ii]
            eval(cmdTxt)
            cmdTxt = "o.update({'%s_r':r_edges})"%stim_params[ii]
            eval(cmdTxt)
                
        dist = self.__process_encoder_signals(ecd.get('ch_a'),ecd.get('ch_b'))
    #    print(dist)
        o.update({'dist':dist})
        with open(self.pd_filename, 'wb+') as f:
            pickle.dump(o, f)
        return o
            
    def __process_encoder_signals(self,cha,chb):
        chat = cha > 2.5
        chat = 1*chat
#        print(chat)
        chbt = chb > 2.5
        chbt = 1*chbt
        encoder_count = 0
        valP = [chat[0],chbt[0]]
#        print(valP)
        dist = [];
        dist.append(encoder_count)
        for ii in range(1,len(cha)):
            valC = [chat[ii],chbt[ii]]
            if valC == valP:
                valP = valC
                dist.append(encoder_count)
                continue
            if valP == [0,0] and valC == [1,0]:
                encoder_count = encoder_count+1
                valP = valC
                dist.append(encoder_count)
                continue
            if valP == [0,0] and valC == [0,1]:
                encoder_count = encoder_count-1
                valP = valC
                dist.append(encoder_count)
                continue
            if valP == [1,0] and valC == [1,1]:
                encoder_count = encoder_count+1
                valP = valC
                dist.append(encoder_count)
                continue
            if valP == [1,0] and valC == [0,0]:
                encoder_count = encoder_count-1
                valP = valC
                dist.append(encoder_count)
                continue
            if valP == [0,1] and valC == [0,0]:
                encoder_count = encoder_count+1
                valP = valC
                dist.append(encoder_count)
                continue
            if valP == [0,1] and valC == [1,1]:
                encoder_count = encoder_count-1
                valP = valC
                dist.append(encoder_count)
                continue
            if valP == [1,1] and valC == [0,1]:
                encoder_count = encoder_count+1
                valP = valC
                dist.append(encoder_count)
                continue
            if valP == [1,1] and valC == [1,0]:
                encoder_count = encoder_count-1
                valP = valC
                dist.append(encoder_count)
                continue
            dist.append(encoder_count)
#            print(ii)
        return np.array(dist)
