# -*- coding: utf-8 -*-
"""
Created on Sun Dec 21 23:39:08 2025

@author: samso
"""

import scipy.io as sio
import pickle


with open("session_data.pkl", "rb") as f:
    loaded_data = pickle.load(f)

sio.savemat("session_data.mat", {"animal": loaded_data["animal"]})
