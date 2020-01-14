# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 15:29:10 2019

@author: samsoon.inayat
"""

from progressbar import ProgressBar


pbar = ProgressBar()

for ii in pbar(range(10)):
    print(ii)
    
for ii in pbar(range(100)):
    print(ii)