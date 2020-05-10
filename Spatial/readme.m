% flow chart for processing 2P data
%%
% first run "mainRunner_2P" to process raw data acquired from 2P microscope
% and obtain ROIs respresenting cells using Suite2P. You have to set the
% animal id, date, and recording number. 
% In addition, you have to set the required Suite2P files make_db.m 
% Also set metaD1.xml which I have setup to tell the program information about the abf files

mainRunner_2P

%%
% Run new_main of Suite2P to manually curate the data and identify rois
% which representative putative cells
new_main

%%
post_processing