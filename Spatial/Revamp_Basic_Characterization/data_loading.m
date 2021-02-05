for ii = 1:60
    pause(1);
end
clear all
% add_to_path
% clc
%%
[f,cName] = getFolders;
colormaps = load('../MatlabCode/colorblind_colormap.mat');
colormaps.colorblind = flipud(colormaps.colorblind);
% mData.colors = mat2cell(colormaps.colorblind,[ones(1,size(colormaps.colorblind,1))]);%
mData.colors = {[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'}; % mData.colors = getColors(10,{'w','g'});
mData.axes_font_size = 6; mData.sigColor = [0.54 0.27 0.06]; 
Uleth_one_drive = 'G:\OneDrives\OneDrive - University of Lethbridge';
mData.pdf_folder = [Uleth_one_drive '\PDFs']; 
mData.pd_folder = [Uleth_one_drive '\ProcessedData\Matlab'];
disp('Done');
%% protocol 10
temp = load('T_10_All.mat');
T10 = temp.T;
selRecs10 = [8 9 10 11];
ET10 = T10(selRecs10,:);
data10 = getData_py(f,ET10(4,:));
binWidths = [0.2,1.5];
data10_1 = load_contexts_responses(data10,'define_contexts_1.m',binWidths);
disp('done')
%% Protocol 15
temp = load('T_15_All.mat');
T15 = temp.T;
d15 = getData_py(f,T15(8,:));
binWidths = [0.2,1.5];
d15_1 = load_contexts_responses(d15,'define_contexts.m',binWidths);
