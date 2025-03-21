% function post_processing
reload = 0;
%% load data table
if reload
%     add_to_path
    clear all
    clc
    [f,cName] = getFolders;
    load('T16.mat');
    animal_ids = cell2mat(T{:,1});
    ow = 0;
%     ei = getData(f,T([3 7 11 13 15 16 18 20 21],:));
    ei = getData_py(f,T(1,:));
%     eip2 = getData_py(f,T([17 19],:));
    disp('Done');
end

binWidths = [0.1,1.5];
ei1 = load_contexts_responses(ei,'define_contexts_1.m',binWidths);

ei1 = find_raster_properties(ei1);

%%
colormaps = load('../MatlabCode/colorblind_colormap.mat');
colormaps.colorblind = flipud(colormaps.colorblind);
% mData.colors = mat2cell(colormaps.colorblind,[ones(1,size(colormaps.colorblind,1))]);%
mData.colors = {[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'}; % mData.colors = getColors(10,{'w','g'});
mData.axes_font_size = 6; mData.sigColor = [0.54 0.27 0.06]; 
Uleth_one_drive = 'G:\OneDrives\OneDrive - University of Lethbridge';
mData.pdf_folder = [Uleth_one_drive '\PDFs']; 
mData.pd_folder = [Uleth_one_drive '\ProcessedData\Matlab'];
disp('Done');