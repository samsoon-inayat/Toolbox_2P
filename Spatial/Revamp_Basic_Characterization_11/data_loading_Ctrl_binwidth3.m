for ii = 1:60
    pause(1);
end
%%
add_to_path
%%
clear all
% clc
%%

data_folder1 = '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data';
data_folder2 = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\Data';
processed_data_folder{1} = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_15';
processed_data_folder{2} = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_15\Matlab_bw3';
animal_list_control = {'183633';'183761';'183745';'183628';'183762'};
date_list_control = {'2019-06-04';'2019-06-06';'2019-06-07';'2019-06-11';'2019-06-11'};
[dS_C,T_C] = get_exp_info_from_folder(data_folder1,processed_data_folder,animal_list_control,date_list_control);
T_C = [T_C];
T_C = reduce_table(T_C,animal_list_control,date_list_control);
disp('Done');
%%
colormaps = load('../MatlabCode/colorblind_colormap.mat');
colormaps.colorblind = flipud(colormaps.colorblind);
mData.colors = mat2cell(colormaps.colorblind,[ones(1,size(colormaps.colorblind,1))]);%{[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'}; % mData.colors = getColors(10,{'w','g'});
mData.axes_font_size = 6; mData.sigColor = [0.54 0.27 0.06];

Uleth_one_drive = 'Z:\homes\brendan.mcallister\2P';
mData.pdf_folder = [Uleth_one_drive '\PDFs']; 
mData.pd_folder = [Uleth_one_drive '\ProcessedDataMatlab'];
disp('Done');

%%
if 0
%     make_db(T_C);
    process_abf(T_C,0);
end
disp('Done');

%%
ei = getData_py_2(T_C);
%%
binwidths = [0.2 3];
for ii = 1:length(ei)
    ei(ii) = make_and_load_rasters(ei(ii),binwidths,[0 0 0]);
end
