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
processed_data_folder{1} = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_AD_paper\Control';
processed_data_folder{2} = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_AD_paper\Matlab_bw3\Control';
animal_list_control = {'173706';'183761';'183745';'183628';'183762'};
date_list_control = {'2019-02-14';'2019-06-19';'2019-06-19';'2019-06-21';'2019-06-24'};
[dS_C,T_C] = get_exp_info_from_folder(data_folder1,processed_data_folder,animal_list_control,date_list_control);
[dS_C1,T_C1] = get_exp_info_from_folder(data_folder2,processed_data_folder,animal_list_control,date_list_control);
T_C = [T_C;T_C1];
T_C = reduce_table(T_C,animal_list_control,date_list_control);
disp('Done');
%%
animal_list_APP = {'183227';'183228';'183329';'001567';'001569'};
date_list_APP = {'2019-06-26';'2019-06-27';'2019-06-25';'2020-04-26';'2020-04-26'};
data_folder1 = '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data';
processed_data_folder{1} = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_AD_paper\APP';
processed_data_folder{2} = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_AD_paper\Matlab_bw3\APP';
[dS_A,T_A] = get_exp_info_from_folder(data_folder1,processed_data_folder,animal_list_APP,date_list_APP);
T_A = reduce_table(T_A,animal_list_APP,date_list_APP);
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
%%
if 0
%     make_db(T_C);
    process_abf(T_C,0);
    process_abf(T_A,0);
end
disp('Done');

%%
ei10_C = getData_py_2(T_C);
ei10_A = getData_py_2(T_A);
%%
nowr = [0 0 0];
owr = [1 1 1];
binwidths = [0.2 3];
for ii = 1:length(ei10_C)
    ei10_C(ii) = make_and_load_rasters(ei10_C(ii),binwidths,nowr);
end
for ii = 1:length(ei10_A)
    ei10_A(ii) = make_and_load_rasters(ei10_A(ii),binwidths,nowr);
end

%%
% ei10_C = loadContextsResponses_ctrl_old(ei10_C,[1 1],[0 0 0]);
% ei10_A = loadContextsResponses_ctrl_old(ei10_A,[1 1],[0 0 0]);
training_data_C = behaviorProcessor;
training_data_A = behaviorProcessor_AD;






%% OLD Code below
% this is just to load behavior

for ii = 1:size(T10.T,1)
    if ismember(ii,[1:size(T10.T,1)]) % select which data to load in the second argument
        eiB(ii) = getBehavior(f,T10.T(ii,:));
    end
end
disp('Done!');

%%
% This is just to visualize behavior graphs
inds = [];
for ii = 1:size(T10.T,1)
    if ismember(ii,[1:size(T10.T,1)]) % select which data to load in the second argument
        if isempty(eiB{ii})
            continue;
        end
         behaviorPlot(eiB(ii))
         ii
         key = getkey;
         if key == 27 % esc
             break;
         end
         if key == 105 %i
             inds = [inds ii];
         end
    end
end
inds
disp('Done!');
