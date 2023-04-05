for ii = 1:60
    pause(1);
end
%% add path
add_to_path
%% clear all if needed 
clear all
% clc
%% Control data folder reading
data_folder1 = 'E:\Data';
data_folder2 = 'E:\Data';
processed_data_folder{1} = 'E:\Processed_Data_AD_paper\Control';
processed_data_folder{2} = 'E:\Processed_Data_AD_paper\Matlab\Control';
processed_data_folder{3} = 'E:\Processed_Data_AD_paper\Matlab_bw3\Control';
animal_list_control = {'173706';'183761';'183745';'183628';'183762'};
date_list_control = {'2019-02-14';'2019-06-19';'2019-06-19';'2019-06-21';'2019-06-24'};
[dS_C,T_C] = get_exp_info_from_folder(data_folder1,processed_data_folder,animal_list_control,date_list_control);
[dS_C1,T_C1] = get_exp_info_from_folder(data_folder2,processed_data_folder,animal_list_control,date_list_control);
T_C = [T_C;T_C1];
T_C = reduce_table(T_C,animal_list_control,date_list_control);
T_C(1:2:end,:) = [];
disp('Done');
%% APP data folder reading
animal_list_APP = {'183227';'183228';'183329';'001567';'001569'};
date_list_APP = {'2019-06-26';'2019-06-27';'2019-06-25';'2020-04-26';'2020-04-26'};
data_folder1 = 'E:\Data';
processed_data_folder{1} = 'E:\Processed_Data_AD_paper\APP';
processed_data_folder{2} = 'E:\Processed_Data_AD_paper\Matlab\APP';
processed_data_folder{3} = 'E:\Processed_Data_AD_paper\Matlab_bw3\APP';
[dS_A,T_A] = get_exp_info_from_folder(data_folder1,processed_data_folder,animal_list_APP,date_list_APP);
T_A = reduce_table(T_A,animal_list_APP,date_list_APP);
disp('Done');
%% meta data setup
colormaps = load('../../Common/Matlab/colorblind_colormap.mat');
colormaps.colorblind = flipud(colormaps.colorblind);
mData.colors = mat2cell(colormaps.colorblind,[ones(1,size(colormaps.colorblind,1))]);%{[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'}; % mData.colors = getColors(10,{'w','g'});
mData.dcolors = mat2cell(distinguishable_colors(20,'w'),[ones(1,20)]);
mData.axes_font_size = 6; mData.sigColor = [0.54 0.27 0.06];
mData.shades = generate_shades(3);
mData.conj_comp_colors = [mData.dcolors(9);mData.colors([3 5])];
% display_colors(mData.colors);
Uleth_one_drive = 'Z:\homes\brendan.mcallister\2P';
Uleth_one_drive = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P';
Uleth_one_drive = 'E:\PostProcessing\AD_paper';
mData.pdf_folder = [Uleth_one_drive '\PDFs']; 
mData.pd_folder = [Uleth_one_drive '\PD_Matlab'];
disp('Done');
mData.magfac = 1;
%%
%%
if 0
%     make_db(T_C);
    process_abf(T_C,0);
    process_abf(T_A,0);
end
disp('Done');

%% Load Primary Data
ei10_C = getData_py_2(T_C);
ei10_A = getData_py_2(T_A);
% Load secondary processed data (raster plots and MI etc.,)
dcfilename = 'define_contexts.m';
%%
nowr = [0 0 0];
owr = [1 1 1];
binwidths = [0.3 3];
for ii = 1:length(ei10_C)
    ei10_C(ii) = make_and_load_rasters(ei10_C(ii),binwidths,nowr);
end
for ii = 1:length(ei10_A)
    ei10_A(ii) = make_and_load_rasters(ei10_A(ii),binwidths,nowr);
end
disp('Done');
% send_email('samsoon.inayat@gmail.com','Done');
%% Load Training Data
% ei10_C = loadContextsResponses_ctrl_old(ei10_C,[1 1],[0 0 0]);
% ei10_A = loadContextsResponses_ctrl_old(ei10_A,[1 1],[0 0 0]);
training_data_C = behaviorProcessor;
training_data_A = behaviorProcessor_AD;

%% load calcium signals into each recording
pl = 1;
for ii = 1:length(ei10_C)
    [~,caSig_C{ii,pl}] = get_calcium_data_raw(ei10_C{ii},pl);
end
for ii = 1:length(ei10_A)
    [~,caSig_A{ii,pl}] = get_calcium_data_raw(ei10_A{ii},pl);
end
%% load calcium signals pl 2
[~,tempcs] = get_calcium_data_raw(ei10_C{1},2);
caSig_C{1} = [caSig_C{1};tempcs];
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

