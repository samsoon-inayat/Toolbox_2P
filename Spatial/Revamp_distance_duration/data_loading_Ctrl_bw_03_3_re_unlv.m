for ii = 1:60
    pause(1);
end
%%
add_to_path_unlv
disp('Done')
%%
clear all
% clc
%%
mainDrive = 'E:\GoogleDrive\InayatSamsoon\ULETH\';
% data_folder1 = '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data';
% data_folder2 = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\Data';
% processed_data_folder{1} = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_15';
% processed_data_folder{2} = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_15\MatlabRe';
% processed_data_folder{3} = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_15\MatlabReD';

data_folder1 = fullfile(mainDrive,'Data');%'\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data';
data_folder2 = 'E:\GoogleDrive\InayatSamsoon\ULETH\Data';%'\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\Data';
processed_data_folder{1} = 'E:\GoogleDrive\InayatSamsoon\ULETH\PData\Processed_Data_15';%'\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_15';
processed_data_folder{2} = 'E:\GoogleDrive\InayatSamsoon\ULETH\PData\Processed_Data_15\MatlabRe';%'\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_15\Matlab';
processed_data_folder{3} = 'E:\GoogleDrive\InayatSamsoon\ULETH\PData\Processed_Data_15\MatlabReD';%'\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_15\Matlab_bw3';

animal_list_control = {'183633';'183761';'183745';'183628';'183762'};
date_list_control = {'2019-06-04';'2019-06-06';'2019-06-07';'2019-06-11';'2019-06-11'};
[dS_C,T_C] = get_exp_info_from_folder(data_folder1,processed_data_folder,animal_list_control,date_list_control);
T_C = [T_C];
T_C1 = reduce_table(T_C,animal_list_control,date_list_control);
disp('Done');
%%
colormaps = load('../../Common/Matlab/colorblind_colormap.mat');
colormaps.colorblind = flipud(colormaps.colorblind);
mData.colors = mat2cell(colormaps.colorblind,[ones(1,size(colormaps.colorblind,1))]);%{[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'}; % mData.colors = getColors(10,{'w','g'});
mData.dcolors = mat2cell(distinguishable_colors(20,'w'),[ones(1,20)]);
mData.axes_font_size = 6; mData.sigColor = [0.54 0.27 0.06];
mData.asterisk_font_size = 9;
mData.shades = generate_shades(3);
mData.conj_comp_colors = [mData.dcolors(9);mData.colors([3 5])];
% display_colors(mData.colors);
% Uleth_one_drive = 'Z:\homes\brendan.mcallister\2P';
Uleth_one_drive = fullfile(mainDrive,'PostProcessing\Revamp_Dist_Dur_Paper');
% Uleth_one_drive = 'D:\OneDrive - University of Lethbridge\PDFs';
mData.pdf_folder = [Uleth_one_drive '\PDFs']; 
mData.pd_folder = [Uleth_one_drive '\PD'];
mData.magfac = 1;
disp('Done');
%%
sel_rec = [1 3 5];
sel_rec = 1:5;
ei = getData_py_2_top_plane(T_C1(sel_rec,:));
disp('Done');
%% delete second plane data
for ii = 1:5
    disp(length(ei{ii}.plane))
%     if length(ei{ii}.plane) > 1
%         ei{ii}.plane(2) = [];
%     end
end

%%
if 0
    ii = 1;
    edit_define_contexts_file(ei{ii});
end
dcfilename = 'define_contexts.m';
%%
if 0
%     make_db(T_C);
    process_abf(T_C,0);
end
disp('Done');
%%
binwidths = [0.3 3];
for ii = 1:length(ei)
    ei(ii) = load_context_info(ei(ii),binwidths,[0 0 0],dcfilename);
    ei(ii) = load_motion_correction_empty(ei(ii),binwidths,[0 0 0]);
end
disp('Done');

%%
for ii = 1:length(ei)
    ei(ii) = cellular_spatial_distances(ei(ii));
end
close all;
disp('Done');

%% new analysis
[udata,udata1] = get_unlv_analysis_data(ei);
% udataT = get_unlv_analysis_data_bin(udata1,0.3,'time');
% udataD = get_unlv_analysis_data_bin(udata1,3,'distance');
% udataT = get_unlv_analysis_data_time_bin(udata1,0.3);
% udataD = get_unlv_analysis_data_dist_bin(udata1,3);
% udataD = get_unlv_analysis_data_dist_bin(udata1,0.5);
disp('Done');
%%
trial_dataT = []; trial_dataD = [];
configurations = {'C3','C4','C5'};
air_phases = {'ON','OFF'};
for an = 1:5
    data_an = udata1{an};
    for cn = 1:length(configurations)
        for ap = 1:length(air_phases)
            [an,cn,ap]
            trial_dataT{an,cn,ap} = get_trial_metrics(data_an,'time',0.3,air_phases{ap},configurations{cn});
            trial_dataD{an,cn,ap} = get_trial_metrics(data_an,'distance',3,air_phases{ap},configurations{cn});
        end
    end
end
disp('Done');
%%
% % binwidths = [0.3 3];
% % for ii = 1:length(ei)
% %     ei(ii) = make_and_load_motion_correction(ei(ii),binwidths,[0 0 0]);
% % end
% 
% 
% %%
% % binwidths = [0.25 3];
% % ctl = [3 4 5]; sml = {'airOnsets55','airOffsets55'}; owr = [0 0 0];
% % for ii = 1:length(ei)
% %     ei(ii) = make_and_load_rasters_selected_contexts(ei(ii),binwidths,owr,ctl,sml);
% % end
% 
% %%
% % tic
% % for ii = 1:length(ei)
% %     ei(ii) = find_MI_MC(ei(ii),[0 0 0]);
% % end
% % toc
% %%
% % tic
% % binwidths = [0.3 3];
% % for ii = 1:length(ei)
% %     ei(ii) = make_and_load_rasters(ei(ii),binwidths,[0 0 0]);
% % end
% % toc
%%
clc
tic
for ii = 1:length(ei)
    ei(ii) = get_speed_response_Zs(ei(ii),[0 0]);
    ii
end
toc

tic
for ii = 1:length(ei)
    ei(ii) = get_accel_response_Zs(ei(ii),[0 0]);
end
toc

disp('Done');

%%
clc
tic
for ii = 1:length(ei)
    ei(ii) = get_motion_onset_response(ei(ii),[0 0 0 0 0]);
end
toc


tic
for ii = 1:length(ei)
    ei(ii) = get_speed_response(ei(ii),[0 0]);
end
toc

tic
for ii = 1:length(ei)
    ei(ii) = get_accel_response(ei(ii),[0 0]);
end
toc


tic
for ii = 1:length(ei)
    ei(ii) = get_speed_response_gauss(ei(ii),[0 0]);
end
toc

disp('Done');
%%
for ii = 1:length(ei)
    min_speed(ii) = min(ei{ii}.b.fSpeed);
    max_speed(ii) = max(ei{ii}.b.fSpeed);
end
%%
%%
clc
tic
for ii = 1:length(ei)
    ei(ii) = get_pca(ei(ii),binwidths,[-1 -1 -1]);
end
toc

%% cell_pose total cells

files = dir(sprintf('%s/mi*.mat',mData.pd_folder));
for ii = 1:length(files)
    tname = files(ii).name;
    anN = str2num(tname(9)); plN = str2num(tname(11));
    mi_cp{ii} = load(fullfile(files(ii).folder,files(ii).name));
    ei{anN}.plane{plN}.total_cells = length(unique(mi_cp{ii}.mask))-1;
    totalCells(ii) = ei{anN}.plane{plN}.total_cells;
end
totalCells(3) = totalCells(3) + 200; % I manually counted the cells as Cellpose didn't do a good job
totalCells(4) = totalCells(4) + 170; % I manually counted the cells as Cellpose didn't do a good job

cellposeCells(1) = totalCells(1);% + totalCells(2); just analyzing plane 1 for dist dur paper
cellposeCells(2) = totalCells(3);% + totalCells(4); just analyzing plane 1 for dist dur paper
cellposeCells(3) = totalCells(5);
cellposeCells(4) = totalCells(6);
cellposeCells(5) = totalCells(7);

disp('Done');

%% training

training_data = behaviorProcessor;
%%
training_data.belt_lengths = [150 150 150 150 150]';
training_data.DOB = {'2018-09-27';'2018-10-11';'2018-10-03';'2018-09-27';'2018-10-11'};
training_data.weight = [34.5000   34.2000   33.7000   33.5000   34.6
   31.6000   31.3000   30.7000   31.5000   30.3
   34.6000   34.5000   33.6000   35.3000   35.7
   26.8000   26.6000   26.7000   27.6000   26.6
   35.8000   35.4000   35.3000   35.1000   35.1];


animal_id_A = [183633,183761,183745,183628,183762];
date_of_rec_A = {'2019-06-04','2019-06-06','2019-06-07','2019-06-11','2019-06-11'};
date_of_surg_A = {'2019-03-27','2019-03-29','2019-04-30','2019-04-03','2019-04-26'};
for rr = 1:length(animal_id_A)
    ind = find(training_data.animalIDs == animal_id_A(rr));
    dob = datetime(training_data.DOB{ind},'InputFormat','yyyy-MM-dd');
%     this_date = datetime(date_of_rec_A{rr},'InputFormat','yyyy-MM-dd');
    this_date = datetime(training_data.training_dates{ind,4},'InputFormat','yyyy-MM-dd');
    dv = datevec(this_date-dob);
    ageAn_rec_A(rr) = 12*dv(1) + dv(2) + dv(3)/30;
    start_date = datetime(training_data.training_dates{ind,3},'InputFormat','yyyy-MM-dd');
    dv = datevec(this_date-start_date);
    rec_day_A(rr) = 12*dv(1) + dv(2) + dv(3)/30;
    surg_date = datetime(date_of_surg_A{rr},'InputFormat','yyyy-MM-dd');
    dv = datevec(start_date-surg_date);
    postsurg_day_A(rr) = 12*dv(1) + dv(2) + dv(3)/30;
    dv = datevec(surg_date - dob);
    age_at_surg(rr) = 12*dv(1) + dv(2) + dv(3)/30;
end
disp('Done');
%%
disp('Hello just executing this cell for testing ctrl enter');
%%