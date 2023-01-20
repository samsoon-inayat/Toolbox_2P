for ii = 1:60
    pause(1);
end
%%
add_to_path
%%
clear all
% clc
%%
data_folder1 = 'E:\Data';%'\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data';
data_folder2 = 'E:\Data';%'\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\Data';
processed_data_folder{1} = 'E:\PData\Processed_Data_15';%'\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_15';
processed_data_folder{2} = 'E:\PData\Processed_Data_15\Matlab';%'\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_15\Matlab';
processed_data_folder{3} = 'E:\PData\Processed_Data_15\Matlab_bw3';%'\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_15\Matlab_bw3';
animal_list_control = {'183633';'183761';'183745';'183628';'183762'};
date_list_control = {'2019-06-04';'2019-06-06';'2019-06-07';'2019-06-11';'2019-06-11'};
[dS_C,T_C] = get_exp_info_from_folder(data_folder1,processed_data_folder,animal_list_control,date_list_control);
T_C = [T_C];
T_C1 = reduce_table(T_C,animal_list_control,date_list_control);
disp('Done');
%% check for video files
% T_C1 = check_for_video_files(T_C1);
%%
colormaps = load('../../Common/Matlab/colorblind_colormap.mat');
colormaps.colorblind = flipud(colormaps.colorblind);
mData.colors = mat2cell(colormaps.colorblind,[ones(1,size(colormaps.colorblind,1))]);%{[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'}; % mData.colors = getColors(10,{'w','g'});
mData.dcolors = mat2cell(distinguishable_colors(20,'w'),[ones(1,20)]);
mData.axes_font_size = 6; mData.sigColor = [0.54 0.27 0.06];
mData.shades = generate_shades(3);
mData.conj_comp_colors = [mData.dcolors(9);mData.colors([3 5])];
% display_colors(mData.colors);
% Uleth_one_drive = 'Z:\homes\brendan.mcallister\2P';
% Uleth_one_drive = 'E:\Users\samsoon.inayat\OneDrive - University of Lethbridge\PDFs';
Uleth_one_drive = 'E:\PostProcessing\Revamp_Conj_Comp_Paper';
mData.pdf_folder = [Uleth_one_drive '\PDFs']; 
mData.pd_folder = [Uleth_one_drive '\PD_Matlab'];
mData.magfac = 1;
disp('Done');
%%
if 0
%     make_db(T_C);
    process_abf(T_C,0);
end
disp('Done');
%%
sel_rec = [1 3 5];
sel_rec = 1:5;
ei = getData_py_2(T_C1(sel_rec,:));

%% tag videos with event related signals
if 0
    ei = check_for_video_files_ei(ei);
    ei = load_video_frame_inds(ei,1);
    tag_videos(ei,0);
end
%%
if 0
    ii = 5;
    edit_define_contexts_file(ei{ii});
end
dcfilename = 'define_contexts_revamp.m';
%% get control light onsets
for ii = 1:length(ei)
    ei(ii) = get_control_light_onsets(ei(ii));
    ei(ii) = get_control_air_onsets(ei(ii));
    ei(ii) = get_control_air_onsets_C(ei(ii));
    ei(ii) = get_control_air_offsets(ei(ii));
end
send_email('samsoon.inayat@gmail.com','Done');
%% get control light onsets
for ii = 1:length(ei)
    ei(ii) = get_control_air_Monsets(ei(ii));
    ei(ii) = get_control_air_Moffsets(ei(ii));
end
send_email('samsoon.inayat@gmail.com','Done');
%%
binwidths = [0.11 3];
for ii = 1:length(ei)
    ei(ii) = make_and_load_motion_correction(ei(ii),binwidths,[0 0 0]);
end
send_email('samsoon.inayat@gmail.com','Done');
%%
binwidths = [0.11 3];
for ii = 1:length(ei)
    ei(ii) = load_context_info(ei(ii),binwidths,[0 0 0],dcfilename);
    ei(ii) = load_motion_correction_empty(ei(ii),binwidths,[0 0 0]);
end
disp('Done');
%%
tic
for ii = 1:length(ei)
    ei(ii) = find_MI_MC(ei(ii),[0 0 0]);
end
toc
%%
clc
binwidths = [0.11 3];
tic
for ii = 1:length(ei)
    ei(ii) = make_and_load_rasters(ei(ii),binwidths,[0 0 0]);
end

for ii = 1:length(ei)
    ei(ii) = get_motion_onset_response(ei(ii),[0 0 0 0 0]);
end
toc
disp('Done');
send_email('samsoon.inayat@gmail.com','Done');
%% cell_pose total cells

files = dir(sprintf('%s/*.mat',mData.pd_folder));
for ii = 1:length(files)
    tname = files(ii).name;
    anN = str2num(tname(9)); plN = str2num(tname(11));
    mi_cp{ii} = load(fullfile(files(ii).folder,files(ii).name));
    ei{anN}.plane{plN}.total_cells = length(unique(mi_cp{ii}.mask))-1;
    totalCells(ii) = ei{anN}.plane{plN}.total_cells;
end
totalCells(3) = totalCells(3) + 200; % I manually counted the cells as Cellpose didn't do a good job
totalCells(4) = totalCells(4) + 170; % I manually counted the cells as Cellpose didn't do a good job

cellposeCells(1) = totalCells(1) + totalCells(2);
cellposeCells(2) = totalCells(3) + totalCells(4);
cellposeCells(3) = totalCells(5);
cellposeCells(4) = totalCells(6);
cellposeCells(5) = totalCells(7);

% for ii = 1:length(files)
%     maskii = mi_cp{ii}.mask;
%     totC = length(unique(maskii))-1;
% %     if ii == 1
% %         ei{1}.plane{1}.totalCells = 
% %     end
% %     if ii == 2
% %         ei{1}.plane{2}.totalCells = length(unique(maskii))-1;
% %     end
% %     rois = regionprops(maskii);
% %     figure(10000);clf;imagesc(maskii);axis equal;colorbar
% %     totalcells(ii) = length(unique(maskii))-1;
% end

%%
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

%% this information is already in the previous function
% tic
% for ii = 1:length(ei)
%     ei(ii) = get_accel_response_gauss(ei(ii),[0 1]);
% end
% toc


%%
for ii = 1:length(ei)
    min_speed(ii) = min(ei{ii}.b.fSpeed);
    max_speed(ii) = max(ei{ii}.b.fSpeed);
    
    t_accel = diff(ei{ii}.b.fSpeed)./diff(ei{ii}.b.ts);
    samplingRate = floor(1/(ei{ii}.b.si*1e-6));
    coeffs = ones(1, samplingRate)/samplingRate;
    ft_accel = filter(coeffs, 1, t_accel);
    min_accel(ii) = min(ft_accel);
    max_accel(ii) = max(ft_accel);
end
%%
%%
clc
tic
for ii = 1:length(ei)
    ei(ii) = get_pca(ei(ii),binwidths,[-1 -1 -1]);
end
toc


%% training

training_data = behaviorProcessor;
training_data.belt_lengths = [150 150 150 150 150]';
training_data.DOB = {'2018-09-27';'2018-10-11';'2018-10-03';'2018-09-27';'2018-10-11'};
training_data.weight = [34.5000   34.2000   33.7000   33.5000   34.6
   31.6000   31.3000   30.7000   31.5000   30.3
   34.6000   34.5000   33.6000   35.3000   35.7
   26.8000   26.6000   26.7000   27.6000   26.6
   35.8000   35.4000   35.3000   35.1000   35.1];



animal_id_A = [183633,183761,183745,183628,183762];
gender_animals = {'M','M','F','F','M'};
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