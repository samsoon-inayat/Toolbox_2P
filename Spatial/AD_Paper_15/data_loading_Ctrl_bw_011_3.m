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
processed_data_folder{1} = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_AD_paper_15\Control';
processed_data_folder{2} = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_AD_paper_15\Matlab\Control';
animal_list_control = {'174374';'183633';'183761';'183745';'183628';'183762'};
date_list_control = {'2019-02-22';'2019-06-04';'2019-06-06';'2019-06-07';'2019-06-11';'2019-06-11'};
[dS_C,T_C] = get_exp_info_from_folder(data_folder1,processed_data_folder,animal_list_control,date_list_control);
T_C = [T_C];
T_C_1 = reduce_table(T_C,animal_list_control,date_list_control);
[dS_C1,T_C1] = get_exp_info_from_folder(data_folder2,processed_data_folder,animal_list_control,date_list_control);
T_C1 = [T_C1];
T_C_2 = reduce_table(T_C1,animal_list_control,date_list_control);
T_Ca = [T_C_2;T_C_1];
disp('Done');

data_folder1 = '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data';
data_folder2 = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\Data';
processed_data_folder{1} = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_AD_paper_15\APP';
processed_data_folder{2} = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_AD_paper_15\Matlab\APP';
animal_list_control = {'001567';'001569';'183224';'183227';'183228';'183329'};
date_list_control = {'2020-04-27';'2020-04-28';'2019-03-08';'2019-03-08';'2019-06-26';'2019-05-10'};
[dS_A,T_A] = get_exp_info_from_folder(data_folder1,processed_data_folder,animal_list_control,date_list_control);
T_A = [T_A];
T_A_1 = reduce_table(T_A,animal_list_control,date_list_control);
[dS_A1,T_A1] = get_exp_info_from_folder(data_folder2,processed_data_folder,animal_list_control,date_list_control);
T_A1 = [T_A1];
T_A_2 = reduce_table(T_A1,animal_list_control,date_list_control);
T_Aa = [T_A_1;T_A_2];
disp('Done');

%%
colormaps = load('../../Common/Matlab/colorblind_colormap.mat');
colormaps.colorblind = flipud(colormaps.colorblind);
mData.colors = mat2cell(colormaps.colorblind,[ones(1,size(colormaps.colorblind,1))]);%{[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'}; % mData.colors = getColors(10,{'w','g'});
mData.axes_font_size = 6; mData.sigColor = [0.54 0.27 0.06];
mData.shades = generate_shades(3);
display_colors(mData.shades.c);
% Uleth_one_drive = 'Z:\homes\brendan.mcallister\2P';
Uleth_one_drive = 'E:\Users\samsoon.inayat\OneDrive - University of Lethbridge\PDFs';
% Uleth_one_drive = 'D:\OneDrive - University of Lethbridge\PDFs';
mData.pdf_folder = [Uleth_one_drive '\PDFs_AD_15']; 
mData.pd_folder = [Uleth_one_drive '\PDFs_AD_15\ProcessedDataMatlab'];
disp('Done');
%%
if 1
%     make_db(T_C);
    process_abf(T_Ca,0);
    process_abf(T_Aa,0);
end
disp('Done');

%%
ei_C = getData_py_2(T_Ca);
ei_A = getData_py_2(T_Aa);
ei_A(1) = [];
%%
ii = 1;
edit_define_contexts_file(ei_C{ii});
plotMarkers(ei_C{ii}.b,ei_C{ii}.b.stim_r,[],100,0)
plotMarkers(ei_C{ii}.b,ei_C{ii}.b.air_puff_r,[],101,0)
%%
ii = 1;
edit_define_contexts_file(ei_A{ii});
plotMarkers(ei_A{ii}.b,ei_A{ii}.b.stim_r,[],100,0)
plotMarkers(ei_A{ii}.b,ei_A{ii}.b.air_puff_r,[],101,0)
%%
binwidths = [0.11 3];
for ii = 1:length(ei_C)
    ei_C(ii) = make_and_load_rasters(ei_C(ii),binwidths,[0 0 0]);
end
for ii = 1:length(ei_A)
    ei_A(ii) = make_and_load_rasters(ei_A(ii),binwidths,[0 0 0]);
end
%%
clc
tic
for ii = 1:length(ei)
    ei(ii) = get_motion_onset_response(ei(ii),[0 0 0 0 0]);
end
toc


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