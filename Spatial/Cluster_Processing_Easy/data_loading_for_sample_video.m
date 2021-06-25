for ii = 1:60
    pause(1);
end
%%
clear all
%%
% add_to_path
clc
%%
[f,cName] = getFolders;
colormaps = load('../MatlabCode/colorblind_colormap.mat');
colormaps.colorblind = flipud(colormaps.colorblind);
% mData.colors = mat2cell(colormaps.colorblind,[ones(1,size(colormaps.colorblind,1))]);%
mData.colors = {[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'}; % mData.colors = getColors(10,{'w','g'});
mData.axes_font_size = 6; mData.sigColor = [0.54 0.27 0.06]; 
Uleth_one_drive = 'G:\OneDrives\OneDrive - University of Lethbridge';
Uleth_one_drive = 'E:\Users\samsoon.inayat\OneDrive - University of Lethbridge';
mData.pdf_folder = [Uleth_one_drive '\PDFs\P15']; 
mData.pd_folder = [Uleth_one_drive '\ProcessedData\Matlab'];
disp('Done');
%%
data_folder = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Data\RSEG_PSEG_more_data';
processed_data_folder{1} = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_Basic_Char';
processed_data_folder{2} = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_Basic_Char\Matlab';
f.data_folder = data_folder; f.processed_data_folder = processed_data_folder;
% [dS,T] = get_exp_info_from_folder(data_folder,processed_data_folder);
% sT = T([2 11 20 3 12 21 4 13 22 5 14 23 6 15 24],:);
animal_list_control = {'3328'};
date_list_control = {'2021-05-11'};
[dS_C,T_C] = get_exp_info_from_folder(data_folder,processed_data_folder,animal_list_control,date_list_control);
T_C = [T_C];
T_C = reduce_table(T_C,animal_list_control,date_list_control);
disp('Done');
%%
if 1
%     make_db(sT);
    process_abf(T_C,0);
end
%%
ei = getData_py_2(T_C);
%%
for ii = 1%:length(ei)
    process_abf(sT(ii,:),0);
%     ei(ii) = loadContextsResponses_ctrl(ei(ii),[1 1],[1 1 1]);
end

%%
binwidths = [0.2 1.5];
for ii = 1:length(ei)
    ei(ii) = make_and_load_rasters(ei(ii),binwidths);
end

%%
video_file = '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data\Videos\mouse_3328_21-05-11_16_13_59_vid_60_fps.mp4';
csv_file = '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data\Videos\mouse_3328_21-05-11_16_13_59_TrigTimes_60_fps.csv';
[raw] = readtable(csv_file);
trigger_frame = raw{1,end};
vo = VideoReader(video_file);

%%
rr = 1; cc = 3;
recording(1).animal_id = dS.exp_list_animal{rr,cc};
recording(1).date = dS.exp_list_date{rr,cc};
recording(1).exp_dir = dS.exp_list_expDir{rr,cc};
pd_rec = load_processed_data(recording(1),data_folder,processed_data_folder);
disp('Done');
%%
trace_plot_all(pd_rec);
%%
rasterNames = {'airD','beltD'};
Rs = get_rasters_data_from_rasterNames(ei,rasterNames);

%%
Rs.rasters = raster_data{1,1}.sp_rasters_nan_corrected1;
plotRasters_simple(Rs,1:350,[])
%%
onsets = pd_rec.b.photo_sensor_f(1:16);
offsets = pd_rec.b.photo_sensor_f(2:17);
rastersP =  getRasters_fixed_bin_width_cluster(pd_rec,onsets,offsets,'dist');
%%
RsP.rasters = rastersP.sp_rasters_nan_corrected1;
plotRasters_simple(RsP,1:350,[])



