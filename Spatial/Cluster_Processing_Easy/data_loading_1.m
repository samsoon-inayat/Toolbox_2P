for ii = 1:60
    pause(1);
end
%%
% clear all
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
processed_data_folder = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_Basic_Char';
f.data_folder = data_folder; f.processed_data_folder = processed_data_folder;
[dS,T] = get_exp_info_from_folder(data_folder,processed_data_folder);
% sT = T([2 11 20 3 12 21 4 13 22 5 14 23 6 15 24],:);
sT = T;
%%
if 0
    make_db(T);
    process_abf(T,0);
end
%%
ei = getData_py_1(f,sT,0);
%%
for ii = 1:length(ei)
%     process_abf(sT(ii,:),0);
    ei(ii) = loadContextsResponses_ctrl(ei(ii),[1 1],[0 0 0]);
end
% parameter_matrices_ctrl('calculate','P2_1',ei(1:3));
%%

selContexts = [3 3];
rasterNames = {'airD','beltD'};
Rs = get_rasters_data(ei,selContexts,rasterNames);

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



