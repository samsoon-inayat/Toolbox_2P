for ii = 1:60
    pause(1);
end
%%
clear all
add_to_path
% clc
%%
[f,cName] = getFolders;
colormaps = load('../MatlabCode/colorblind_colormap.mat');
mData.colors = mat2cell(colormaps.colorblind([5 7 8 6 1 2 11],:),[ones(1,7)]);%
mData.line_styles = {'-',':','-.','-',':','-.'};
mData.marker_styles = {'s','d','^','s','d','^','^','*'};
mData.axes_font_size = 6; mData.sigColor = [0.54 0.27 0.06]; 
Uleth_one_drive = 'G:\OneDrives\OneDrive - University of Lethbridge';
Uleth_one_drive = 'E:\Users\samsoon.inayat\OneDrive - University of Lethbridge';
% Uleth_one_drive = 'D:\Dropbox\OneDrive\Documents\Manuscripts\Air_Puff_Characterization';
mData.pdf_folder = [Uleth_one_drive '\PDFs']; 
mData.pd_folder = [Uleth_one_drive '\ProcessedData'];
disp('Done');
%% Protocol 15
temp = load('T_15_All.mat');
T15 = temp.T;
T15_old = T15;
T15_c = T15;
sel15 = [2 4 6 8 12];
new_pd_folder = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_15';
for ii = 1:length(sel15)
    tempF = T15_old{sel15(ii),7}{1};
    pos = strfind(tempF,'Processed_Data');
    T15_c{sel15(ii),7} = {[new_pd_folder '\' tempF((pos+length('Processed_Data\')):end)]};
end
% for ii = 1:length(sel15)
%     process_abf(T15(sel15(ii),:));
% end
disp('done')
%%
ET15 = T15_c(sel15(1:5),:);
d15 = getData_py(f,T15_c(sel15(1:5),:),0);
binWidths = [0.15,1.75];
% d15_1 = load_contexts_responses(d15,'define_contexts.m',binWidths);
d15 = loadContextsResponses_ctrl(d15,[1 1],[0 0 0]);
selContexts = [1 2 3 3 4 4 5 5 6 7];
rasterNames = {'light22T','air55T','air77T','airD','air77T','airD','air77T','airD','light22T','air55T'};
raster_data = get_rasters_data(d15,selContexts,rasterNames);
parameter_matrices_ctrl('calculate','15_Charac',d15);
% testing why the number of columns for animal 1 and 3 are different
% d15_1 = loadContextsResponses_ctrl(d15(3),[1 1],[0 0 0]);
% d15{5} = loadContextsResponses_ctrl(d15(5),[1 1],[0 0 0]);

d15 = loadContextsResponses_ctrl(d15,[1 1],[0 0 0]);

%%


%%
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