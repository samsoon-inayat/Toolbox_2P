for ii = 1:60
    pause(1);
end
clear all
%%
% add_to_path
% clc
[f,cName] = getFolders;
T10 = load('T_10_All.mat');
selRecs10 = [8 9 10 11];
ET10 = T10.T(selRecs10,:);
disp('done')
%%
colormaps = load('../MatlabCode/colorblind_colormap.mat');
colormaps.colorblind = flipud(colormaps.colorblind);
% mData.colors = mat2cell(colormaps.colorblind,[ones(1,size(colormaps.colorblind,1))]);%
mData.colors = {[0 0 0],[0.1 0.7 0.3],'r','b','m','c','g','y'}; % mData.colors = getColors(10,{'w','g'});
mData.axes_font_size = 6; mData.sigColor = [0.54 0.27 0.06]; 
mData.pdf_folder = 'E:\Users\samsoon.inayat\OneDrive - University of Lethbridge\PDFs'; 
mData.pd_folder = 'E:\Users\samsoon.inayat\OneDrive - University of Lethbridge\ProcessedData\Matlab';
disp('Done');
%%
for ii = 1:size(ET10,1)
    data(ii) = getData_py(f,ET10(ii,:));
end

data(4) = loadContextsResponses_ctrl(data(4),[1 1],[0 0 0]);
% training_data_C1 = behaviorProcessor;
% training_data_C1.belt_lengths = [150 142 142 142 142 142 150 150 150 150 150]';
% weight_day1 = [52.1 NaN NaN NaN 33.7 30.9 26.8 34.5 34.6 31.6 35.8]';
% weight_day2 = [52.8 NaN NaN NaN 33.3 30.4 26.6 34.2 34.5 31.3 35.4]';
% weight_day3 = [51.6 NaN NaN NaN 33.3 30.5 26.7 33.7 33.6 30.7 35.3]';
% training_data_C1.weight = [weight_day1 weight_day2 weight_day3];
parameter_matrices_ctrl('calculate','10_CD_Ctrl',ei10_C);

disp('All Done!');
%%
