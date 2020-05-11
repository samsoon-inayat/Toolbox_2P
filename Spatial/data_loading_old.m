% function post_processing
%%
add_to_path
%% load data table
clear all
clc
[f,cName] = getFolders;
load('T.mat');
animal_ids = cell2mat(T{:,1});
%%
ow = 0;
% ei = getData(f,T(18,:));
eip = getData_py(f,T([3 7 11 13 15 16 18 20 21],:));
disp('Done');

%%
sei = [];
sei = eip;%([1 3 5 6]);
%% if contexts are already processed go to the step of loading contexts
% open this file and define contexts

selR = 6;
behaviorPlot(sei(1:selR))

%% load the contexts into ei variable
sei = loadContextsResponses(sei);
disp('Done!!!');
%% do gaussian fit on means (for now on distance rasters and for air trials only)
sei = gaussfitOnMeans(sei,'air');

%% find place cell properties
sei = placeCellProperties(sei,'air');

%% get data of responses only that is only the variable rasters from ei
[data,mData] = getRasterData(sei,[1 2 3 4],'air');
mData.colors = {[0 0 0],'b','r',[0 0.7 0.3],'m','c'};
mData.sigColor = [0.54 0.27 0.06];
