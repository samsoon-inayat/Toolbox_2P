% function post_processing
%%
add_to_path

clear all
clc
%%
[f,cName,ExpInfos] = getFolders;
%% load data

animal_id = '183761'; exp_date = '2019-06-21'; owr = 0; protocol = 10;
numberInList = find(strcmp(ExpInfos.animalList,animal_id));
animal_FS = ExpInfos.animal(numberInList)


plane_number = 1; ei{1} = getAllData_pyS2p(animal_FS,animal_id,exp_date,plane_number,'',owr); 
ei{1}.protocol = protocol;
ei{1}.ops1{1} = ei{1}.tP.ops;
ei{1}.b.belt_length = 150;
% plane_number = 2; ei{2} = getAllData(animal_id,exp_date,plane_number,'',owr); 

%% if contexts are already processed go to the step of loading contexts
% open this file and define contexts
behaviorPlot
edit('defineContexts.m');
%%
defineContexts(ei);
processContextDefinitions(ei);

%%
processContexts(ei);

%% load the contexts into ei variable
ei = loadContextsResponses(ei);

%% do gaussian fit on means (for now on distance rasters and for air trials only)
ei = gaussfitOnMeans(ei,'air');

%% find place cell properties
ei = placeCellProperties(ei,'air');

%% get data of responses only that is only the variable rasters from ei
[data,mData] = getRasterData(ei);
mData.colors = {[0 0 0],'b','r',[0 0.7 0.3],'m','c'};
mData.sigColor = [0.54 0.27 0.06];
