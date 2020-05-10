%%
add_to_path
%% load data
clear all
clc
% Thy1-GCaMP
ei{1} = getAllData('171179','2018-02-16',1,'RHL_4');%RSC
ei{2} = getAllData('171179','2018-02-16',2,'RHL_4');%RSC
ei{3} = getAllData('172110','2018-02-13',1,'LHL_6');
ei{4} = getAllData('172110','2018-02-13',2,'LHL_6');

% AD Thy1-GCaMP
ei{5} = getAllData('171465','2018-04-05',1,'RH_2');
ei{6} = getAllData('171465','2018-04-05',2,'RH_2');
ei{7} = getAllData('171466','2018-04-05',1,'RH_2');
ei{8} = getAllData('171466','2018-04-05',2,'RH_2');

%%
tei = ei{7};
onsets = tei.b.air_puff_r;
offsets = tei.b.air_puff_f;
plotMarkers(tei,onsets,offsets,101);


%% get rasters and MI
ei{1} = getRastersAndMI(ei{1},13:23,0,0);
ei{2} = getRastersAndMI(ei{2},13:23,0,0);
ei{3} = getRastersAndMI(ei{3},14:23,0,0);
ei{4} = getRastersAndMI(ei{4},14:23,0,0);
ei{5} = getRastersAndMI(ei{5},[1:11,13:15,17:20],0,0);
ei{6} = getRastersAndMI(ei{6},[1:11,13:15,17:20],0,0);
ei{7} = getRastersAndMI(ei{7},[1:5,8:15,17:20],0,0);
ei{8} = getRastersAndMI(ei{8},[1:5,8:15,17:20],0,0);
