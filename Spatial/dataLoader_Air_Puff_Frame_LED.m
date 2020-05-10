%%
add_to_path
%% load data
clear all
clc

owr = 0;
ei{1} = getAllData('183633','2019-06-04',1,'',owr); ei{2} = getAllData('183633','2019-06-04',2,'',owr);


for ii = [1:2]
    ei(ii) = loadContextsResponses(ei(ii));
end

%% View Data
sei = ei{1};
ccs = sei.areCells;
plotRastersMulti(data,ccs,0,1,1)

mean_raster_fits{ii} = getFits(data{1},ccs,3:10);
%%

[data,mData] = loadDataFrom_ei_light(ei(1));
mData.colors = {[0 0 0],'b','r',[0 0.7 0.3],'m','c'};
mData.sigColor = [0.54 0.27 0.06];
% mData.colors = {'b','r',([35 142 35]/255),

[data1,mData1] = loadDataFrom_ei(ei(7:8));
data = data1;
mData = mData1;
mData.colors = {[0 0 0],'b','r',[0 0.7 0.3],'m','c'};
mData.sigColor = [0.54 0.27 0.06];
