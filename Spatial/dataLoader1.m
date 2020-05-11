%%
add_to_path
%% load data
clear all
clc

owr = 0;
ei{1} = getAllData('173062','2018-07-11',1,''); ei{2} = getAllData('173062','2018-07-11',2,'');
ei{3} = getAllData('173511','2018-07-11',1,''); ei{4} = getAllData('173511','2018-07-11',2,'');
ei{5} = getAllData('173198','2018-07-11',1,''); ei{6} = getAllData('173198','2018-07-11',2,'');
ei{7} = getAllData('174374','2019-02-12',1,''); ei{8} = getAllData('174374','2019-02-12',2,'');
ei{9} = getAllData('173706','2019-02-14',1,''); ei{10} = getAllData('173706','2019-02-14',2,'');

for ii = [1:10]
    ei(ii) = loadContextsResponses(ei(ii));
end

[data,mData] = loadDataFrom_ei(ei(1:6));
mData.colors = {[0 0 0],'b','r',[0 0.7 0.3],'m','c'};
mData.sigColor = [0.54 0.27 0.06];
% mData.colors = {'b','r',([35 142 35]/255),

[data1,mData1] = loadDataFrom_ei(ei(7:8));
data = data1;
mData = mData1;
mData.colors = {[0 0 0],'b','r',[0 0.7 0.3],'m','c'};
mData.sigColor = [0.54 0.27 0.06];
