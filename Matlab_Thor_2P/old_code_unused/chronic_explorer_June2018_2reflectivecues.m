%%
add_to_path
%% load data
clear all
clc

owr = 0;

ei{1} = getAllData('172653','2018-06-25',1,''); ei{2} = getAllData('172653','2018-06-25',2,'');
trials = {[1:15],[16:31],[32:47]};
recs = [1 2];
for ii = 1:length(trials)
    temp = getRasters(ei(recs),trials{ii},owr);
    temp1 = find_mutual_information(ei(recs),temp,0);
end
temp = getRastersAir(ei(recs),trials{1},owr);
temp1 = find_mutual_information(ei(recs),temp,0);

ei{3} = getAllData('172649','2018-06-25',1,''); ei{4} = getAllData('172649','2018-06-25',2,'');
recs = [3 4];
for ii = 1:length(trials)
    temp = getRasters(ei(recs),trials{ii},owr);
    temp1 = find_mutual_information(ei(recs),temp,0);
end
temp = getRastersAir(ei(recs),trials{1},owr);
temp1 = find_mutual_information(ei(recs),temp,0);

ei{5} = getAllData('173062','2018-06-25',1,''); ei{6} = getAllData('173062','2018-06-25',2,'');
recs = [5 6];
for ii = 1:length(trials)
    temp = getRasters(ei(recs),trials{ii},owr);
    temp1 = find_mutual_information(ei(recs),temp,0);
end
temp = getRastersAir(ei(recs),trials{1},owr);
temp1 = find_mutual_information(ei(recs),temp,0);


