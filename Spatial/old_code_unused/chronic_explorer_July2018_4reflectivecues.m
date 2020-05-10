%%
add_to_path
%% load data
clear all
clc

owr = 0;

ei{1} = getAllData('172653','2018-06-22',1,''); ei{2} = getAllData('172653','2018-06-22',2,'');
trials = {[1:15],[16:31],[32:47]};
recs = [1 2];
for ii = 1:length(trials)
    temp = getRasters(ei(recs),trials{ii},owr);
    temp1 = find_mutual_information(ei(recs),temp,0);
end
temp = getRastersAir(ei(recs),trials{1},owr);
temp1 = find_mutual_information(ei(recs),temp,0);

ei{3} = getAllData('172649','2018-06-22',1,''); ei{4} = getAllData('172649','2018-06-22',2,'');
recs = [3 4];
for ii = 1:length(trials)
    temp = getRasters(ei(recs),trials{ii},owr);
    temp1 = find_mutual_information(ei(recs),temp,0);
end
temp = getRastersAir(ei(recs),trials{1},owr);
temp1 = find_mutual_information(ei(recs),temp,0);

ei{1} = getAllData('173062','2018-06-22',1,''); ei{2} = getAllData('173062','2018-06-22',2,'');

ei{3} = getAllData('173062','2018-06-28',1,''); ei{4} = getAllData('173062','2018-06-28',2,'');
recs = [5 6];
for ii = 1:length(trials)
    temp = getRasters(ei(recs),trials{ii},owr);
    temp1 = find_mutual_information(ei(recs),temp,0);
end
temp = getRastersAir(ei(recs),trials{1},owr);
temp1 = find_mutual_information(ei(recs),temp,0);

ei{5} = getAllData('173062','2018-07-02',1,'L1'); ei{6} = getAllData('173062','2018-07-02',2,'L1');
ei{7} = getAllData('173062','2018-07-03',1,'L1'); ei{8} = getAllData('173062','2018-07-03',2,'L1');
recs = [7 8];
trials = {[1:13],[13:25],[26:38]};
ii = 1;
temp = getRasters(ei(recs),trials{ii},owr);
temp1 = find_mutual_information(ei(recs),temp,0);
for ii = 2:length(trials)
    temp = getRasters_AL(ei(recs),trials{ii},owr);
    temp1 = find_mutual_information(ei(recs),temp,0);
end
temp = getRastersAir(ei(recs),trials{1},owr);
temp1 = find_mutual_information(ei(recs),temp,0);

ei{1} = getAllData('173062','2018-07-05',1,'L2'); ei{2} = getAllData('173062','2018-07-05',2,'L2');
trials = {[1:15],[16:31],[32:47]};
recs = [1 2];
for ii = 1:length(trials)
    temp = getRasters(ei(recs),trials{ii},owr);
    temp1 = find_mutual_information(ei(recs),temp,0);
end
temp = getRastersAir(ei(recs),trials{1},owr);
temp1 = find_mutual_information(ei(recs),temp,0);
