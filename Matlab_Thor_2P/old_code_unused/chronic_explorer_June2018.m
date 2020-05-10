%%
add_to_path
%% load data
clear all
clc

owr = 0;

ei{1} = getAllData('172653','2018-06-21',1,''); ei{2} = getAllData('172653','2018-06-21',2,'');
forFindingTrials(ei{1},27:42);
trials = {[1:10],[11:26],[27:42]};
recs = [1 2];
for ii = 1:length(trials)
    temp = getRasters(ei(recs),trials{ii},owr);
    temp1 = find_mutual_information(ei(recs),temp,0);
end
temp = getRastersAir(ei(recs),trials{1},owr);
temp1 = find_mutual_information(ei(recs),temp,0);


ei{3} = getAllData('172649','2018-06-21',1,''); ei{4} = getAllData('172649','2018-06-21',2,'');
forFindingTrials(ei{3},27:42)
recs = [3 4];
for ii = 1:length(trials)
    temp = getRasters(ei(recs),trials{ii},owr);
    temp1 = find_mutual_information(ei(recs),temp,0);
end
temp = getRastersAir(ei(recs),trials{1},owr);
temp1 = find_mutual_information(ei(recs),temp,0);

ei{5} = getAllData('173062','2018-06-21',1,''); ei{6} = getAllData('173062','2018-06-21',2,'');
forFindingTrials(ei{5},1:10)
recs = [5 6];
for ii = 1:length(trials)
    temp = getRasters(ei(recs),trials{ii},owr);
    temp1 = find_mutual_information(ei(recs),temp,0);
end
temp = getRastersAir(ei(recs),trials{1},owr);
temp1 = find_mutual_information(ei(recs),temp,0);



ei{7} = getAllData('173062','2018-06-20',1,''); ei{8} = getAllData('173062','2018-06-20',2,'');
trials = {[1:15],[16:31],[32:47]};
recs = [7 8];
for ii = 1:length(trials)
    temp = getRasters(ei(recs),trials{ii},owr);
    temp1 = find_mutual_information(ei(recs),temp,0);
end
temp = getRastersAir(ei(recs),trials{1},owr);
temp1 = find_mutual_information(ei(recs),temp,0);

% 
% 
% ei{1} = getAllData('172653','2018-06-16',1,'');
% ei{2} = getAllData('172653','2018-06-16',2,'');
% ei{3} = getAllData('172649','2018-06-16',1,'');
% ei{4} = getAllData('172649','2018-06-16',2,'');
% ei{5} = getAllData('173062','2018-06-16',1,'');
% ei{6} = getAllData('173062','2018-06-16',2,'');
% 
% 
% 
% owr = 1;
% recs = [1 2];
% temp = getRasters(ei(recs),1:11,owr);
% temp = getRasters(ei(recs),12:23,owr);
% recs = [3 4];
% temp = getRasters(ei(recs),1:11,owr);
% temp = getRasters(ei(recs),12:23,owr);
% recs = [5 6];
% temp = getRasters(ei(recs),1:11,owr);
% temp = getRasters(ei(recs),12:23,owr);
% 
% 
% 
% 
% 
% 
% 
