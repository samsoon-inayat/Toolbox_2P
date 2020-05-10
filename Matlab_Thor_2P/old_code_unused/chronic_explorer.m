%%
add_to_path
%% load data
clear all
clc

ei{1} = getAllData('171465','2018-04-05',1,'RH_2');
ei{2} = getAllData('171465','2018-04-05',2,'RH_2');

ei{3} = getAllData('171466','2018-04-05',1,'RH_2');
ei{4} = getAllData('171466','2018-04-05',2,'RH_2');

ei{1} = getAllData('171179','2018-03-27',1,'RH_9');
ei{2} = getAllData('171179','2018-03-27',2,'RH_9');

ei{1} = getAllData('171179','2018-02-13',1,'RHL_2');%M2
ei{2} = getAllData('171179','2018-02-13',2,'RHL_2');%M2
ei{3} = getAllData('171179','2018-02-16',1,'RHL_4');%RSC
ei{4} = getAllData('171179','2018-02-16',2,'RHL_4');%RSC
ei{5} = getAllData('172110','2018-02-13',1,'LHL_6');
ei{6} = getAllData('172110','2018-02-13',2,'LHL_6');
ei{7} = getAllData('172376','2018-02-09',1,'Right_Hemisphere_Location3');
ei{8} = getAllData('172376','2018-02-09',2,'Right_Hemisphere_Location3');

owr = 1;
recs = [1 2];
temp = getRasters(ei(recs),1:11,owr);
temp = getRasters(ei(recs),13:23,owr);
recs = [3 4];
temp = getRasters(ei(recs),1:11,owr);
temp = getRasters(ei(recs),13:23,owr);
recs = [5 6];
temp = getRasters(ei(recs),1:11,owr);
temp = getRasters(ei(recs),13:23,owr);


%% looking at tone light data
ei{1} = getAllData('172110','2018-01-30',1,'Left_Hemisphere_Location1',1);





