temp = load('\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_AD_paper\Control\183745\2019-06-19\1_002\avgstack.mat');
temp = load('\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data_AD_paper\APP\183329\2019-06-25\1\avgstack.mat')
avgimg = temp.avg_stack;
m = min(avgimg(:));
M = max(avgimg(:));
figure(600);clf;imagesc(avgimg,[0.7*m 0.5*M]);axis equal;axis off;colormap gray
