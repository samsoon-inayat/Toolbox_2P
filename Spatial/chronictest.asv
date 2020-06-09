clear all
rF{1} = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Data\001569\2020-05-19\1';
rF{2} = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Data\001569\2020-05-20\1';
rF{3} = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Data\001569\2020-05-21\1';

for ii = 1:length(rF)
    ii
    ei{ii} = thorGetExperimentInfo(rF{ii});
    frames{ii} = getFramesFromRaw(ei{ii},1:7000);
end

for ii = 1:length(rF)
    aframes{ii} = mean(frames{ii},3);
    mframes{ii} = max(frames{ii},[],3);
end


figure(100);clf;
subplot 131;imagesc(aframes{1});axis equal;axis off;
subplot 132;imagesc(aframes{2});axis equal;axis off;
subplot 133;imagesc(aframes{3});axis equal;axis off;


figure(101);clf;
subplot 131;imagesc(mframes{1});axis equal;axis off;
subplot 132;imagesc(mframes{2});axis equal;axis off;
subplot 133;imagesc(mframes{3});axis equal;axis off;

% for ii = 1:length(rF)
%     imwrite(aframes{ii},sprintf('avg_chronic_%d.tif',ii));
%     imwrite(mframes{ii},sprintf('max_chronic_%d.tif',ii));
% end




% 
% % 
% 
% for ii = 1:size(frames{1},3)
%     figure(100);clf;
%     imagesc(frames{1}(:,:,ii));
%     pause(0.1);
% end