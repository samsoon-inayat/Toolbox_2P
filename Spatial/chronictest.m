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


rgb = cat(3,aframes{1}/max(aframes{1}(:)),aframes{2}/max(aframes{2}(:)));
rgb = cat(3,rgb,zeros(size(aframes{1})));
figure(105);imagesc(rgb);colorbar


fixed = aframes{1};%dicomread('knee1.dcm');
moving = aframes{2};%dicomread('knee2.dcm');
figure(102)
imshowpair(fixed, moving,'Scaling','joint')

[optimizer, metric] = imregconfig('multimodal');
movingRegistered = imregister(moving, fixed, 'affine', optimizer, metric);

figure(103)
imshowpair(fixed, movingRegistered,'Scaling','joint')

combined = cat(3,fixed,moving);
minV = min(combined(:));
maxV = max(combined(:));

figure(104);clf;
subplot 121; imagesc(fixed,[minV maxV]);
subplot 122; imagesc(moving,[minV maxV]);

%%
fixed = aframes{1};%dicomread('knee1.dcm');
moving = aframes{2};%dicomread('knee2.dcm');

[optimizer, metric] = imregconfig('multimodal');
movingRegistered = imregister(moving, fixed, 'affine', optimizer, metric);

rgb = cat(3,fixed/max(fixed(:)),moving/max(moving(:)));
rgb = cat(3,rgb,zeros(size(aframes{1})));
figure(105);imagesc(rgb);colorbar

rgb = cat(3,fixed/max(fixed(:)),movingRegistered/max(movingRegistered(:)));
rgb = cat(3,rgb,zeros(size(aframes{1})));
figure(106);imagesc(rgb);colorbar


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