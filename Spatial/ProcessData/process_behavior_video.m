function process_behavior_video

file_path = '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data\001569\2020-04-25'
file_name = 'mouse_001569_20-04-25_17_43_57_vid_60_fps.mp4';

fileName = fullfile(file_path,file_name);

v = VideoReader(fileName);

n = 0;
%%
roi{1} = [150,350,100,550];
thisROI = roi{1};
fns = 1050:8550;
for ii = 1:length(fns)
    frame = read(v,fns(ii));
    trig_detect(ii) = mean(frame(1:50,:),'All');
    rows = thisROI(1):thisROI(2); cols = thisROI(3):thisROI(4);
    figure(1000);clf;imagesc(frame(rows,cols));colorbar
    title(fns(ii));
    pause(0.1);
end
%%
figure(1001);clf;
plot(trig_detect)