function tag_videos(ei,ow)

for ii = 1:length(ei)
    tei = ei{ii};
    if isempty(tei.video_file)
        continue;
    end
    vfile = fullfile(tei.recordingFolder,tei.video_file_mp4);
    b = tei.b;
    b.air_puff_rawD = b.air_puff_raw > 0.5;
    b.stim_rawD = b.stim_raw > 0.5;
    vo = VideoReader(vfile);
    vts = b.video.times;
    cam_frame_inds = b.video.cam_frame_inds;
    st = b.air_puff_r(1)-1000;
    en = b.air_puff_f(2)+1000;
    fts = vts(cam_frame_inds > st & cam_frame_inds < en);
    [~,frames,~] = load_file_frames(vo,fts);
    cfii = cam_frame_inds(cam_frame_inds > st & cam_frame_inds < en);
    for fi = 1:length(fts)
        figure(100);clf;
        imagesc(frames{fi});
        title(fi);
        if b.air_puff_rawD(cfii(fi))
            text(100,100,'*','FontSize',20,'color','r');
        else
            text(100,100,'*','FontSize',20,'color','g');
        end
        pause(0.01);
    end
    figure(200);clf;plot(b.ts(cfii),b.air_puff_rawD(cfii))
    n = 0;
end

function [success,frames,frame_times] = load_file_frames(video_object,frameNums)
success = 1;
if isinteger(frameNums(1))

% video_object = VideoReader(fullfile(file_path,file_name));
number_of_frames = (ceil(video_object.FrameRate*video_object.Duration)-1);
frameNumber = 1;
for ii = 1:5
    t_frames{ii} = readFrame(video_object);
    t_frame_times(ii) = video_object.CurrentTime;
end
dt = t_frame_times(2)-t_frame_times(1);
startTime = dt * frameNums(1);
endTime = (dt * frameNums(end));
video_object.CurrentTime = startTime;
while video_object.CurrentTime <= endTime
    frame_times(frameNumber) = video_object.CurrentTime;
    try
        frames{frameNumber} = readFrame(video_object);
    catch
        break;
    end
    frameNumber = frameNumber + 1;
%     str = sprintf('File: %s ... loading frame %d/%d',file_name,(frameNumber),length(frameNums));
%     disp(str)
end
if ~success
    frames = [];
    frame_times = [];
    return;
end
else
    for ii = 1:length(frameNums)
        video_object.CurrentTime = frameNums(ii);
        frames{ii} = readFrame(video_object);
    end
    frame_times = frameNums;
end
