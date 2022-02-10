function ei = load_video_frame_inds(ei,ow)

for ii = 1:length(ei)
    tei = ei{ii};
    if isempty(tei.video_file)
        continue;
    end
    file_name = fullfile(tei.matlab_folder,'camera_frame_indices.mat');
    if exist(file_name,'file') && ow == 0
        temp = load(file_name);
        tei.b.video.cam_frame_inds = temp.cam_frame_inds;
        tei.b.video.times = temp.vts;
        ei{ii} = tei;
        continue;
    end
    disp('Processing video indexing');
    vfile = fullfile(tei.recordingFolder,tei.video_file_mp4);
    b = tei.b;
    b.air_puff_rawD = b.air_puff_raw > 0.5;
    b.stim_rawD = b.stim_raw > 0.5;
    vo = VideoReader(vfile);
    number_of_frames = (ceil(vo.FrameRate*vo.Duration)-1);
    [success,frames,frame_times] = load_file_frames(vo,1:100);
    dt = frame_times(2) - frame_times(1);
    vts = frame_times(1):dt:vo.Duration;
    bts = b.ts;
    cam_frame_inds = NaN(size(vts));
    parfor vi = 1:length(vts)
        cam_frame_inds(vi) = find(bts - vts(vi) < 0,1,'last');
    end
%     figure(1000);clf;
%     plot(b.ts,b.air_puff_raw);hold on;
%     plot(b.ts(cam_frame_inds),b.air_puff_rawD(cam_frame_inds),'r');
    save(file_name,'cam_frame_inds','vts');
    tei.b.video.cam_frame_inds = cam_frame_inds;
    tei.b.video.times = vts;
    ei{ii} = tei;
    n = 0;
end


function [success,frames,frame_times] = load_file_frames(video_object,frameNums)

success = 1;
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
