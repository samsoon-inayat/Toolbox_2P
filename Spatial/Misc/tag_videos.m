function tag_videos(ei,T,ow)

for ii = 1:length(ei)
    if isempty(cell2mat(T{ii,10}))
        continue;
    end
    tei = ei{ii};
    vfile = fullfile(tei.recordingFolder,cell2mat(T{ii,10}));
    vfile(end) = [];
    vfile((length(vfile)-2):end) = 'mp4';
    n = 0;
    b = tei.b;
    vo = VideoReader(vfile);
    number_of_frames = (ceil(vo.FrameRate*vo.Duration)-1);
    [success,frames,frame_times] = load_file_frames(vo,1:10);
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
