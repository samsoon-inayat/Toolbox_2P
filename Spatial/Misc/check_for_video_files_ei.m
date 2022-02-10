function ei = check_for_video_files(ei)


for ii = 1:length(ei)
    tei = ei{ii};
    dataFolder = tei.recordingFolder;
    files_list = dir(sprintf('%s\\*.h264',dataFolder));
    if isempty(files_list)
        tei.video_file = [];
    else
        tei.video_file = files_list(1).name;
    end
%     winopen(dataFolder);
    files_list = dir(sprintf('%s\\*.mp4',dataFolder));
    if isempty(files_list)
        tei.video_file_mp4 = [];
    else
        tei.video_file_mp4 = files_list(1).name;
    end
    ei{ii} = tei;
end