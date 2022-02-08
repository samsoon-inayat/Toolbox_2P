function T = check_for_video_files(T)

sz = size(T,2)+1;
T = [T table(cell(size(T,1),1))]
for ii = 1:size(T,1)
    dataFolder = cell2mat(T{ii,6});
    files_list = dir(sprintf('%s\\*.h264',dataFolder));
    if isempty(files_list)
        continue;
    end
    T(ii,sz) = {files_list(1).name};
%     winopen(dataFolder);
end