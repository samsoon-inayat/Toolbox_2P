function all_tif_file_names = get_tif_file_names (tif_folder_name,wildcard)

wild_card = sprintf('%s\\%s*.tif',tif_folder_name,wildcard);
fileNames = dir(wild_card);
for ii = 1:length(fileNames)
    all_tif_file_names{ii} = fullfile(tif_folder_name,fileNames(ii).name);
end