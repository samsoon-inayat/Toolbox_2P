function visualizeWithImageJ (animal,nums,stimType,trial)

for an = nums%1:length(animal.data_folders)
    dataFolder = animal.data_folders{an};
    mask = getMask(dataFolder);
    if isempty(mask)
        display('Mask does not exist');
        break;
    end
    peDataFolder = animal.processed_data_folders{an}
    fileName = makeName('folders_data.mat',peDataFolder);
    temp = load(fileName);
    folders_data = temp.folders_data; clear temp;
    for ii = 1:length(folders_data)
        if ~strcmp(folders_data(ii).stimulus.stimulus_type,stimType)
            continue;
        end
        evFolder = makeName(folders_data(ii).name,peDataFolder);
        fileName = makeName(sprintf('hi%d.mat',trial),evFolder);
        temp = load(fileName);
        temp.df = applySpatialFilter(temp.dfbyfo);
%         temp.dfbyfo = applyMask(temp.dfbyfo,mask,'same');
        fileName = mat2tif(temp.df,'tempDeleteTheFile',pwd);
        winopen(fileName);
        n = 0;
    end
end
