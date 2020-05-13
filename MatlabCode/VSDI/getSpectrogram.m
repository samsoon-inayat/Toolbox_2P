function sp = getSpectrogram(allAnimals,animalType,nums,stimType,roi)
% %% Spectofram for Evoked Data

% animalType = 5;
% nums = 1;
% stimType = 'FL';

coord = get_ROI_coordinates(allAnimals,animalType,nums);
pixels = get_ROI_pixels(coord,[3 3]);
names = fieldnames(coord);names(1) = [];
regionNum = roi;

animal = allAnimals(animalType);
% for an = nums%1:length(animal.data_folders)
an = nums;
dataFolder = animal.data_folders{an};
mask = getMask(dataFolder);

if isempty(mask)
    display('Mask does not exist');
    dfbyf = [];
    return;
end
peDataFolder = animal.processed_data_folders{an}
fileName = makeName('folders_data.mat',peDataFolder);
temp = load(fileName);
folders_data = temp.folders_data; clear temp;
for ii = 1:length(folders_data)
    if strcmp(stimType,folders_data(ii).stimulus.stimulus_type)
        break;
    end
end
evFolder = makeName(folders_data(ii).name,peDataFolder);

fileName = sprintf('spectrogram_%s_%s.mat',stimType,names{roi});
fileName = makeName(fileName,evFolder);
sp = load(fileName);
