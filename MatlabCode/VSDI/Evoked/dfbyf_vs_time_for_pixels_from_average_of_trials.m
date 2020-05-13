function [dfbyf,xics] = dfbyf_vs_time_for_pixels_from_average_of_trials (ags,pags,gg,aa,stimType,frames_info,maskFactor,pixels)
animal = ags(gg).animals{aa};
panimal = pags(gg).animals{aa};
mask = getMask(ags(gg).animals{aa}.root_folder,maskFactor);
masko = getMask(ags(gg).animals{aa}.root_folder,1);
% for aa = aa%1:length(animal.data_folders)

dataFolder = animal.root_folder
mask = getMask(dataFolder);
if isempty(mask)
    display('Mask does not exist');
    dfbyf = [];
    return;
end
peDataFolder = panimal.root_folder;
eRecordings = panimal.eRecordings;
for ii = 1:length(eRecordings)
    if strcmp(stimType,eRecordings{ii}.stimulus.stimulus_type)
        break;
    end
end
evFolder = eRecordings{ii}.root_folder;
try
    fileName = makeName('dfbyfo.mat',evFolder);
    temp = load(fileName);
catch
    fileName = makeName('mean_hi_lo.mat',evFolder);
    temp = load(fileName);
end
if ~isfield(temp,'dfbyfo')
    temp.dfbyfo = temp.df;
end
dfbf1 = temp.dfbyfo;

dfbf = zeros(size(mask,1),size(mask,2),size(dfbf1,2));
dfbf = reshape(dfbf,size(mask,1)*size(mask,2),size(dfbf1,2));
maskI = find(masko);
dfbf(maskI,:) = dfbf1;
dfbf = reshape(dfbf,size(mask,1),size(mask,2),size(dfbf1,2));
df = applyMask(dfbf,mask);

df = applySpatialFilter(df);
%     df = applyTemporalFilter(df);
stimulusFrame = find((frames_info(4):frames_info(5)) == frames_info(1));
srF = stimulusFrame;
erF = srF + 14;
ebF = stimulusFrame - 1;
sbF = ebF - 14;
for jj = 1:length(pixels)
    thisPixels = pixels{jj};
    thisMask = makeMask(size(df(:,:,1)),thisPixels);
    thisdfbyf = mean(getMaskedValues(df,thisMask));
    dfbyf(jj,:) = thisdfbyf';
end
xics = getXicsFromResponse(dfbyf,frames_info);
xics.avg_response = mean(df(:,:,srF:erF),3)-mean(df(:,:,sbF:ebF),3);
% end