%%
clear all
clc
% folderName = 'F:\Sam\Data\CalciumImaging\Test\rsc_prism_20160417_001';
% folderName = 'F:\Sam\Data\AnimalID_xxx\04_27_2016\R1\TSeries\HPC_axons_RSC002\images';
% pFolderName = 'F:\Sam\ProcessedData\AnimalID_xxx\04_27_2016\R1\TSeries\HPC_axons_RSC002\images';

folderName = 'F:\Sam\Data\AnimalID_xxx\04_27_2016\R1\TSeries\HPC_axons_RSC';
pFolderName = 'F:\Sam\Data\AnimalID_xxx\04_27_2016\R1\TSeries\HPC_axons_RSC';

fileNameTag = 'Projection';
wildCard = sprintf('%s*.tif',makeName(fileNameTag,folderName));
fileNames = dir(wildCard);
for ii = 1:length(fileNames)
    fns{ii} = makeName(fileNames(ii).name,folderName);
end
% seeData(fns); 

% avgImage = double(imread(fns{1},'tif'));
% for ii = 2:length(fileNames)
%     ii
%     fns{ii} = makeName(fileNames(ii).name,folderName);
% %     imageStack(:,:,ii) = double(imread(fns{ii},'tif'));
%     avgImage = avgImage + double(imread(fns{ii},'tif'));
% end
% %%
% avgImage = avgImage/length(fileNames);
% averageImage = mean(imageStack,3);
avgImage = imread(makeName('Average_Intensity_Projection of HPC_axons_RSC.tif',folderName),'tif');
averageImage = avgImage;
figure(1);clf;imagesc(averageImage);
colormap gray
axis equal
axis off
fn = makeName('averageImage.mat',pFolderName);
save(fn,'averageImage');


figure(10);clf;
imagesc(averageImage);
colormap gray;
axis equal
axis off
%%
% plot of time series
xlFile = makeName('Results.xlsx',folderName);
[num,txt,raw] = xlsread(xlFile);
indicesOfMeanCalSig = [];
for ii = 1:size(txt,2)
    if strfind(txt{ii},'Mean')
        indicesOfMeanCalSig = [indicesOfMeanCalSig ii];
    end
end
calSig = num(:,indicesOfMeanCalSig);

allMean = mean(averageImage(:));
sCalSig = calSig - allMean;


frameRate = 24.69;
totalTime = length(fileNames)/frameRate;
timeAxis = linspace(0,totalTime,length(fileNames));

%%
minCS = min(sCalSig(:))-50;
maxCS = max(sCalSig(:))+50;
% interesting = 1:size(calSig,2);
interesting = [1 2 5 6 7 14];
figure(1);clf;
for ii = 1:length(interesting)
    subplot(length(interesting),1,ii);
    plot(timeAxis,sCalSig(:,interesting(ii)));
    set(gca,'FontSize',12);
    ylabel(sprintf('%d',interesting(ii)));
    if ii < length(interesting)
        set(gca,'XTickLabel','','XTick',[]);
    end
    set(gca,'YTick',[0 round(maxCS)])
    xlim([0 totalTime]);
    ylim([minCS maxCS]);
    box off
%     grid on
end
xlabel('secs');

