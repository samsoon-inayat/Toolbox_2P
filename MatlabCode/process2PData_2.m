%%
clear all
clc
% folderName = 'F:\Sam\Data\CalciumImaging\Test\rsc_prism_20160417_001';
% folderName = 'F:\Sam\Data\AnimalID_xxx\04_27_2016\R1\TSeries\HPC_axons_RSC002\images';
% pFolderName = 'F:\Sam\ProcessedData\AnimalID_xxx\04_27_2016\R1\TSeries\HPC_axons_RSC002\images';

folderName = 'F:\Sam\Data\AnimalID_xxx\04_27_2016\R1\TSeries\HPC_axons_RSC\images1';
pFolderName = 'F:\Sam\Data\AnimalID_xxx\04_27_2016\R1\TSeries\HPC_axons_RSC\images1';

fileNameTag = 'Reduced HPC_axon';
wildCard = sprintf('%s*.tif',makeName(fileNameTag,folderName));
fileNames = dir(wildCard);
for ii = 1:length(fileNames)
    fns{ii} = makeName(fileNames(ii).name,folderName);
end
% seeData(fns); 

avgImage = double(imread(fns{1},'tif'));
for ii = 2:length(fileNames)
    ii
    fns{ii} = makeName(fileNames(ii).name,folderName);
%     imageStack(:,:,ii) = double(imread(fns{ii},'tif'));
    avgImage = avgImage + double(imread(fns{ii},'tif'));
end
%%
avgImage = avgImage/length(fileNames);
% averageImage = mean(imageStack,3);
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
xlFile = 'F:\Sam\Data\AnimalID_xxx\04_27_2016\R1\TSeries\HPC_axons_RSC\ROIs1\calSig.xlsx';
calSig = xlsread(xlFile,'Sheet1','A2:N5001');

allMean = mean(averageImage(:));
sCalSig = calSig - allMean;


frameRate = 24.69;
totalTime = length(fileNames)/frameRate;
timeAxis = linspace(0,totalTime,length(fileNames));

%%
minCS = min(sCalSig(:))-50;
maxCS = max(sCalSig(:))+50;
interesting = [1 2 8 13 7 11];
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

