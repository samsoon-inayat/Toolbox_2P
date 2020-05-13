%%clear all
clc
% folderName = 'F:\Sam\Data\CalciumImaging\Test\rsc_prism_20160417_001';
% folderName = 'F:\Sam\Data\AnimalID_xxx\04_27_2016\R1\TSeries\HPC_axons_RSC002\images';
% pFolderName = 'F:\Sam\ProcessedData\AnimalID_xxx\04_27_2016\R1\TSeries\HPC_axons_RSC002\images';

folderName = 'F:\Sam\Data\AnimalID_xxx\04_27_2016\R1\TSeries\HPC_axons_RSC\images';
pFolderName = 'F:\Sam\Data\AnimalID_xxx\04_27_2016\R1\TSeries\HPC_axons_RSC\images';

fileNameTag = 'Image';
wildCard = sprintf('%s*.tif',makeName(fileNameTag,folderName));
fileNames = dir(wildCard);
for ii = 1:length(fileNames)
    fns{ii} = makeName(fileNames(ii).name,folderName);
end
% seeData(fns); 

% avgImage = double(imread(fns{1},'tif'));
for ii = 1:length(fileNames)
    ii
    fns{ii} = makeName(fileNames(ii).name,folderName);
    imageStack(:,:,ii) = double(imread(fns{ii},'tif'));
%     avgImage = avgImage + double(imread(fns{ii},'tif'));
end
%%
% avgImage = avgImage/length(fileNames);
averageImage = mean(imageStack,3);
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
xlFile = 'F:\Sam\Data\AnimalID_xxx\04_27_2016\R1\TSeries\HPC_axons_RSC\ROIs\calSig.xlsx';
calSig = xlsread(xlFile,'Sheet1','A2:U1667');

allMean = mean(imageStack(:));
sCalSig = calSig - allMean;


frameRate = 24.69/3;
totalTime = length(fileNames)/frameRate;
timeAxis = linspace(0,totalTime,length(fileNames));

%%
minCS = min(sCalSig(:))-50;
maxCS = max(sCalSig(:))+50;
interesting = [18 19 2 9 12 17 20];
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

