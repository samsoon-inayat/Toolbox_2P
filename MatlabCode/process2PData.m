clear all
clc
% folderName = 'F:\Sam\Data\CalciumImaging\Test\rsc_prism_20160417_001';
folderName = 'F:\Sam\Data\AnimalID_xxx\04_27_2016\R1\TSeries\HPC_axons_RSC002\images';
pFolderName = 'F:\Sam\ProcessedData\AnimalID_xxx\04_27_2016\R1\TSeries\HPC_axons_RSC002\images';
fileNameTag = 'Image';
wildCard = sprintf('%s*.tif',makeName(fileNameTag,folderName));
fileNames = dir(wildCard);
for ii = 1:length(fileNames)
    fns{ii} = makeName(fileNames(ii).name,folderName);
end

avgImage = double(imread(fns{1},'tif'));
for ii = 2:length(fileNames)
    ii
    fns{ii} = makeName(fileNames(ii).name,folderName);
    avgImage = avgImage + double(imread(fns{ii},'tif'));
end
avgImage = avgImage/length(fileNames);
figure(1);clf;imagesc(avgImage);
colormap gray
axis equal
axis off
averageImage = avgImage;
fn = makeName('averageImage.mat',pFolderName);
save(fn,'averageImage');

for ii = 1:length(fileNames)
    ii
    fns{ii} = makeName(fileNames(ii).name,folderName);
    imageStack(:,:,ii) = double(imread(fns{ii},'tif'));
end

roifns = {'0086-0598.roi'
    '0087-0475.roi'
    '0137-0501.roi'
    '0187-0391.roi'
    '0195-0362.roi'
    '0237-0214.roi'
    '0258-0236.roi'
    '0290-0269.roi'
    '0320-0355.roi'
    '0340-0215.roi'
    '0353-0109.roi'
    '0387-0256.roi'
    '0392-0330.roi'
    '0402-0358.roi'
    '0483-0347.roi'};

roiFolder = 'F:\Sam\Data\AnimalID_xxx\04_27_2016\R1\TSeries\HPC_axons_RSC002\images\forMajid\rois';

for ii = 1:15
    roicstr{ii} = makeName(roifns{ii},roiFolder);
end

[cvsROIs] = ReadImageJROI(roicstr) ;
pFolderName = 'F:\Sam\ProcessedData\AnimalID_xxx\04_27_2016\R1\TSeries\HPC_axons_RSC002\images';
fn = makeName('averageImage.mat',pFolderName);
load(fn);

figure(10);clf;
imagesc(averageImage);
colormap gray;
axis equal
axis off

for ii = 1:15
    rect = cvsROIs{ii}.vnRectBounds;
    rectangle('Position',rect,'Curvature',0.2)
end


% plot of time series
xlFile = 'F:\Sam\Data\AnimalID_xxx\04_27_2016\R1\TSeries\HPC_axons_RSC002\images\forMajid\Results1.xlsx';
calSig = xlsread(xlFile,'MeansSheet','A2:O5001');

allMean = mean(averageImage(:));
sCalSig = calSig - allMean;


frameRate = 24.54;
totalTime = 5000/24.54;
timeAxis = linspace(0,totalTime,5000);

interesting = 1:15;
figure(1);clf;
for ii = 1:length(interesting)
    subplot(length(interesting),1,ii);
    plot(timeAxis,sCalSig(:,interesting(ii)));
    set(gca,'FontSize',12);
    ylabel('A.U.');
    if ii < length(interesting)
        set(gca,'XTickLabel','');
    end
    xlim([0 totalTime]);
    grid on
end
xlabel('secs');

