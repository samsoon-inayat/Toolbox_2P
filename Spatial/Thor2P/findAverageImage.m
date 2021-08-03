function findAverageImage (ei)

% %%
% clear all
% clc

% dataFolder = 'F:\Sam\Data\Animal_153465\05_17_2016\freeRun';
% dataFolder = 'F:\Sam\Data\Animal_160513\CoverSlip940\05_29_2016\freeRun';
% ei = thorGetExperimentInfo(dataFolder);
if exist(ei.averageImageFile)
    display('Average Image File already exists');
    return;
end

row = ei.pixelY;
col = ei.pixelX;
frameNumbers = 1:ei.totalFrames;
rawFileN = uigetfile(sprintf('%s//*.raw',ei.dataFolder))
rawFile = makeName(rawFileN,ei.dataFolder);
% rawFile = makeName('Image_0001_0001_MC.raw',dataFolder);
% rawFileMC = sprintf('%s_MC.raw',rawFile(1:(length(rawFile)-4)));
hWaitBar = waitbar(0,sprintf('Processed Frames -'));
fid = fopen(rawFile,'r');
sumI = zeros(row*col,1);
maxI = zeros(row*col,1);
% allImages = [];
for ii = 1:length(frameNumbers)
    frameNumber = frameNumbers(ii);
    fseek(fid, (frameNumber-1)*row*col*2, 'bof');
    I1 = fread(fid,row*col,'uint16=>double',0,'l'); 
    waitbar(ii/length(frameNumbers),hWaitBar,sprintf('Processing Frame %d',frameNumber));
    sumI = sumI + I1;
    maxI = max(maxI,I1);
%     allImages(:,ii) = I1;
end
fclose(fid);
close(hWaitBar) 
sumI = sumI/length(frameNumbers);
Z = reshape(sumI,row,col); Z = Z';
figure(1);clf;
k=imagesc(Z);
colormap gray;
titleText = sprintf('Frame Number = %d',frameNumber); title(titleText);
axis equal
axis off;
averageImage = Z;
m = matfile(ei.averageImageFile,'Writable',true);
tmpIdx = strfind(rawFileN,'.');
theName = rawFileN(1:(tmpIdx-1));
cmdText = sprintf('m.%s_averageImage = averageImage;',theName);
eval(cmdText);
imwrite(uint16(averageImage),ei.averageImageTifFile,'tif');

maxIm = reshape(maxI,row,col);
maxIm =maxIm';
cmdText = sprintf('m.%s_maxImage = maxIm;',theName);
eval(cmdText);
imwrite(uint16(maxIm),ei.maxImageTifFile,'tif');


