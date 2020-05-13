function raw2tif (ei,varargin)

% %%
% clear all
% clc

% dataFolder = 'F:\Sam\Data\Animal_153465\05_17_2016\freeRun';
% dataFolder = 'F:\Sam\Data\Animal_160513\CoverSlip940\05_29_2016\freeRun';
% ei = thorGetExperimentInfo(dataFolder);

row = ei.pixelY;
col = ei.pixelX;
if nargin < 2
    frameNumbers = 1:ei.totalFrames;
else
    ssF = varargin{1}(1);
    seF = varargin{1}(2);
    frameNumbers = ssF:seF;
end
rawFileN = uigetfile(sprintf('%s//*.raw',ei.dataFolder))
rawFile = makeName(rawFileN,ei.dataFolder);
tifFolder = uigetdir(ei.dataFolder);
% rawFile = makeName('Image_0001_0001_MC.raw',dataFolder);
% rawFileMC = sprintf('%s_MC.raw',rawFile(1:(length(rawFile)-4)));
hWaitBar = waitbar(0,sprintf('Processed Frames -'));
fid = fopen(rawFile,'r');
for ii = 1:length(frameNumbers)
    frameNumber = frameNumbers(ii);
    fseek(fid, (frameNumber-1)*row*col*2, 'bof');
    I1 = fread(fid,row*col,'uint16=>double',0,'l'); 
    Z = reshape(I1,row,col); Z1 = Z';
    fileName = makeName(sprintf('frame_%0.7d.tif',frameNumber),tifFolder);
    imwrite(uint16(Z1),fileName,'tif');
    waitbar(ii/length(frameNumbers),hWaitBar,sprintf('Processing Frame %d',frameNumber));
end
fclose(fid);
close(hWaitBar) 



