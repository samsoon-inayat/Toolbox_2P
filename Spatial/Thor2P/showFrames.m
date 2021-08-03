function showFrames (fileName,savePath,ei,varargin)

row = ei.pixelY;
col = ei.pixelX;
frameNumbers = 1:ei.totalFrames;
if nargin > 3
    ssF = varargin{1}(1);
    seF = varargin{1}(2);
    frameNumbers = ssF:seF;
end
% rawFileN = uigetfile(sprintf('%s//*.raw',ei.dataFolder))
% rawFile = makeName(rawFileN,ei.dataFolder);
rawFile = fileName;
% tifFolder = uigetdir(ei.dataFolder);
tifFolder = savePath;
% rawFile = makeName('Image_0001_0001_MC.raw',dataFolder);
% rawFileMC = sprintf('%s_MC.raw',rawFile(1:(length(rawFile)-4)));
figure(1000);clf;
% hWaitBar = waitbar(0,sprintf('Processed Frames -'));
fid = fopen(rawFile,'r');
for ii = 1:length(frameNumbers)
    frameNumber = frameNumbers(ii);
    fseek(fid, (frameNumber-1)*row*col*2, 'bof');
    I1 = fread(fid,row*col,'uint16=>double',0,'l'); 
    Z = reshape(I1,row,col); Z1 = Z';
    fileName = makeName(sprintf('frame_%0.7d.tif',frameNumber),tifFolder);
    imagesc(Z1);axis equal;
    title(frameNumber);
    pause(0.01);
%     imwrite(uint16(Z1),fileName,'tif');
%     waitbar(ii/length(frameNumbers),hWaitBar,sprintf('Processing Frame %d',frameNumber));
end
fclose(fid);
% close(hWaitBar) 



