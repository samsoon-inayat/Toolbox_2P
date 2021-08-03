function frames = getFramesFromRaw (ei,varargin)

row = ei.pixelY;
col = ei.pixelX;
frameNumbers = 1:ei.totalFrames;
if nargin > 1
    frameNumbers = varargin{1};
end
% rawFileN = uigetfile(sprintf('%s//*.raw',ei.dataFolder))
% rawFile = makeName(rawFileN,ei.dataFolder);
rawFile = ei.rawFile;
% tifFolder = uigetdir(ei.dataFolder);

% rawFile = makeName('Image_0001_0001_MC.raw',dataFolder);
% rawFileMC = sprintf('%s_MC.raw',rawFile(1:(length(rawFile)-4)));
% figure(1000);clf;
hWaitBar = waitbar(0,sprintf('Getting Frames -'));
fid = fopen(rawFile,'r');
if fid<0
    close(hWaitBar);
    frames = getFramesFromTifs(ei,frameNumbers);
    return;
end
for ii = 1:length(frameNumbers)
    frameNumber = frameNumbers(ii);
    fseek(fid, (frameNumber-1)*row*col*2, 'bof');
    I1 = fread(fid,row*col,'uint16=>double',0,'l'); 
    Z = reshape(I1,row,col); frames(:,:,ii) = Z';
%     imwrite(uint16(Z1),fileName,'tif');
    waitbar(ii/length(frameNumbers),hWaitBar,sprintf('Getting Frame %d//%d',ii,length(frameNumbers)));
end
fclose(fid);
close(hWaitBar) 



