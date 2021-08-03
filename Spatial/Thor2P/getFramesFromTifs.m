function frames = getFramesFromTifs (ei,varargin)

rawDataFolder = ei.folders.rawDataFolder;

frameNumbers = 1:ei.totalFrames;
if nargin > 1
    frameNumbers = varargin{1};
end

hWaitBar = waitbar(0,sprintf('Getting Frames -'));
for ii = 1:length(frameNumbers)
    frameNumber = frameNumbers(ii);
    if frameNumber < 10000
        fileName = sprintf('ChanA_0001_0001_0001_%.4d.tif',frameNumber);
    else
        fileName = sprintf('ChanA_0001_0001_0001_%.5d.tif',frameNumber);
    end
    fileName = makeName(fileName,rawDataFolder);
    Z = double(imread(fileName)); 
    frames(:,:,ii) = Z;
    waitbar(ii/length(frameNumbers),hWaitBar,sprintf('Getting Frame %d//%d',frameNumber,length(frameNumbers)));
end
close(hWaitBar) 



