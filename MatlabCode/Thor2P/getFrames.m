function frames = getFrames (savePath,frameNumbers)
tifFolder = savePath;
hWaitBar = waitbar(0,sprintf('Reading Frames -'));
for ii = 1:length(frameNumbers)
    frameNumber = frameNumbers(ii);
    fileName = makeName(sprintf('frame_%0.7d.tif',frameNumber),tifFolder);
    frames(:,:,ii) = double(imread(fileName));
    waitbar(ii/length(frameNumbers),hWaitBar,sprintf('Reading Frame %d',frameNumber));
end
close(hWaitBar) 



