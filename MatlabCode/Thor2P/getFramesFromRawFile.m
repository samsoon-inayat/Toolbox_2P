function frames = getFramesFromRawFile (rawFile,pixelX,pixelY,totalFrames)

row = pixelY;
col = pixelX;
frameNumbers = 1:totalFrames;
hWaitBar = waitbar(0,sprintf('Getting Frames -'));
fid = fopen(rawFile,'r');
for ii = 1:length(frameNumbers)
    frameNumber = frameNumbers(ii);
    fseek(fid, (frameNumber-1)*row*col*2, 'bof');
    I1 = fread(fid,row*col,'uint16=>double',0,'l'); 
    Z = reshape(I1,row,col); 
    frames(:,:,ii) = Z';
%     imwrite(uint16(Z1),fileName,'tif');
    waitbar(ii/length(frameNumbers),hWaitBar,sprintf('Getting Frame %d//%d',ii,length(frameNumbers)));
end
fclose(fid);
close(hWaitBar) 



