function frames = getFramesFromRaw_simple (rawFile,rows,cols,frameNumbers)

hWaitBar = waitbar(0,sprintf('Getting Frames -'));
fid = fopen(rawFile,'r');
for ii = 1:length(frameNumbers)
    frameNumber = frameNumbers(ii);
    fseek(fid, (frameNumber-1)*rows*cols*2, 'bof');
    I1 = fread(fid,rows*cols,'uint16=>double',0,'l'); 
    Z = reshape(I1,rows,cols); frames(:,:,ii) = Z';
%     imwrite(uint16(Z1),fileName,'tif');
    waitbar(ii/length(frameNumbers),hWaitBar,sprintf('Getting Frame %d//%d',frameNumber,length(frameNumbers)));
end
fclose(fid);
close(hWaitBar) 



