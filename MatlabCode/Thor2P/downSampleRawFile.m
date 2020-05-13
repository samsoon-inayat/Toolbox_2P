function downSampleRawFile(ei)

% clear all
% clc
% dataFolder = 'F:\Sam\Data\Animal_160513\CoverSlip940\05_29_2016\freeRun';
% ei = thorGetExperimentInfo(dataFolder);
if exist(ei.dRawFile)
    display('Down sampled raw file already exists');
    return;
else
    delete(ei.dRawFile);
end

oddFrames = 1:2:ei.totalFrames;
evenFrames = 2:2:ei.totalFrames;
rows = ei.pixelY;
cols = ei.pixelX;

hWaitBar = waitbar(0,sprintf('Processed Frames -'));
fid = fopen(ei.mcRawFile,'r');
fidd = fopen(ei.dRawFile,'w');

for jj = 1:length(oddFrames)
    waitbar(jj/length(oddFrames),hWaitBar,sprintf('Processing Frame %d',jj));
    fseek(fid, (oddFrames(jj)-1)*rows*cols*2, 'bof');
    I1 = fread(fid,rows*cols,'uint16=>double',0,'l');
    I1 = reshape(I1,rows,cols)';
    
    fseek(fid, (evenFrames(jj)-1)*rows*cols*2, 'bof');
    I2 = fread(fid,rows*cols,'uint16=>double',0,'l');
    I2 = reshape(I2,rows,cols)';
    
    mI = (I1 + I2)/2;
    mI = mI';
    mI = uint16(reshape(mI,numel(mI),1));
    fwrite(fidd,mI,'uint16',0,'l');
end

fclose(fid);
fclose(fidd);
close(hWaitBar) ;