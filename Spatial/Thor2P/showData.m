clear all
clc

dataFolder = 'F:\Sam\Data\Animal_153465\05_17_2016\freeRun';
expInfo = thorGetExperimentInfo(dataFolder);
abfFiles = dir(sprintf('%s\\*.abf',dataFolder));
[abfd,abfsi,abfh]=abf2load(makeName(abfFiles(1).name,dataFolder));
frameTrigger = abfd(:,1);
encoderA = abfd(:,1);
encoderB = abfd(:,2);
encoderR = abfd(:,3);
reward = abfd(:,4);

timeAxis = linspace(0,abfsi*length(reward)*1e-6,length(reward));


row = expInfo.pixelY;
col = expInfo.pixelX;
frameNumbers = 1:3:expInfo.totalFrames;
rawFile = makeName('Image_0001_0001.raw',dataFolder);
rawFile1 = makeName('Image_0001_0001_1.raw',dataFolder);
hWaitBar = waitbar(0,sprintf('Processed Frames -'));
fid = fopen(rawFile,'r');
fid1 = fopen(rawFile1,'w');
for ii = 1:length(frameNumbers)
    frameNumber = frameNumbers(ii);
    fseek(fid, (frameNumber-1)*row*col*2, 'bof');
    I1 = fread(fid,row*col,'uint16=>double',0,'l'); 
    
    frameNumber = frameNumbers(ii)+1;
    fseek(fid, (frameNumber-1)*row*col*2, 'bof');
    I2 = fread(fid,row*col,'uint16=>double',0,'l'); 
    
    frameNumber = frameNumbers(ii)+2;
    fseek(fid, (frameNumber-1)*row*col*2, 'bof');
    I3 = fread(fid,row*col,'uint16=>double',0,'l'); 
    
    Im = (I1+I2+I3)/3;
    
    fwrite(fid1,uint16(Im),'uint16',0,'l');
    waitbar(ii/length(frameNumbers),hWaitBar,sprintf('Processing Frame %d',frameNumber));
    
%     Z = reshape(Im,row,col);
%     Z = Z';
%     figure(1);clf;
%     k=imagesc(Z);
%     colormap gray;
%     titleText = sprintf('Frame Number = %d',frameNumber); title(titleText);
%     pause(0.01);
end
fclose(fid);
fclose(fid1);
close(hWaitBar) 