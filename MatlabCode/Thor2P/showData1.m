clear all
clc

plotsEncoder = 1;

% dataFolder = 'F:\Sam\Data\Animal_153465\05_17_2016\freeRun';
dataFolder = 'F:\Sam\Data\Animal_160513\CoverSlip940\05_29_2016\freeRun';
expInfo = thorGetExperimentInfo(dataFolder);
abfFiles = dir(sprintf('%s\\*.abf',dataFolder));
[abfd,abfsi,abfh]=abf2load(makeName(abfFiles(1).name,dataFolder));
frameTrigger = abfd(:,1);
encoderA = abfd(:,2);
encoderB = abfd(:,3);
encoderR = abfd(:,4);
reward = abfd(:,5);
timeAxis = linspace(0,size(abfd,1)*abfsi*1e-6,size(abfd,1));
timeAxis = timeAxis';
abfd(:,6) = timeAxis;
% figure(1);clf;
% plot(timeAxis,encoderA);
% hold on;
% plot(timeAxis,encoderB,'r');

% makeFigure(abfd);


 
timeAxis = linspace(0,abfsi*length(reward)*1e-6,length(reward));

try
    row = expInfo.pixelY;
    col = expInfo.pixelX;
    startN = 1;
    endN = 1000;
    frameNumbers = startN:endN;
    rawFile = makeName('Image_0001_0001_MC.raw',dataFolder);
    fid = fopen(rawFile,'r');
    for ii = 1:length(frameNumbers)
        frameNumber = frameNumbers(ii);
        fseek(fid, (frameNumber-1)*row*col*2, 'bof');
        I1 = fread(fid,row*col,'uint16=>double',0,'l');     
        Z = reshape(I1,row,col)';
        figure(1);clf;
        k=imagesc(Z);
        colormap gray;
        titleText = sprintf('Frame Number = %d',frameNumber); title(titleText);
        axis equal
        axis off;
    %     pause(0.01);
    end
    fclose(fid);
catch
    msgstr = lasterr
    error(msgstr);
    fclose(fid);
end