function motionCorrection (ei)

% 
% dataFolder = 'F:\Sam\Data\Animal_160513\CoverSlip940\05_29_2016\freeRun';
% ei = thorGetExperimentInfo(dataFolder);
% abfFiles = dir(sprintf('%s\\*.abf',dataFolder));
% [abfd,abfsi,abfh]=abf2load(makeName(abfFiles(1).name,dataFolder));
% frameTrigger = abfd(:,1);
% encoderA = abfd(:,1);
% encoderB = abfd(:,2);
% encoderR = abfd(:,3);
% reward = abfd(:,4);
% 
% timeAxis = linspace(0,abfsi*length(reward)*1e-6,length(reward));

if ~exist(ei.registerationImageFile)
    display('Registration Image File missing');
    return;
end

if exist(ei.mcRawFile)
    display('Motion corrected raw file already exists');
    return;
end


%% parallel for loop in the nest to speed up processing
framesPerSet = 500;
setsOfFrames = floor(ei.totalFrames/framesPerSet); % find how many sets of frames to be processed  with par for
startOfFramesSet = 1:framesPerSet:ei.totalFrames;
endOfFramesSet = (1:framesPerSet:ei.totalFrames)-1;
endOfFramesSet(1) = []; endOfFramesSet(end+1) = ei.totalFrames;

m = matfile(ei.registerationImageFile);
templateImage = m.registerationImage;
row = ei.pixelY;
col = ei.pixelX;
try
    if exist(ei.mcRawFile)
        delete(ei.mcRawFile);
    end
    fidw = fopen(ei.mcRawFile,'w');
    for jj = 1:setsOfFrames
        frameNumbers = startOfFramesSet(jj):endOfFramesSet(jj);
        hWaitBar = waitbar(0,sprintf('Processed Frames -'));
        fid = fopen(ei.rawFile,'r');
        cImages = [];
        oImages = [];
        for ii = 1:length(frameNumbers)
            frameNumber = frameNumbers(ii);
            fseek(fid, (frameNumber-1)*row*col*2, 'bof');
            I1 = fread(fid,row*col,'uint16=>double',0,'l'); 
            I1 = reshape(I1,row,col);
            oImages(:,:,ii) = I1';
            waitbar(ii/length(frameNumbers),hWaitBar,sprintf('Processing Frame %d',frameNumber));
        end
        fclose(fid);
        close(hWaitBar) ;
        display('Par for start');
        parfor ii = 1:length(frameNumbers)
            frameNumber = frameNumbers(ii);
            I1 = oImages(:,:,ii);
             [output Greg] = dftregistration(fft2(templateImage),fft2(I1),1);   % translational registration  
            I1c = abs(ifft2(Greg));
           cImages(:,:,ii) = I1c;
        end
        display('Par for end');
        hWaitBar = waitbar(0,sprintf('Processed Frames -'));
        for ii = 1:length(frameNumbers)
            frameNumber = frameNumbers(ii);
            I1 = cImages(:,:,ii)';
            I1 = reshape(I1,row*col,1);
            I1 = fwrite(fidw,uint16(I1),'uint16',0,'l'); 
            waitbar(ii/length(frameNumbers),hWaitBar,sprintf('Processing Frame %d',frameNumber));
        end
        close(hWaitBar) ;
        display(sprintf('Done with Set # %d of %d',jj,setsOfFrames));
    end
    fclose(fidw);
catch
    msgstr = lasterr
    error(msgstr);
    fclose(fidw);
end
%% 
% for ii = 1:500
%     figure(1);clf;
%     imagesc(cImages(:,:,ii));
%     axis off; axis equal;
% end

%%
% acImage = mean(cImages,3);
% figure(1);clf;
% subplot 121;
% imagesc(acImage);
% axis off; axis equal;
% 
% aoImage = mean(oImages,3);
% subplot 122
% imagesc(aoImage);
% axis off; axis equal;
% colormap gray