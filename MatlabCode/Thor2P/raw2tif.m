function raw2tif (fileName,savePath,ei,varargin)
p = inputParser;
default_frameNumbers = 1:ei.totalFrames;
default_downsamplespace = 1;
default_downsampletime = 1;
default_plane = 1;
default_bidishift = 0;
addRequired(p,'fileName',@ischar);addRequired(p,'savePath',@ischar);addRequired(p,'ei',@isstruct);
addOptional(p,'frameNumbers',default_frameNumbers,@isnumeric);
addOptional(p,'downsamplespace',default_downsamplespace,@isnumeric);
addOptional(p,'downsampletime',default_downsampletime,@isnumeric);
addOptional(p,'bidishift',default_bidishift,@isnumeric);
addOptional(p,'plane',default_plane,@isnumeric);
parse(p,fileName,savePath,ei,varargin{:});

row = ei.pixelY;
col = ei.pixelX;
raw_averageImage = zeros(row,col);
raw_maxImage = zeros(row,col);
tempImage = zeros(row,col);
frameNumbers = p.Results.frameNumbers;
downsamplespace = p.Results.downsamplespace;
downsampletime = p.Results.downsampletime;
plane = p.Results.plane;

if ei.zFastEnable
    frameNumbers = plane:(ei.zSteps+1):ei.streaming_frames;
end


if downsamplespace < 1
    raw_averageImage = imresize(raw_averageImage,downsamplespace);
    raw_maxImage = imresize(raw_maxImage,downsamplespace);
    tempImage = imresize(tempImage,downsamplespace);
end

bidishift = p.Results.bidishift;
rawFile = fileName;
tifFolder = savePath;
rawDataFolder = ei.folders.rawDataFolder;
hWaitBar = waitbar(0,sprintf('Processed Frames -'));
fid = fopen(rawFile,'r');

if downsampletime == 1
    for ii = 1:length(frameNumbers)
        frameNumber = frameNumbers(ii);
        if fid > 0
            fseek(fid, (frameNumber-1)*row*col*2, 'bof');
            I1 = fread(fid,row*col,'uint16=>double',0,'l'); 
            Z = reshape(I1,row,col); 
            Z1 = Z';
        else
            if frameNumber < 10000
                fileName = sprintf('ChanA_0001_0001_0001_%.4d.tif',frameNumber);
            else
                fileName = sprintf('ChanA_0001_0001_0001_%.5d.tif',frameNumber);
            end
            fileName = makeName(fileName,rawDataFolder);
            Z1 = double(imread(fileName));
        end

        if abs(bidishift) > 0 
            Z1 = ShiftBiDi(bidishift, Z1, row, col);
        end
        if downsamplespace < 1
            Z1 = imresize(Z1,downsamplespace);
        end
        raw_averageImage = raw_averageImage + Z1;
        tempImage(:,:,1) = raw_maxImage;
        tempImage(:,:,2) = Z1;
        raw_maxImage = max(tempImage,[],3);
        fileName = makeName(sprintf('frame_%0.7d.tif',ii),tifFolder);
        imwrite(uint16(Z1),fileName,'tif');
        waitbar(ii/length(frameNumbers),hWaitBar,sprintf('Processing Frame %d/%d',ii,length(frameNumbers)));
    end
    raw_averageImage = raw_averageImage./length(frameNumbers);
end

if downsampletime < 1
    jumpFactor = 1/downsampletime;
    ss = frameNumbers(1):jumpFactor:frameNumbers(end);
    se = (frameNumbers(1)+jumpFactor-1):jumpFactor:frameNumbers(end);
    [ml,idx] = min([length(ss) length(se)]);
    ss = ss(1:ml);
    se = se(1:ml);
    for iii = 1:length(ss)
        frameNumbers = ss(iii):se(iii);
        cumZ1 = zeros(row,col);
        for ii = 1:length(frameNumbers)
            frameNumber = frameNumbers(ii);
            
            if fid > 0
                fseek(fid, (frameNumber-1)*row*col*2, 'bof');
                I1 = fread(fid,row*col,'uint16=>double',0,'l'); 
                Z = reshape(I1,row,col); 
                Z1 = Z';
            else
                if frameNumber < 10000
                    fileName = sprintf('ChanA_0001_0001_0001_%.4d.tif',frameNumber);
                else
                    fileName = sprintf('ChanA_0001_0001_0001_%.5d.tif',frameNumber);
                end
                fileName = makeName(fileName,rawDataFolder);
                Z1 = double(imread(fileName));
            end
            if abs(bidishift) > 0 
                Z1 = ShiftBiDi(bidishift, Z1, row, col);
            end
            cumZ1 = cumZ1 + Z1;
        end
        Z1 = cumZ1/jumpFactor;
        if downsamplespace < 1
            Z1 = imresize(Z1,downsamplespace);
        end
        
        
        raw_averageImage = raw_averageImage + Z1;
        tempImage(:,:,1) = raw_maxImage;
        tempImage(:,:,2) = Z1;
        raw_maxImage = max(tempImage,[],3);
        
        fileName = makeName(sprintf('frame_%0.7d.tif',iii),tifFolder);
        imwrite(uint16(Z1),fileName,'tif');
        waitbar(iii/length(ss),hWaitBar,sprintf('Processing frame %d/%d',iii,length(ss)));
    end
    raw_averageImage = raw_averageImage./length(ss);
end

if fid > 0
    fclose(fid);    
end
close(hWaitBar);
try
pSaveDataFolder = ei.folders.thispFolder;
fileName = makeName('raw_averageImage.mat',pSaveDataFolder);
save(fileName,'raw_averageImage','raw_maxImage');
catch
    n=0;
end



