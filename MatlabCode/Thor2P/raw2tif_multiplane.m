function raw2tif_multiplane (fileName,savePath,ei,varargin)
p = inputParser;
addRequired(p,'fileName',@ischar);addRequired(p,'savePath',@ischar);addRequired(p,'ei',@isstruct);
addOptional(p,'bidishift',0,@isnumeric);
addOptional(p,'plane',1,@isnumeric);
addOptional(p,'number_of_planes',1,@isnumeric);
parse(p,fileName,savePath,ei,varargin{:});

row = ei.pixelY;
col = ei.pixelX;
raw_averageImage = zeros(row,col);
raw_maxImage = zeros(row,col);
tempImage = zeros(row,col);
plane = p.Results.plane;
numberOfPlanes = p.Results.number_of_planes;

frameNumbers = plane:numberOfPlanes:(ei.streaming_frames-ei.totalFrames);

bidishift = p.Results.bidishift;
rawFile = fileName;
tifFolder = savePath;
rawDataFolder = ei.folders.rawDataFolder;
hWaitBar = waitbar(0,sprintf('Processed Frames -'));
fid = fopen(rawFile,'r');


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
    waitbar(ii/length(frameNumbers),hWaitBar,sprintf('Processing Frame %d/%d',frameNumber,length(frameNumbers)));
end
raw_averageImage = raw_averageImage./length(frameNumbers);

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



