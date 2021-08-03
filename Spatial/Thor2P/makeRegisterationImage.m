function makeRegisterationImage (ei,varargin)
        
if nargin == 2
    fns(1) = varargin{1}(1);
    fns(2) = varargin{1}(2);
else
    if ~exist(ei.frameNumbersForRegistrationImageFile)
        display('Select the frame numbers for making registration image and then save a file named frameNumbersForRegistrationImage.txt with frame numbers in it');
        return;
    else
        fns = load(ei.frameNumbersForRegistrationImageFile);
    end
end

if exist(ei.registerationImageFile)
    display('Registration Image File already exists');
    return;
end

row = ei.pixelY;
col = ei.pixelX;
frameNumbers = fns(1):fns(2);
hWaitBar = waitbar(0,sprintf('Processed Frames -'));
fid = fopen(ei.rawFile,'r');
sumI = zeros(row*col,1);
for ii = 1:length(frameNumbers)
    frameNumber = frameNumbers(ii);
    fseek(fid, (frameNumber-1)*row*col*2, 'bof');
    I1 = fread(fid,row*col,'uint16=>double',0,'l'); 
    waitbar(ii/length(frameNumbers),hWaitBar,sprintf('Processing Frame %d',frameNumber));
    sumI = sumI + I1;
end
fclose(fid);
close(hWaitBar) 
sumI = sumI/length(frameNumbers);
Z = reshape(sumI,row,col); Z = Z';
m = matfile(ei.registerationImageFile,'Writable',true);
m.registerationImage = Z;
delete(m);
figure(1);clf;
k=imagesc(Z);
colormap gray;
titleText = sprintf('Frame Number = %d',frameNumber); title(titleText);
axis equal
axis off;
display('Done making registration image');
% imwrite(uint16(Z),makeName('registerImage.tif',ei.pDataFolder),'tif');