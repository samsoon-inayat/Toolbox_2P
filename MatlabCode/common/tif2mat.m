function varargout = tif2mat (PathName,varargin)

if nargin == 2
    tifFileName = varargin{1};
    fileName = makeName(tifFileName,PathName);
else
    fileName = PathName;
end

info = imfinfo(fileName);
warning('off','MATLAB:imagesci:rtifc:notPhotoTransformed');
num_images = numel(info);
for k = 1:num_images
    ImgSeq(:,:,k) = double(imread(fileName,'tif',k));
end
warning('on','MATLAB:imagesci:rtifc:notPhotoTransformed');
if nargout == 1
    varargout{1} = ImgSeq;
else
    matFileName = makeName(tifFileName(1:(end-4)),PathName)
    cmdTxt = sprintf('%s = ImgSeq;',tifFileName(1:(end-4)));
    eval(cmdTxt);
    if strcmp(tifFileName(1:(end-4)),'Mask')
        mask = Mask;
        save(matFileName,sprintf('%s',tifFileName(1:(end-4))),'mask');
    else
        save(matFileName,sprintf('%s',tifFileName(1:(end-4))));
    end
    
end