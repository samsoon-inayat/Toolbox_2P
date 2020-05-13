function PreProcessEV_HiLo(file,filepath,savePath,Fs,Rotate,theta)
trial = file(3:end);

% savePath = filepath;
% if nargin > 2
%     savePath = varargin{1};
% end

% Rotate = 0;
% theta = 0;
% if nargin > 3
%     Rotate = varargin{2};
% end
% if Rotate && nargin > 4
%     theta = varargin{3};
% end

% Fs = 150; % sampling rate
% 2D Gaussian
param.GaussWindow = 5;
param.Sigma       = 1;
Center  = fix([param.GaussWindow/2,param.GaussWindow/2])+1;
[R,C]   = ndgrid(1:param.GaussWindow, 1:param.GaussWindow);
Gauss2d = gauss2dC(R,C,param.Sigma,Center);
% Load ImgSeq_raw
filetype = '.tif';
filenameImgSeq = [filepath, file, filetype];
ImgSeq_raw = imreadalltiff(filenameImgSeq); ImgSeq_raw = double(ImgSeq_raw);
% ImgSeq_raw(:,:,32) = (ImgSeq_raw(:,:,31)+ImgSeq_raw(:,:,33))/2;
[width,height,nFrames] = size(ImgSeq_raw);
% Load no trials and use it to calc F0
filetype = '.tif';
filenameNo = [filepath, 'no',trial, filetype];
F0 = imreadalltiff(filenameNo); F0 = double(F0);
for r = 1:width
    for c = 1:height
        sig = F0(r,c,:);sig =sig(:);
        yhat = sig - locdetrend(sig,Fs,[0.3,0.1]);
        F0(r,c,:) = yhat;
    end
end
% Calc DF/F0
ImgSeq = 100*(ImgSeq_raw-F0)./F0; ImgSeq(isnan(ImgSeq))=0; ImgSeq(isinf(ImgSeq))=0;
% LPF
Gauss1d = LPF_Gauss_25Hz; Lhalf = round((length(Gauss1d)-1)/2);
for r = 1:width
    for c = 1:height
        sig = ImgSeq(r,c,:); sig = sig(:);
        ImgSeq(r,c,:) = conv(sig,Gauss1d,'same');
    end
end
% Gaussian spatial filtering
for idx = 1:nFrames
    ImgSeq(:,:,idx) = conv2(ImgSeq(:,:,idx),Gauss2d, 'same');
    ImgSeq(:,:,idx) = conv2(ImgSeq(:,:,idx),Gauss2d, 'same');
end
% Trimm the filtering edge effect
ImgSeq = ImgSeq(:,:,Lhalf+1:nFrames-Lhalf);
if Rotate
    n90 = round(theta/90);
    for idx = 1:size(ImgSeq,3)
        im = ImgSeq(:,:,idx);
        im = rot90(im,n90);
        ImgSeq(:,:,idx) = im;
    end
end
dfbyfo = ImgSeq;
% save the results
% imwriteallraw([savePath, file, '_DF_F0_LPF_25Hz_G1_G1.raw'],ImgSeq,'*float32');
save([savePath, 'dfbyfo.mat'],'dfbyfo','-v7.3');
end
