function PreProcessSpon(filepath, file, varargin)

Rotate = 0;
theta = 0;
if nargin > 2
    Rotate = varargin{1};
end
if Rotate && nargin > 3
    theta = varargin{2};
end
    
Fs = 150; % sampling rate
% 2D Gaussian
param.GaussWindow = 5;
param.Sigma       = 1;
Center  = fix([param.GaussWindow/2,param.GaussWindow/2])+1;
[R,C]   = ndgrid(1:param.GaussWindow, 1:param.GaussWindow);
Gauss2d = gauss2dC(R,C,param.Sigma,Center);
% Load ImgSeq
filetype = '.tif';
filename = [filepath, file, filetype];
ImgSeq = imreadalltiff(filename); ImgSeq = double(ImgSeq); 
[width,height,nFrames] = size(ImgSeq);
% calc F0 using locdetrend
F0 = zeros(size(ImgSeq));
for r = 1:width
    for c = 1:height
        sig = ImgSeq(r,c,:);sig =sig(:);
        yhat = sig - locdetrend(sig,Fs,[2,1]);
        F0(r,c,:) = yhat;
    end
end
DF_F0 = 100*(ImgSeq-F0)./F0; DF_F0(isnan(DF_F0))=0; DF_F0(isinf(DF_F0))=0;
% LPF
Gauss1d = LPF_Gauss_25Hz;
for r = 1:width
    for c = 1:height
        sig = DF_F0(r,c,:); sig = sig(:);
        DF_F0(r,c,:) = conv(sig,Gauss1d,'same');
    end
end
% Gaussian spatial filtering
for idx = 1:nFrames
    DF_F0(:,:,idx) = conv2(DF_F0(:,:,idx),Gauss2d, 'same');
    DF_F0(:,:,idx) = conv2(DF_F0(:,:,idx),Gauss2d, 'same');
end
if Rotate
    n90 = round(theta/90);
    for idx = 1:size(DF_F0,3)
        im = DF_F0(:,:,idx);
        im = rot90(im,n90);
        DF_F0(:,:,idx) = im;
    end
end
% save the results
% imwriteallraw([filepath, file, '_DF_F0_LPF_25Hz_G1_G1.raw'],DF_F0,'*float32');
save([filepath, file, '_DF_F0_LPF_25Hz_G1_G1.mat'],'DF_F0','-v7.3');
