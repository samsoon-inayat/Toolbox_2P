function ImgSeq = applySpatialFilter(ImgSeq)

param.GaussWindow = 5;
param.Sigma       = 1;
Center  = fix([param.GaussWindow/2,param.GaussWindow/2])+1;
[R,C]   = ndgrid(1:param.GaussWindow, 1:param.GaussWindow);
Gauss2d = gauss2dC(R,C,param.Sigma,Center);
nFrames = size(ImgSeq,3);
parfor idx = 1:nFrames
    ImgSeq(:,:,idx) = conv2(ImgSeq(:,:,idx),Gauss2d, 'same');
    ImgSeq(:,:,idx) = conv2(ImgSeq(:,:,idx),Gauss2d, 'same');
end