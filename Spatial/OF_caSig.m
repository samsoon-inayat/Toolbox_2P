function OF_frames = OF_caSig(ei,coi)

mimg = ei.ops1{1}.mimg1;
cimg = zeros(size(mimg));
centroids = [];
for ii = 1:length(coi)
    temp = zeros(size(mimg));
    xpix = ei.tP.stat(coi(ii)).xpix;
    ypix = ei.tP.stat(coi(ii)).ypix;
    thisPixelIdxs = sub2ind(size(mimg),ypix,xpix);
    cimg(thisPixelIdxs) = ones(size(thisPixelIdxs));
    temp(thisPixelIdxs) = ones(size(thisPixelIdxs));
    s = regionprops(temp,'centroid');
    centroids = [centroids;round(s.Centroid)];
end
nFrames = size(ei.signals,2);
frames = zeros(size(cimg,1),size(cimg,2),nFrames);
for ii = 1:length(coi)
    frames(centroids(ii,1),centroids(ii,2),:) = ei.signals(coi(ii),:);
end

alpha = 0.012;    ratio = 0.75;    minWidth = 20;    nOuterFPIterations = 7;    nInnerFPIterations = 1;    nSORIterations = 30;
para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];
tic;
for ii = 2:nFrames
    ii
    previousFrame = frames(:,:,ii-1);
    thisFrame = frames(:,:,ii);
    [u,v,~] = Coarse2FineTwoFrames(previousFrame,thisFrame,para);
    OF_Frames(:,:,ii-1) = u.*(i*v);
end
toc
n = 0;