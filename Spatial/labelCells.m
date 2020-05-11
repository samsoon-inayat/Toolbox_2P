function labelCells (ei,coi,SIs)
mimg = ei.plane{1}.ops1{1}.meanImg;
cimg = zeros(size(mimg));
centroids = [];
for ii = 1:length(coi)
    temp = zeros(size(mimg));
    stat = ei.plane{1}.tP.stat(coi(ii))
    xpix = stat{1}.xpix;
    ypix = stat{1}.ypix;
    thisPixelIdxs = sub2ind(size(mimg),ypix,xpix);
    cimg(thisPixelIdxs) = SIs(ii)*ones(size(thisPixelIdxs));
    temp(thisPixelIdxs) = ones(size(thisPixelIdxs));
    s = regionprops(temp,'centroid');
    centroids = [centroids;round(s.Centroid)];
end
figure(1000);clf;
% subplot(1,2,1);
% imagesc(mimg);axis equal
% subplot(1,2,2);
imagesc(cimg);axis equal
for ii = 1:length(coi)
    text(centroids(ii,1),centroids(ii,2),num2str(coi(ii)));
end
colorbar