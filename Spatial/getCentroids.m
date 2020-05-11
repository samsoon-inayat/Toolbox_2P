function centroids = getCentroids (ei,coi)
mimg = ei.ops1{1}.mimg1;
centroids = [];
for ii = 1:length(coi)
    temp = zeros(size(mimg));
    xpix = ei.tP.stat(coi(ii)).xpix;
    ypix = ei.tP.stat(coi(ii)).ypix;
    thisPixelIdxs = sub2ind(size(mimg),ypix,xpix);
    temp(thisPixelIdxs) = ones(size(thisPixelIdxs));
    s = regionprops(temp,'centroid');
    centroids = [centroids;round(s.Centroid)];
end
