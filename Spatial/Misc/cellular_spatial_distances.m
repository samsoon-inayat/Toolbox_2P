function ei = cellular_spatial_distances(ei)
tei = ei{1};
micronsPerPixel = tei.thorExp.widthUM/tei.thorExp.pixelX;

n = 0;

pl = tei.plane{1};
tP = pl.tP;

[ops,astat,arecells] = get_ops(tei,1);
cellinds = find(arecells(:,1)==1);
xrange = ops.xrange;
yrange = ops.yrange;
mimg = ops.meanImgE;
maskZ = zeros(size(mimg));
%%
figure(100);clf;imagesc(mimg);hold on;
for ii = 1:length(cellinds)
    cn = cellinds(ii);
    tstat = astat{cn};
    mask = maskZ;
    mask(tstat.ipix) = 1;
    mask = mask';
    rp = regionprops(mask,'centroid');
    centroids(ii,:) = rp.Centroid;
    plot(centroids(ii,1),centroids(ii,2),'.r');
%     figure(100);clf;imagesc(mimg+mask);
%     pause;
end
%%
combs = nchoosek(1:length(cellinds),2);
parfor ii = 1:size(combs,1)
    edist(ii,1) = sqrt(sum((centroids(combs(ii,2),:) - centroids(combs(ii,1),:)).^2));
end
%%
tei.plane{1}.tP.centroids = centroids;
tei.plane{1}.tP.dists_centroids = [combs edist];
ei{1} = tei;