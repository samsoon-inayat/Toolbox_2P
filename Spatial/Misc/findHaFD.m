function fd = findHaFD(raster,trials)
% raster = raster(trials,:);
raster(isnan(raster)) = 0;
bw = imbinarize(raster(trials,:));
fd = BoxCountfracDim(bw);