function fd = findHiFD(raster,trials)
% raster = raster(trials,:);
raster(isnan(raster)) = 0;
rasterg = mat2gray(raster(trials,:));
rraster = reshape(rasterg',1,numel(rasterg));
fd = Higuchi_FD(rraster,20);