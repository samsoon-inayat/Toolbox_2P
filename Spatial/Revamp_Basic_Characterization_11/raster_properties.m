function data = raster_properties(rasters)

for ii = 1:size(rasters,3)
    thisRaster = rasters(:,:,ii);
    HaFD(ii) = findHaFD(thisRaster,1:size(thisRaster,1));
    HiFD(ii) = findHiFD(thisRaster,1:size(thisRaster,1));
    Ent(ii) = entropy(mat2gray(thisRaster));
    temp = info_metrics_S_onlyMI(thisRaster,[],4,[],0);
    MI(ii) = temp.ShannonMI;
end
n = 0;
data.MI = MI;
data.Ent = Ent;
data.HaFD = HaFD;
data.HiFD = HiFD;


