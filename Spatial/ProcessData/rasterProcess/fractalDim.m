function out = fractalDim(thispFolder,contextNumber,stimMarker,Rs,thisRasterType,trials,ow)


fileName = makeName(sprintf('factal_dim_%d_%s_%s.mat',contextNumber,stimMarker,thisRasterType),thispFolder);
if exist(fileName,'file') && ow == 0
    out = load(fileName);
    return;
end
hWaitBar = waitbar(0,sprintf('Finding fractal dim ... plz wait'));

rasters = Rs.sp_rasters_nan_corrected;
numCells = size(rasters,3);

for ii = 1:numCells
    waitbar(ii/numCells,hWaitBar,sprintf('Finding fractal dim - Processing cell %d/%d',ii,numCells));
    raster = rasters(:,:,ii);
    try
        out.HaFD(ii) = findHaFD(raster,trials);
    catch
%         disp('Error finding fractal dim Ha')
        out.HaFD(ii) = NaN;
    end
    try
        out.HiFD(ii) = findHiFD(raster,trials);
    catch
%         disp('Error finding fractal dim Ha');
        out.HiFD(ii) = NaN;
    end
end
close(hWaitBar);
save(fileName,'-struct','out','-v7.3');







