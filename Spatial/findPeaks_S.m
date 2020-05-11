function out = findPeaks_S(Rs,selCells,nPeaks)
if isempty(selCells)
    selCells = 1:size(Rs.rasters,3);
end
for ii = 1:length(selCells)
    thisCell = selCells(ii);
    raster = Rs.rasters(:,:,thisCell);
    
    mSig = nanmean(raster);
    if nPeaks == 1
        props = place_cells_single_peak(mSig,raster,'cm_per_bin',Rs.dist(end)/size(raster,2));
    end
%     if nPeaks == 2
%         props = place_cells_two_peaks(mSig,raster,'cm_per_bin',Rs.dist(end)/size(raster,2));
%     end
    pfs{ii} = props;
end
for ii = 1:length(pfs)
    t_widths(ii) = pfs{ii}.widths;
    t_centers(ii) = pfs{ii}.centers;
    t_peaks(ii) = pfs{ii}.peaks;
end
inds = find(~isnan(t_widths));
out.widths = t_widths(inds);
out.centers = t_centers(inds);
out.peaks = t_peaks(inds);
out.cells = selCells(inds);
n = 0;