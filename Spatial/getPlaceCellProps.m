function [pcws,pccs] = getPlaceCellProps(Rs,selCells,trials)

if ~exist('trials','var')
    trials = 1:size(Rs.rasters,1);
end

for ii = 1:length(selCells)
    thisCell = selCells(ii);
    raster = Rs.rasters(trials,:,thisCell);
    mSig = nanmean(raster);
    props = place_cell_properties(mSig,raster,'cm_per_bin',142/50);
    pcws(ii) = props.width;
    pccs(ii) = props.center;
end


n = 0;