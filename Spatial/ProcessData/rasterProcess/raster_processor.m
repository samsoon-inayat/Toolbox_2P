function ei = raster_processor(ei,selContexts,stimMarker,trials,ow)

if ~exist('ow','var')
    ow = 0;
end

for ii = 1:length(ei)
    % extract rasters from ei
    cmdTxt = sprintf('Rs = ei{ii}.contexts(selContexts).rasters.%sD;',stimMarker);
    eval(cmdTxt);
    fileName = fullfile(ei{ii}.folders.thispFolder,sprintf('mean_raster_fits_%d_%d_%d_%s.mat',selContexts,stim));
    mean_raster_fits = getFits(Rs,1:size(Rs.rasters,3),trials);
end
