function ei = findSpeedRasters(ei,contextNs,ow)

if ~exist('ow','var')
    ow = 0;
end

for ii = 1:length(ei)
    aa = ii;
    % extract rasters from ei
    for pp = 1%:length(ei{aa}.plane)
        thispFolder = ei{aa}.plane{pp}.folder;
        disp(sprintf('Processing context %s',thispFolder));
        contexts = ei{ii}.plane{pp}.contexts;
        for jjj = 1:length(contextNs)
            jj = contextNs(jjj);
            thisContext = contexts(jj);
            onsets = thisContext.markers.air_onsets - round(0e6 * 1/ei{ii}.b.si);
            offsets = thisContext.markers.air_offsets + round(10e6 * 1/ei{ii}.b.si);
            rasters =  getSpeedRasters(ei{ii},onsets,offsets,ow,0);
            contexts(jj).rastersSpeed = rasters;
%             contexts(jj) = thisContext;
        end
        ei{ii}.plane{pp}.contexts = contexts;
    end
end
display('Found Speed Rasters ');