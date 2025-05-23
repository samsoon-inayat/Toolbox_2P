function ei = get_control_air_offsets(ei)
binWidths = evalin('base','binwidths');
binwidths = binWidths;
tei = ei{1};
b = tei.b;
numplanes = length(tei.plane);
thisRasterType = 'time';
owr = [0 0 0];
thisStimMarker = 'airOffsets22_C';
trials = 1:10;
for pp = 1:numplanes
    tplane = tei.plane{pp};
    thispFolder = tei.plane{pp}.folder;
    for ci = [2 3 4 5 7]
        context3 = tplane.contexts(ci);
        disp(sprintf('%s - %s',tplane.folder,context3.name));
        contextMC = tplane.contextsMC(ci);
        air_offsets = context3.markers.airOffsets11_onsets + round(1e6 * 1/tei.b.si);
        if ismember(ci,[2 7])
            onsets = air_offsets + round(1e6 * 5/tei.b.si);
        else
            onsets = air_offsets + round(1e6 * 7.5/tei.b.si);
        end
        offsets = onsets;
        toffsets = offsets + round(1e6 * 1/tei.b.si);
        tonsets = onsets - round(1e6 * 1/tei.b.si);
        markersOn = tonsets;
        markersOff = toffsets;

        rasters = make_rasters(tei,pp,markersOn,markersOff,thisRasterType,binWidths);
        rasters = findRasterProperties_1(thispFolder,ci,thisStimMarker,rasters,thisRasterType,1:10,[0 0 0]);
        
        rastersMC = make_rasters_motion_correction(tei,pp,markersOn,markersOff,'time',binwidths);
        rastersMC = findRasterProperties_1(thispFolder,ci,sprintf('%sMC',thisStimMarker),rastersMC,thisRasterType,trials,owr);
        contextMC.rasters.airOffsets11_C = rastersMC;
        tplane.contextsMC(ci) = contextMC;
        
        rasters = findRasterProperties_1_MC(thispFolder,ci,thisStimMarker,{rasters,rastersMC},thisRasterType,trials,owr);
        context3.rasters.airOffsets11_C = rasters;
        tplane.contexts(ci) = context3;
    end
    tei.plane{pp} = tplane;
end
ei{1} = tei;