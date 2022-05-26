function ei = get_control_air_onsets(ei)
binWidths = evalin('base','binwidths');
binwidths = binWidths;
tei = ei{1};
b = tei.b;
numplanes = length(tei.plane);
thisRasterType = 'time';
owr = [0 0 0];
thisStimMarker = 'airOnsets22P';
trials = 1:10;
for pp = 1:numplanes
    tplane = tei.plane{pp};
    thispFolder = tei.plane{pp}.folder;
    for ci = [2 3 4 5 7]
        context3 = tplane.contexts(ci);
        disp(sprintf('%s - %s',tplane.folder,context3.name));
        contextMC = tplane.contextsMC(ci);
        air_onsets = context3.markers.airOnsets11_onsets + round(1e6 * 1/tei.b.si);
        onsets = air_onsets - round(1e6 * 4.5/tei.b.si);
        offsets = onsets;
        toffsets = offsets + round(1e6 * 1/tei.b.si);
        tonsets = onsets - round(1e6 * 1/tei.b.si);
        markersOn = tonsets;
        markersOff = toffsets;

        rasters = make_rasters(tei,pp,markersOn,markersOff,thisRasterType,binWidths);
        rasters = findRasterProperties_1(thispFolder,ci,thisStimMarker,rasters,thisRasterType,1:10,[0 0 0]);
        
        rastersMC = make_rasters_motion_correction(tei,pp,markersOn,markersOff,'time',binwidths);
        rastersMC = findRasterProperties_1(thispFolder,ci,sprintf('%sMC',thisStimMarker),rastersMC,thisRasterType,trials,owr);
        contextMC.rasters.airOnsets11P = rastersMC;
        tplane.contextsMC(ci) = contextMC;
        
        rasters = findRasterProperties_1_MC(thispFolder,ci,thisStimMarker,{rasters,rastersMC},thisRasterType,trials,owr);
        context3.rasters.airOnsets11P = rasters;
        tplane.contexts(ci) = context3;
    end
    tei.plane{pp} = tplane;
end
ei{1} = tei;