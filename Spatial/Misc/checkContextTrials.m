function checkContextTrials(ei,contextNumber)

for ii = 1:length(ei)
    tei = ei{ii};
    contexts = tei.plane{1}.contexts(contextNumber);
    plotMarkers(tei.b,contexts.markers.air_onsets,contexts.markers.air_offsets,101,0);
    title(tei.recordingFolder);
    pause;
    n = 0;
end

