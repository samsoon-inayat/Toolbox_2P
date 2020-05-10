function ei = loadContextsResponses(ei,owr,owrp)

allContexts = contextDefinitions;
for aa = 1:length(ei)
    ptei = ei{aa};
    tei = ptei;
    for pp = 1:length(ei{aa}.plane)
        thispFolder = tei.plane{pp}.folder;
        disp(thispFolder);
        tei.folders.thispFolder = thispFolder;
%         tei.areCells = find(tei.plane{pp}.tP.areCells);
        tei.deconv = tei.plane{pp}.tP.deconv;
        tei.ops1 = tei.plane{pp}.tP.ops;
        tei.tP = tei.plane{pp}.tP;
        tplane = tei.plane{pp};
        fileName = makeName('contexts.mat',thispFolder);
        temp = load(fileName);
        if ~isfield(tplane,'contexts')
            fileName = makeName('contexts_trials_markers.mat',thispFolder);
            S = load(fileName);
            contexts = S.contexts;
        else
            contexts = tplane.contexts;
        end
%         typesOfMarkers = S.typesOfMarkers;
        for contextNumber = 1:length(contexts)
            thisContext = contexts(contextNumber);
            stimMarkers = thisContext.stimMarkers;
            trials = thisContext.trials;
            if isempty(trials)
                continue;
            end
            if length(stimMarkers) == 0
                continue;
            end
            for jj = 1:length(stimMarkers)
                thisStimMarker = stimMarkers{jj};
                cmdTxt = sprintf('markersOn = contexts(contextNumber).markers.%s_onsets;',thisStimMarker);
                eval(cmdTxt);
                cmdTxt = sprintf('markersOff = contexts(contextNumber).markers.%s_offsets;',thisStimMarker);
                eval(cmdTxt);
                if isempty(markersOn) & isempty(markersOff)
                    continue;
                end
                typesOfRasters = getTypesOfRasters(allContexts,thisStimMarker);
                for kk = 1:length(typesOfRasters)
                    thisRasterType = typesOfRasters{kk};
                    if strcmp(thisRasterType,'dist') && owr(1) == 0
                        continue;
                    end
                    if strcmp(thisRasterType,'time') && owr(2) == 0
                        continue;
                    end
                    disp(sprintf('%s--- %s -- %s',thisContext.name,stimMarkers{jj},thisRasterType));
                    rasters = getRasters_fixed_bin_width(tei,pp,markersOn,markersOff,thisRasterType);
                    trials = 4:size(rasters.sp_rasters,1);
                    rasters = findRasterProperties(thispFolder,contextNumber,thisStimMarker,rasters,thisRasterType,trials,owrp);
                    if strcmp(thisRasterType,'dist')
                        cmdTxt = sprintf('contexts(contextNumber).rasters.%sD = rasters;',thisStimMarker); eval(cmdTxt);
                    end
                   if strcmp(thisRasterType,'time')
                        cmdTxt = sprintf('contexts(contextNumber).rasters.%sT = rasters;',thisStimMarker); eval(cmdTxt);
                    end 
                end
            end
        end
        ptei.plane{pp}.contexts = contexts;
        ptei.b.belt_length = temp.belt_length;
    end
    ei{aa} = ptei;
end

