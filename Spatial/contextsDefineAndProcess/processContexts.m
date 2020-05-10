function processContexts (ei,ow)
% ow = 1;

if ~exist('ei','var')
    ei = evalin('base','ei([1:2])');
end

allContexts = contextDefinitions;

for aa = 1:length(ei) % process each recording
    tei = ei{aa};
    for pp = 1:length(ei{aa}.plane)
        thispFolder = tei.plane{pp}.folder;
        tei.folders.thispFolder = thispFolder;
%         tei.areCells = find(tei.plane{pp}.tP.areCells);
        tei.deconv = tei.plane{pp}.tP.deconv;
        tei.ops1 = {tei.plane{pp}.tP.ops};
        tei.tP = tei.plane{pp}.tP;
        tei.b.frames_f = tei.plane{pp}.b.frames_f;
        disp(sprintf('Processing context %s',thispFolder));
        getDistanceBetweenCells(tei);
        fileName = makeName('contexts_trials_markers.mat',thispFolder);
        S = load(fileName);
        contexts = S.contexts;
        typesOfMarkers = S.typesOfMarkers;
        for ii = 1:length(contexts)
            thisContext = contexts(ii);
            stimMarkers = thisContext.stimMarkers;
            trials = thisContext.trials;
            if isempty(trials)
                continue;
            end
            if length(stimMarkers) == 0
                continue;
            end
            for jj = 1:length(stimMarkers)
                display(sprintf('%s---%s',thisContext.name,stimMarkers{jj}));
                thisStimMarker = stimMarkers{jj};
                cmdTxt = sprintf('markersOn = contexts(ii).markers.%s_onsets;',thisStimMarker);
                eval(cmdTxt);
                cmdTxt = sprintf('markersOff = contexts(ii).markers.%s_offsets;',thisStimMarker);
                eval(cmdTxt);
                if isempty(markersOn) & isempty(markersOff)
                    continue;
                end
                typesOfRasters = getTypesOfRasters(allContexts,thisStimMarker);
                for kk = 1:length(typesOfRasters)
                    thisRaster = typesOfRasters{kk};
                    if strcmp(thisRaster,'dist')
                        dRasters = getDistRasters(tei,markersOn,markersOff,ow,0);
                    end
                    if strcmp(thisRaster,'time')
                        tRasters = getTimeRasters(tei,markersOn,markersOff,ow,0);
                    end
                end
                n = 0;
            end
        end
        n = 0;
    end
end
% January 25, 2020, old code might delete later 
%                 if (strcmp(thisStimMarker,'air') && strcmp(contexts(ii).name,'Air - Brake')) || strcmp(thisStimMarker,'light') || strcmp(thisStimMarker,'tone')
%                     tRasters = getTimeRasters(tei,markersOn,markersOff,ow,0);
%                     continue;
%                 end
%                 if strcmp(thisStimMarker,'air') || strcmp(thisStimMarker,'belt') || strcmp(thisStimMarker,'motionI') || strcmp(thisStimMarker,'motionOnsetsOffsets')
%                     if strcmp(thisStimMarker,'motionI')
%                         n = 0;
%                     end
%                     dRasters = getDistRasters(tei,markersOn,markersOff,ow,0);
%                     if strcmp(thisContext.name,'Air - Brake');
%                         tRasters = getTimeRasters(tei,markersOn,markersOff,ow,0);
%                     end
%                 end
%                 if strcmp(thisStimMarker,'motionOnsets22') || strcmp(thisStimMarker,'motionOffsets22') || strcmp(thisStimMarker,'airOnsets22') ...
%                         || strcmp(thisStimMarker,'airOffsets22') || strcmp(thisStimMarker,'airI') || strcmp(thisStimMarker,'motionOffsetAirOnset')...
%                         || strcmp(thisStimMarker,'airOnsets27') || strcmp(thisStimMarker,'airOffsets27') ...
%                         || strcmp(thisStimMarker,'airOnsets11') || strcmp(thisStimMarker,'airOffsets11') ...
%                         || strcmp(thisStimMarker,'airOnsets01') || strcmp(thisStimMarker,'airOffsets01') ...
%                         || strcmp(thisStimMarker,'airOnsets010') || strcmp(thisStimMarker,'airOffsets010') 
%                     if strcmp(thisStimMarker,'motionOffsetAirOnset')
%                         n = 0;
%                     end
%                     tRasters = getTimeRasters(tei,markersOn,markersOff,ow,0);
%                 end