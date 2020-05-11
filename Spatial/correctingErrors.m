function ei = correctingErrors(ei,owr,owrp)

if ~exist('ei','var')
    ei = evalin('base','ei10');
%     ei1 = evalin('base','ei');
    owr = 0; owrp = [0 0 0];
end


allContexts = contextDefinitions;
for aa = 1:length(ei)
    ptei = ei{aa};
    tei = ptei;
%     b = tei.b;
%     b1 = ei1{aa}.b;
%     figure(100);clf;plot(b.ts,b.fSpeed);hold on;plot(b1.ts,b1.fSpeed);
%     pause;
%     continue;
    for pp = 1:length(ei{aa}.plane)
        thispFolder = tei.plane{pp}.folder;
        files = dir(sprintf('%s/tRaster*.mat',thispFolder));
        n = 0;
        for ff = 1:length(files)
            fileName = fullfile(thispFolder,files(ff).name);
%             delete(fileName);
            disp(fileName);
        end
        n = 0;
        continue;
        disp(thispFolder);
        tei.folders.thispFolder = thispFolder;
%         tei.areCells = find(tei.plane{pp}.tP.areCells);
        tei.deconv = tei.plane{pp}.tP.deconv;
        tei.ops1 = tei.plane{pp}.tP.ops;
        tei.tP = tei.plane{pp}.tP;
        fileName = makeName('contexts.mat',thispFolder);
        temp = load(fileName);
        fileName = makeName('contexts_trials_markers.mat',thispFolder);
        S = load(fileName);
        contexts = S.contexts;
        typesOfMarkers = S.typesOfMarkers;
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
                display(sprintf('%s---%s',thisContext.name,stimMarkers{jj}));
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
                    thisRaster = typesOfRasters{kk};
                    if strcmp(thisRaster,'dist')
                        [dRasters,fileName] = getDistRasters(tei,markersOn,markersOff,owr,0);
%                         disp(fileName);
%                         delete(fileName);
                        n = 0;
                    end
                    if strcmp(thisRaster,'time')
                        [tRasters,fileName] = getTimeRasters(tei,markersOn,markersOff,owr,0);
%                         disp(fileName);
%                         delete(fileName);
                        n = 0;
                    end
                end
            end
        end
        ptei.plane{pp}.contexts = contexts;
        ptei.b.belt_length = temp.belt_length;
    end
%     ei{aa} = ptei;
end

