function find_trial_distance_time_distributions(ei,owr)

if ~exist('ei','var')
    ei = evalin('base','ei10');
    owr = 0;
end
ind = 1;
allContexts = contextDefinitions;
for aa = 1:length(ei)
    ptei = ei{aa};
    tei = ptei;
    b = tei.b;
    for pp = 1:length(ei{aa}.plane)
        thispFolder = tei.plane{pp}.folder;
        b.frames_f = tei.plane{pp}.b.frames_f;
        disp(thispFolder);
        tei.folders.thispFolder = thispFolder;
%         tei.areCells = find(tei.plane{pp}.tP.areCells);
        tei.deconv = tei.plane{pp}.tP.deconv;
        tei.ops1 = tei.plane{pp}.tP.ops;
        tei.tP = tei.plane{pp}.tP;
        frameRates(aa,pp) = 1/(tei.b.ts(tei.plane{1}.b.frames_f(2)) - tei.b.ts(tei.plane{1}.b.frames_f(1)));
        tplane = tei.plane{pp};
        fileName = makeName('contexts.mat',thispFolder);
        temp = load(fileName);
        if owr == 0 || owr == 1
            fileName = makeName('contexts_trials_markers.mat',thispFolder);
            S = load(fileName);
            contexts = S.contexts;
        else
            contexts = tplane.contexts;
        end
%         typesOfMarkers = S.typesOfMarkers;
        if pp > 1
            continue;
        end
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
                if strcmp(thisStimMarker,'air')
                    [minDist(ind),minTime(ind)] =  get_min_dist_time(b,markersOn,markersOff);
                    ind = ind + 1;
                end
                n = 0;
%                 typesOfRasters = getTypesOfRasters(allContexts,thisStimMarker);
%                 for kk = 1:length(typesOfRasters)
%                     thisRaster = typesOfRasters{kk};
%                     if owr == 0
%                         if strcmp(thisRaster,'dist')
%                             dRasters = getDistRasters(tei,pp,markersOn,markersOff,owr,0);
%                         end
%                         if strcmp(thisRaster,'time')
%                             tRasters = getTimeRasters(tei,pp,markersOn,markersOff,owr,0);
%                         end
%                     end
%                 end
            end
        end
        ptei.plane{pp}.contexts = contexts;
        ptei.b.belt_length = temp.belt_length;
    end
    ei{aa} = ptei;
end

n = 0;



function [minDist,minTime] =  get_min_dist_time(b,onsets,offsets)
trialNumbers = 1:length(onsets);
for iii = 1:length(trialNumbers)
    ii = trialNumbers(iii);
    st = onsets(ii);
    se = offsets(ii);
%     frames = (find(b.frames_f >= st & b.frames_f <= se)); % find frame numbers
%     oSignals{iii} = caSig(frames);
%     spoSignals{iii} = spSignal(frames);
%     maxSignals(iii) = max(oSignals{iii});
%     minSignals(iii) = min(oSignals{iii});
%     ts{iii} = b.ts(b.frames_f(frames))-b.ts(st);
    temp = b.fSpeed(st:se);
    temp(temp<0) = 0;
    speed{iii} = temp;
    tss{iii} = b.ts(st:se)-b.ts(st);
    dist{iii} = b.dist(st:se)-b.dist(st);
    encoderCount(iii) = b.encoderCount(se) - b.encoderCount(st);
    distVal(iii) = dist{iii}(end);
    timeVal(iii) = tss{iii}(end);
%     if isfield(b,'tone_light_stim')
%         tempTL = (find(b.tone_light_stim_r >= st & b.tone_light_stim_r <= se)); % find frame numbers
%         try
%             tone_light{iii} = tempTL(3:4);
%         catch
%             tone_light{iii} = tempTL(1:2);
%         end
%         tstl{iii} = b.ts(b.tone_light_stim_r(tone_light{iii}))-b.ts(st);
%         dstl(iii,:) = b.dist(b.tone_light_stim_r(tone_light{iii}))-b.dist(st);
%     end
end
minDist = min(distVal);
minTime = min(timeVal);


