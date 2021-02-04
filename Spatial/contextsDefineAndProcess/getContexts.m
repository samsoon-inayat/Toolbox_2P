function contexts = getContexts(ei,dcfilename)

%% old development notes
% There are three or four contexts
% In each context, I want to define air onset and offset, belt markers,
% motion onset and offsets in the intertrial intervals
% so the type of data set I will have will be context definition and then
% in that context I will have different types of markers i.e. air onset and
% offset, belt markers and motion onsets and offsets

%%

if ~exist('dcfilename','var')
    dcfilename = 'define_contexts.m';
end

plotFlag = 0;
allContexts = contextDefinitions;

aa = 1
tei = ei{aa};
contexts = defineContexts(tei,dcfilename);
[motion.markersOn,motion.markersOff] = getMarkers(tei,NaN,NaN,'motion');
[motionI.markersOn,motionI.markersOff] = getMarkers(tei,NaN,NaN,'motionI',motion);
for ii = 1:length(contexts) % process all contexts
    if isempty(contexts(ii).trials)
        continue;
    end
    stimMarkers = contexts(ii).stimMarkers;
    if length(stimMarkers) == 0
        continue;
    end
    for jj = 1:length(stimMarkers)
%             [ii jj]
%             if isequal([ii jj],[4 5])
%                 n = 0;
%             end
        trials = contexts(ii).trials;
        if iscell(trials)
            trials = trials{jj};
        end
        display(sprintf('%s---%s',contexts(ii).name,stimMarkers{jj}));
        if strcmp(stimMarkers{jj},'air') | strcmp(stimMarkers{jj},'airI') | strcmp(stimMarkers{jj},'motionOnsetsOffsets') | ...
                strcmp(stimMarkers{jj},'motionOffsetAirOnset')
        end
        if strcmp(stimMarkers{jj},'motionI')
            markersOn = motionI.markersOn;
            markersOff = motionI.markersOff;
        else
            if strcmp(stimMarkers{jj},'motionOnsets22') | strcmp(stimMarkers{jj},'motionOffsets22')
                [markersOn,markersOff] = getMarkers(tei,ii,trials,stimMarkers{jj},motion);
            elseif strcmp(stimMarkers{jj},'motionOnsetsOffsets')
                [markersOn,markersOff] = getMarkers(tei,ii,trials,stimMarkers{jj},motion);
                motionOnsetsOffsets.markersOn = markersOn; motionOnsetsOffsets.markersOff = markersOff;
            elseif strcmp(stimMarkers{jj},'motionOffsetAirOnset')
                [markersOn,markersOff] = getMarkers(tei,ii,trials,stimMarkers{jj},motionOnsetsOffsets);
            else
                [markersOn,markersOff] = getMarkers(tei,contexts(ii).name,trials,stimMarkers{jj});
            end
        end
        if ~strcmp(stimMarkers{jj},'airI') & ~strcmp(stimMarkers{jj},'motionOffsetAirOnset') & ~strcmp(contexts(ii).name,'Light - Brake') & ~strcmp(contexts(ii).name,'Air - Brake')...
                & isempty(strfind(stimMarkers{jj},'light')) & isempty(strfind(stimMarkers{jj},'tone')) & ~strcmp(stimMarkers{jj},'airOffsets22') & ~strcmp(stimMarkers{jj},'airOnsets22')...
                & ~strcmp(stimMarkers{jj},'airOffsets27') & ~strcmp(stimMarkers{jj},'airOnsets27')...
                & ~strcmp(stimMarkers{jj},'airOffsets11') & ~strcmp(stimMarkers{jj},'airOnsets11')...
                & ~strcmp(stimMarkers{jj},'airOffsets01') & ~strcmp(stimMarkers{jj},'airOnsets01')...
                & ~strcmp(stimMarkers{jj},'airOffsets010') & ~strcmp(stimMarkers{jj},'airOnsets010')
%                     [markersOn,markersOff] = validateForFrames(tei,markersOn,markersOff);
        end
        cmdTxt = sprintf('contexts(ii).markers.%s_onsets = markersOn;',stimMarkers{jj});
        eval(cmdTxt);
        cmdTxt = sprintf('contexts(ii).markers.%s_offsets = markersOff;',stimMarkers{jj});
        eval(cmdTxt);
%             if strcmp(stimMarkers{jj},'belt')
%                 plotMarkers(tei.b,markersOn,markersOff,101,0);
%             else
%                 plotMarkers(tei.b,markersOn,markersOff,101,1);
%             end
        if strcmp(stimMarkers{jj},'motionOffsetAirOnset')
            plotMarkers(tei.b,markersOn,markersOff,101,1);
        end
        if plotFlag
            plotMarkers(tei.b,markersOn,markersOff,101,0);
            title(tei.recordingFolder);
            pause(0.1);
        end
% %             ylim([0 40]);
        if jj == 3 && ii == 4
            n = 0;
        end
    end
end
n = 0;
S.contexts = contexts;
%         S.typesOfMarkers = allContexts.typesOfMarkers;


function trials = adjustTrials(trials,adjust)
trials(1) = trials(1) + adjust(1);
trials(end) = trials(end) + adjust(2);

function [markersOn,markersOff] = validateForFrames(tei,markersOn,markersOff)
b = tei.b;
allFrames = b.frames_f;
inds = [];
for ii = 1:length(markersOn)
    st = markersOn(ii);
    se = markersOff(ii);
    ns = find(allFrames > st & allFrames < se);
    if isempty(ns)
        inds = [inds ii];
        continue;
    end
    d = b.dist(st:se)-b.dist(st);
    if d(end) < 20
        inds = [inds ii];
        continue;
    end
end
markersOn(inds) = [];
markersOff(inds) = [];
n = 0;