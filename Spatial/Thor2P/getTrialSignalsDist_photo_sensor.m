function [thisDist,oSignals,trials,mInfo] =  getTrialSignalsDist_photo_sensor(signal,b,onset,offset)
trialNumbers = 1:length(b.sTrials);
nBins = 100;
increment = 2300/nBins;
sBins = 1:increment:2300;
eBins = increment:increment:2300;
thisDist = (1:nBins)*increment*2*pi*5/500;
GF = gausswin(3); % 3 bins means around 4.5 cm
tns = 0;
trials = [];
mInfo = [];
for iii = 1:length(trialNumbers)
    if iii == 12
        n = 0;
    end
    ii = trialNumbers(iii);
    try
        st = onset(ii);
        se = offset(ii);
    catch
        error;
    end
    frames = (find(b.frames_f >= st & b.frames_f <= se)); % find frame numbers
    dists_pulses = b.ch_a_r(find(b.ch_a_r>= st & b.ch_a_r<= se)); % find encoder pulses
    ldp(iii) = length(dists_pulses);
    if length(dists_pulses) < 2300
        continue;
    end
    tns = tns + 1;
    trials = [trials iii];
    thisSignal = []; thisSpeed = [];
    for jj = 1:nBins
        sB = dists_pulses(sBins(jj));
        eB = dists_pulses(eBins(jj));
        bFrames = find(b.frames_f >= sB & b.frames_f<= eB);
        dDur = b.ts(eB) - b.ts(sB);
        if isempty(bFrames)
            thisSignal(jj) = 0;
        else
            thisSignal(jj) = max(signal(bFrames))/dDur;
        end
        thisSpeed(jj) = nBins/((eB-sB)*b.si*1e-6);
        mInfo.trialDurations(tns,jj) = dDur;
        if isnan(thisSignal(jj))
            n = 0;
        end
    end
    nanPos = find(isnan(thisSignal));
    for nn = 1:length(nanPos)
        nnn = nanPos(nn);
        if nnn == 1 & ~isnan(thisSignal(2))
            if isnan(thisSignal(2))
                error;
            end
            thisSignal(nnn) = thisSignal(2);
            continue;
        end
        if nnn == length(thisSignal) & ~isnan(thisSignal(length(thisSignal)-1))
            if isnan(thisSignal(length(thisSignal)-1))
                error;
            end
            thisSignal(nnn) = thisSignal(length(thisSignal)-1);
            continue;
        end
        if ~isnan(thisSignal(nnn+1))
            thisSignal(nnn) = mean([thisSignal(nnn-1) thisSignal(nnn+1)]);
        else
            thisSignal(nnn) = mean([thisSignal(nnn-1) thisSignal(nnn+2)]);
        end
    end
    if sum(isnan(thisSignal)) > 0
        n = 0;
    end
    thisSignal = thisSignal-min(thisSignal);
    thisSignal = [0 0 thisSignal];
    thisSignal = filter(GF,1,thisSignal);
    oSignals(tns,:) = thisSignal(3:end);
    mInfo.speeds(tns,:) = thisSpeed;
end


