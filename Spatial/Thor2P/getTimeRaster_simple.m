function out =  getTimeRaster_simple(b,caSig,spSignal,onsets,offsets,plotFlag)
trialNumbers = 1:length(onsets);
for iii = 1:length(trialNumbers)
    ii = trialNumbers(iii);
    st = onsets(ii);
    se = offsets(ii);
    frames = (find(b.frames_f >= st & b.frames_f <= se)); % find frame numbers
    oSignals{iii} = caSig(frames);
    spoSignals{iii} = spSignal(frames);
    lspoSigs(iii) = length(frames);
    maxSignals(iii) = max(oSignals{iii});
    minSignals(iii) = min(oSignals{iii});
    ts{iii} = b.ts(b.frames_f(frames))-b.ts(st);
    temp = b.speed(st:se);
    temp(temp<0) = 0;
    speed{iii} = temp;
    tss{iii} = b.ts(st:se)-b.ts(st);
    dist{iii} = b.dist(st:se)-b.dist(st);
    distVal(iii) = dist{iii}(end);
    speed{iii} = b.speed(st:se);
%     tempTL = (find(b.tone_light_stim_r >= st & b.tone_light_stim_r <= se)); % find frame numbers
%     try
%         tone_light{iii} = tempTL(3:4);
%     catch
%         tone_light{iii} = tempTL(1:2);
%     end
%     tstl{iii} = b.ts(b.tone_light_stim_r(tone_light{iii}))-b.ts(st);
%     dstl(iii,:) = b.dist(b.tone_light_stim_r(tone_light{iii}))-b.dist(st);
end
minL = min(lspoSigs);
for iii = 1:length(trialNumbers)
    thisspSign = spoSignals{iii};
    thisSpeed = speed{iii};
    sigRaster(iii,:) = thisspSign(1:minL);
    speedRaster(iii,:) = thisSpeed(1:minL);
end

out.sigRaster = sigRaster;
out.speedRaster = speedRaster;
temp = ts{1};
out.times = temp(1:minL);

