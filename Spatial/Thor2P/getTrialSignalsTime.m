function [thisTs oSignals varargout] = getTrialSignalsTime (signal,b,onset,offset)
trialNumbers = 1:(length(b.sTrials)-1);
for iii = 1:length(trialNumbers)
    ii = trialNumbers(iii);
    st = onset(ii);
    se = offset(ii);
    frames = (find(b.frames_f >= st & b.frames_f <= se)); % find frame number
    if isempty(frames)
        iii
        error;
        break;
    end
    thisSignal = signal(frames);
%     figure(1);clf;
%     plot(thisSignal);
    numSamplesSignal(ii) = length(thisSignal);
end
numSamples = min(numSamplesSignal);
for iii = 1:length(trialNumbers)
    ii = trialNumbers(iii);
    st = onset(ii);
    se = offset(ii);
    frame_times = b.frames_f(find(b.frames_f >= st & b.frames_f <= se));
    frames = (find(b.frames_f >= st & b.frames_f <= se)); % find frame number
%     if isempty(frames)
%         break;
%     end
    thisSpeed1 = b.speed(frame_times);
    thisTs = (b.ts(frame_times)-b.ts(st));
    thisSignal1 = signal(frames);
    xs = 1:length(thisSignal1);
    thisSignal = interp1(xs,thisSignal1,1:numSamples);
    xs = 1:length(thisSpeed1);
    thisSpeed = interp1(xs,thisSpeed1,1:numSamples);
    oSignals(iii,:) = thisSignal-min(thisSignal);
    speed(iii,:) = thisSpeed;
end
thisTs = interp1(xs,thisTs,1:numSamples);
if nargout > 2
    varargout{1} = speed;
end