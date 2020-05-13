function out = alignPeaks(signals)

for ii = 1:length(signals)
    thisSig = signals{ii};
    maxPos(ii) = find(thisSig == max(thisSig),1,'first');
    lenSig(ii) = length(thisSig);
end
centerFrame = min(maxPos);
framesBefore = centerFrame - 1;
framesAfter = min(lenSig - maxPos);
for ii = 1:length(signals)
    thisSig = signals{ii};
    alignedSignals(ii,:) = thisSig((maxPos(ii)-framesBefore):(maxPos(ii)+framesAfter));
    indices(ii,:) = (maxPos(ii)-framesBefore):(maxPos(ii)+framesAfter);
end
out.alignedSignals = alignedSignals;
out.indices = indices;
out.maxFrame = centerFrame;
out.framesBefore = framesBefore;
out.framesAfter = framesAfter;
