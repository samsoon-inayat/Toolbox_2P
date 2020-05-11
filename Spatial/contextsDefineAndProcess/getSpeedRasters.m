function rasters =  getSpeedRasters(ei,onsets,offsets,overwrite,plotFlag)

fileName = makeName(sprintf('sRaster_%d_%d_%d_%d.mat',onsets(1),onsets(end),offsets(1),offsets(end)),ei.plane{1}.folder);
if exist(fileName,'file') && overwrite == 0
    load(fileName);
    return;
end

b = ei.b;
out = getDistRaster_forSpeed(b,onsets,offsets);
rasters.speed = out.distSpeedRaster;
rasters.duration = out.distDurRaster;
rasters.dist = out.dist;
rasters.onsets = onsets;
rasters.offsets = offsets;
save(fileName,'rasters','-v7.3');

function out =  getDistRaster_forSpeed(b,onsets,offsets)
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
    distVal(iii) = dist{iii}(end);
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
nbins = 50;
bins = linspace(0,minDist,(nbins+1));
sbins = bins(1:nbins);
ebins = bins(2:end);

%%% Method %%%
% 1) Divide the minDist into nbins

for iii = 1:length(trialNumbers)
    thisDist = dist{iii};
    thisT = tss{iii};
%     thisTSig = ts{iii};
%     thisspSign = spoSignals{iii};
    thisSpeed = speed{iii};
    for jj = 1:nbins
        st = sbins(jj);  se = ebins(jj);
        dVals = find(thisDist >= st & thisDist <se);
        tVals = thisT(dVals);
        stt = tVals(1);
        sett = tVals(end);
%         tValsInd = find(thisTSig >= stt & thisTSig < sett);
%         if isnan(mean(thisspSign(tValsInd))/(sett-stt))
%             error('NaN in raster calculation');
%         end
%         distSigRaster(iii,jj) = mean(thisspSign(tValsInd));
%         distSpikesRaster(iii,jj) = sum(thisspSign(tValsInd))*(sett-stt);
%         distSigRaster(iii,jj) = mean(thisspSign(tValsInd))/(sett-stt);
        distDurRaster(iii,jj) = sett-stt;
        distSpeedRaster(iii,jj) = mean(thisSpeed(dVals));
    end
%     dSigF(iii,:) = applyGaussFilt(distSigRaster(iii,:),3);
end

out.distSpeedRaster = distSpeedRaster;
out.distDurRaster = distDurRaster;
out.dist = sbins;








