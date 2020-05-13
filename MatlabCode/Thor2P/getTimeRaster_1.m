function out =  getTimeRaster_1(b,caSig,spSignal,onsets,offsets,plotFlag)
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
    if isfield(b,'tone_light_stim')
        tempTL = (find(b.tone_light_stim_r >= st & b.tone_light_stim_r <= se)); % find frame numbers
        try
            tone_light{iii} = tempTL(3:4);
        catch
            tone_light{iii} = tempTL(1:2);
        end
        tstl{iii} = b.ts(b.tone_light_stim_r(tone_light{iii}))-b.ts(st);
        dstl(iii,:) = b.dist(b.tone_light_stim_r(tone_light{iii}))-b.dist(st);
    end
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

if plotFlag
n=0;
mmSig = max(maxSignals);
mnSig = min(minSignals);
ff = makeFigureRowsCols(556,[6.2 -0.5 6 9],'RowsCols',[length(trialNumbers) 1],'spaceRowsCols',[0.021 0.0009],'rightUpShifts',[0.05 0.05],'widthHeightAdjustment',[-70 -25]);
for ii = 1:length(trialNumbers)
    axes(ff.h_axes(ii,1));
    plot(ts{ii},oSignals{ii});hold on;
    plot(ts{ii},spoSignals{ii});
    tempSpeed = speed{ii}*0.45;
    plot(tss{ii},tempSpeed);
    plot(tstl{ii},ones(size(tstl{ii})) * mmSig/2,'.');
    ylim([mnSig mmSig]);
    box off;
    ylabel(ii);
    set(ff.h_axes(ii,1),'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
end

% mdstl = mean(dstl);
% toneDist = mdstl(1);
% lightDist = mdstl(2);
% binI = bins-toneDist;
% tBinI = find(binI>0,1,'first');
% binI = bins-lightDist;
% lBinI = find(binI>0,1,'first');

% figure(666);clf;
% imagesc(sigRaster);
% colorbar;
% xVs = get(gca,'XTick');
% set(gca,'XTickLabel',round(bins(xVs+1)));
% text(tBinI,1,'T','color','r');
% text(lBinI,1,'L','color','r');
% figure(667);clf;
% mdSig = nanmean(sigRaster);
% dSigT = applyGaussFilt(mdSig,3);
% plot(dSigT(3:end));
end