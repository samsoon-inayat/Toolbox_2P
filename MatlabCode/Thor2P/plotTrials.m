function oSignals =  plotTrials(b,signal,spSignal,onset,offset)

figure(100);clf;
plot(b.ts,b.air_puff_raw);hold on;
plot(b.ts(onset),0.5*ones(size(onset)),'o');

trialNumbers = 1:length(onset);
% si = 1/ei.frameRate;
for iii = 1:length(trialNumbers)
    if iii == 15
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
    oSignals{iii} = signal(frames);
%     spoSignals{iii} = signal(frames);
    spoSignals{iii} = spSignal(frames);
    maxSignals(iii) = max(oSignals{iii});
    minSignals(iii) = min(oSignals{iii});
    ts{iii} = b.ts(b.frames_f(frames))-b.ts(st);
    temp = b.speed(st:se);
    temp(temp<0) = 0;
    speed{iii} = temp;
%     [speed{iii},~] = envelope(b.speed(st:se));
    tss{iii} = b.ts(st:se)-b.ts(st);
    dist{iii} = b.dist(st:se)-b.dist(st);
    distVal(iii) = dist{iii}(end);
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
%     
%     tempTL = (find(b.tone_light_stim_r >= st & b.tone_light_stim_r <= se)); % find frame numbers
%     tone_light{iii} = tempTL(1:2);
%     tstl{iii} = b.ts(b.tone_light_stim_r(tone_light{iii}))-b.ts(st);
%     dstl(iii,:) = b.dist(b.tone_light_stim_r(tone_light{iii}))-b.dist(st);
end
minDist = min(distVal);
nbins = 50;
bins = linspace(0,minDist,(nbins+1));
sbins = bins(1:nbins);
ebins = bins(2:end);
GF = gausswin(3);
for iii = 1:length(trialNumbers)
    thisDist = dist{iii};
    thisT = tss{iii};
    thisTSig = ts{iii};
    thisspSign = spoSignals{iii};
    for jj = 1:nbins
        st = sbins(jj);  se = ebins(jj);
        dVals = find(thisDist >= st & thisDist <se);
        tVals = thisT(dVals);
        stt = tVals(1);
        sett = tVals(end);
        tValsInd = find(thisTSig >= stt & thisTSig < sett);
        dSig(iii,jj) = mean(thisspSign(tValsInd))/(sett-stt);
        dTim(iii,jj) = sett-stt;
    end
    dSigT = [0 0 dSig(iii,:)];
    dSigT = filter(GF,1,dSigT);
    dSigF(iii,:) = dSigT(3:end);
end

n=0;
mmSig = max(maxSignals);
mnSig = min(minSignals);
ff = makeFigureRowsCols(555,[1 1 6 9],'RowsCols',[length(trialNumbers) 1],'spaceRowsCols',[0.021 0.0009],'rightUpShifts',[0.05 0.05],'widthHeightAdjustment',[-70 -25]);
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

mdstl = mean(dstl);
toneDist = mdstl(1);
lightDist = mdstl(2);
binI = bins-toneDist;
tBinI = find(binI>0,1,'first');
binI = bins-lightDist;
lBinI = find(binI>0,1,'first');

figure(666);clf;
imagesc(dSigF);
colorbar;
xVs = get(gca,'XTick');
set(gca,'XTickLabel',round(bins(xVs+1)));
text(tBinI,1,'T','color','r');
text(lBinI,1,'L','color','r');
figure(667);clf;
mdSig = mean(dSig);
dSigT = [0 0 mdSig];
dSigT = filter(GF,1,dSigT);
plot(dSigT(3:end));


% ff = makeFigureRowsCols(556,[5 1 4 9],'RowsCols',[length(trialNumbers) 1],'spaceRowsCols',[0.021 0.0009],'rightUpShifts',[0.05 0.05],'widthHeightAdjustment',[-70 -25]);
% for ii = 1:length(trialNumbers)
%     axes(ff.h_axes(ii,1));
%     plot(tss{ii},speed{ii});
%     box off;
%     ylabel(ii);
%     set(ff.h_axes(ii,1),'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
% end

