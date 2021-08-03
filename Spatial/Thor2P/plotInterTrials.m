function oSignals =  plotInterTrials(ei,signal)
b = ei.b;
onset = b.air_puff_f(b.sTrials(1:(end-1)));
offset = b.air_puff_r(b.sTrials(2:end))+100;

trialNumbers = 1:length(onset);
si = 1/ei.frameRate;
for iii = 1:length(trialNumbers)
    ii = trialNumbers(iii);
    try
        st = onset(ii);
        se = offset(ii);
    catch
        error;
    end
    frames = (find(b.frames_f >= st & b.frames_f <= se)); % find frame numbers
    oSignals{iii} = signal(frames);
    maxSignals(iii) = max(oSignals{iii});
    minSignals(iii) = min(oSignals{iii});
    ts{iii} = b.ts(b.frames_f(frames))-b.ts(st);
    speed{iii} = b.speed(st:se);
%     [speed{iii},~] = envelope(b.speed(st:se));
    tss{iii} = b.ts(st:se)-b.ts(st);
    tone_light{iii} = (find(b.tone_light_stim_r >= st & b.tone_light_stim_r <= se)); % find frame numbers
    tstl{iii} = b.ts(b.tone_light_stim_r(tone_light{iii}))-b.ts(st);
end
n=0;
mmSig = max(maxSignals);
mnSig = min(minSignals);
ff = makeFigureRowsCols(556,[7 1 3 9],'RowsCols',[length(trialNumbers) 1],'spaceRowsCols',[0.021 0.0009],'rightUpShifts',[0.05 0.05],'widthHeightAdjustment',[-70 -25]);
for ii = 1:length(trialNumbers)
    axes(ff.h_axes(ii,1));
    plot(ts{ii},oSignals{ii});hold on;
%     plot(tss{ii},speed{ii});
    plot(tstl{ii},ones(size(tstl{ii})) * mmSig/2,'.');
    ylim([mnSig mmSig]);
    box off;
    ylabel(ii);
    set(ff.h_axes(ii,1),'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
end

% ff = makeFigureRowsCols(556,[5 1 4 9],'RowsCols',[length(trialNumbers) 1],'spaceRowsCols',[0.021 0.0009],'rightUpShifts',[0.05 0.05],'widthHeightAdjustment',[-70 -25]);
% for ii = 1:length(trialNumbers)
%     axes(ff.h_axes(ii,1));
%     plot(tss{ii},speed{ii});
%     box off;
%     ylabel(ii);
%     set(ff.h_axes(ii,1),'XTick',[],'XTickLabel',[],'YTick',[],'YTickLabel',[]);
% end

