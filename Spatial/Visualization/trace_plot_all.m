function trace_plot_all(pd_rec)

signals = pd_rec.deconv';
frames_f = pd_rec.b.frames_f;
b = pd_rec.b;

traceTime = b.ts(frames_f);

% onsets = b.photo_sensor_f;
% offsets = b.photo_sensor_r;
% st = find(traceTime>(b.ts(onsets(1))-5),1,'first'); en = find(traceTime<(b.ts(offsets(end))+5),1','last');
% xStart = traceTime(st);
% xEnd = traceTime(en);
% xlims = [xStart xEnd];


%%
figure(100);clf;
spSigAllN = normalizeSignal(signals,2);
[maxVal,maxLoc] = max(spSigAllN,[],2);
[sorted,locs] = sort(maxLoc);
hi = imagesc([traceTime(1), traceTime(end)],[1 size(spSigAllN,1)],spSigAllN(locs,:),[0 0.75]);colorbar;
hold on;
% plot(b.ts(frames_f),spSigAll(:,:));
ylims = ylim;
hold on;
if isfield(b,'air_puff_raw')
   plot(b.ts,b.air_puff_raw*ylims(2),'g','linewidth',0.5)
end
if isfield(b,'stim_raw')
   plot(b.ts,b.stim_raw*ylims(2),'r','linewidth',0.5)
end
% plot(b.ts,b.photo_sensor_raw*ylims(2),'c','linewidth',0.5)
plot(b.ts,b.fSpeed,'y');
set(gca,'Ydir','Normal');
xlabel('Time (sec)');
ylabel('Cell Number (arranged by peaks)');
set(gca,'FontSize',12,'FontWeight','Bold');
title('Deconvolved Signal of All Cells');