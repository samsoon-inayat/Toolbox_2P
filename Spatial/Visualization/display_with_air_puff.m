function diplay_with_air_puff(b,st,se)
figure(1000);clf;hold on;
plot(b.ts/60,b.air_puff_raw,'k','linewidth',0.5);hold on;
plot(b.ts/60,0.35*(b.fSpeed/max(b.fSpeed)),'linewidth',0.75,'color','c');
zm = nan(size(b.air_puff_raw)); zm(st-1) = 0; zm(st) = 1;
plot(b.ts/60,zm,'b-','linewidth',0.5);
zm = nan(size(b.air_puff_raw)); zm(se-1) = 0; zm(se) = 1;
plot(b.ts/60,zm,'m-','linewidth',0.5);
