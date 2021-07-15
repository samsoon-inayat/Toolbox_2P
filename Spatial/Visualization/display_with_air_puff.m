function diplay_with_air_puff(b,st,se)
figure(1000);clf;
plot(b.ts,b.air_puff_raw);hold on;
plot(b.ts,0.25*(b.fSpeed/max(b.fSpeed)));
zm = zeros(size(b.air_puff_raw)); zm(st) = 1;
plot(b.ts,zm);
zm = zeros(size(b.air_puff_raw)); zm(se) = 1;
plot(b.ts,zm);
