function plotMotion(ei,motionOnsets,motionOffsets)
figure(101);clf
plot(ei.b.ts,ei.b.fSpeed);hold on;
tMOn = zeros(size(ei.b.fSpeed));
tMOn(motionOnsets) = 15;
plot(ei.b.ts,tMOn,'m');
tMOn = zeros(size(ei.b.fSpeed));
tMOn(motionOffsets) = 15;
plot(ei.b.ts,tMOn,'c');
plot(ei.b.ts,ei.b.air_puff_raw,'k');
