function plotTrialsForSelection(b)
if ~exist('b')
ei = evalin('base','ei{7}');
b = ei.b;
end

figure(2000);clf;
frames = zeros(size(b.ts));
frames(b.frames_f) = 0.5;
subplot 211;
plot(b.ts,frames,'linewidth',0.25);hold on;
plot(b.ts,b.air_puff_raw,'linewidth',4);
for ii = 1:length(b.air_puff_f)
    text(b.ts(b.air_puff_f(ii))+5,0.75,num2str(ii));
end
subplot 212;
frames = zeros(size(b.ts));
frames(b.ch_a_r) = 0.5;
plot(b.ts,frames,'linewidth',0.25);hold on;
plot(b.ts,b.air_puff_raw,'linewidth',4);