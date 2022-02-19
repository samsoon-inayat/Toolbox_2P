function trying

ei = evalin('base','ei');
n = 0;

for ii = 1:length(ei)
    b = ei{ii}.b;
    ston = b.stim_r;
    stoff = b.stim_f;
    
    dT(ii) = median(b.ts(stoff) - b.ts(ston));
    adT = (b.ts(ston(2:80)) - b.ts(stoff(1:79)));
    
    ston = b.air_puff_r;
    stoff = b.air_puff_f;
    
    dTA(ii) = median(b.ts(stoff) - b.ts(ston));
    FR(ii) = ei{ii}.thorExp.frameRate;
end
