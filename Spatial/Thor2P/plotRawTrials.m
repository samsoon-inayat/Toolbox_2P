function plotRawTrials (xvals,signals,separation,figNum)
xs = [];
signal = [];
stim = [];
for ii = 1:size(signals,1)
    if ii == 1
        xs = [xvals NaN];
    else
        xs = [xs xvals+max(xs)+separation NaN];
    end
    signal = [signal signals(ii,:) NaN];
    stim = [stim 0 ones(1,size(signals,2)-2) 0 NaN];
end
figure(figNum);clf;
plot(xs,signal);hold on;
plot(xs,stim * max(signal)/2,'r');
end