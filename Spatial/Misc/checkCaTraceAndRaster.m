function checkContextTrials(ei,plane,contextNumber)

context = ei{1}.plane{plane}.contexts(contextNumber);
cellsSigs = ei{1}.plane{plane}.tP.deconv.spSigAll;

hf = figure(1000);clf;
hf1 = figure(1001);clf;
for ii = 1:length(cellsSigs)
    figure(hf1);clf;
    thisSig = cellsSigs{ii};
    plot(thisSig);
    pause;
end

