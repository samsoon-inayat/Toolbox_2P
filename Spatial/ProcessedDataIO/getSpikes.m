function [caSigAll,spSigAll] = getSpikes(ei,tP,ow)
ow = 0;
fileName = makeName('deconvolved.mat',ei.folder);
if exist(fileName,'file') & ~ow
    load(fileName);
    return;
end

signals = tP.signals;
% ccs = find(tP.areCells);
ccs = 1:size(signals,1);
caSigAll = cell(1,length(ccs));
spSigAll = caSigAll;
for cc = 1:length(ccs)
    disp(cc);
    caSignal = signals(ccs(cc),:);
    try
        [ca_sig,tsp,~] = deconvolveCa(caSignal, 'ar1',...
        'constrained','optimize_b', 'optimize_pars', 'optimize_smin');
    catch
        disp('Error occured while deconvolving');
        ca_sig = nan(size(caSignal));
        tsp = ca_sig;
    end
    if length(tsp) < length(caSignal)
        len_diff = length(caSignal) - length(tsp);
        tsp((length(tsp)+1):length(caSignal),1) = zeros(len_diff,1);
    end
    caSigAll{cc} = ca_sig;
    spSigAll{cc} = tsp;    
end

save(fileName,'caSigAll','spSigAll');