function [caSigAll,spSigAll] = getSpikes(ei,tP,ow)
% ow = 0;
fileName = makeName('deconvolved.mat',ei.folder);
if exist(fileName,'file') & ~ow
    load(fileName);
    return;
end

signals = tP.signals;
ccs = find(tP.iscell(:,1));
% ccs = 1:size(signals,1);
% caSigAll = cell(1,length(ccs));
caSigAll = zeros(length(ccs),size(signals,2));
spSigAll = caSigAll;
parfor cc = 1:length(ccs)
%     disp(cc);
%     caSignal = signals(ccs(cc),:);
    try
        [caSigAll(cc,:),spSigAll(cc,:),~] = deconvolveCa(signals(ccs(cc),:), 'ar1',...
        'constrained','optimize_b', 'optimize_pars', 'optimize_smin');
    catch
%         disp('Error occured while deconvolving');
%         ca_sig = nan(size(caSignal));
%         tsp = ca_sig;
    end
%     if length(tsp) < length(caSignal)
%         len_diff = length(caSignal) - length(tsp);
%         tsp((length(tsp)+1):length(caSignal),1) = zeros(len_diff,1);
%     end
%      = ca_sig;
%      = tsp;    
end

save(fileName,'caSigAll','spSigAll');