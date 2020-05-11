function [ptc,CRc,cellNums] = findPopulationVectorPlot(ptc,ccs,cellNums)


% if ~isempty(trials)
%     ptc = findMeanRasters(Rs,trials);
% else
%     ptc = findMeanRasters(Rs);
% end

if isempty(ccs)
    ccs = 1:size(ptc,1);
end

ptc = ptc(ccs,:);
ptc = normalizeSignal(ptc,2);
if ~exist('cellNums','var')
    [~,peakPos] = max(ptc,[],2);
    [~,cellNums] = sort(peakPos);
    ptc = ptc(cellNums,:);
else
    ptc = ptc(cellNums,:);
end
CRc = corrcoef(ptc);