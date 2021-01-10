function [ptco,CRc,cellNums] = findPopulationVectorPlot(ptc,ccs,cellNums)


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
% ptc = fillmissing(ptc,'linear',2,'EndValues','nearest');
% if sum(isnan(ptc),'all')> 0
%     n = 0;
% end
% GF = gausswin(3);
% ptc = filter(GF,1,ptc,[],2);
% [ptco,window] = smoothdata(ptc,2,'gaussian',[2 2]);
ptco = ptc;
CRc = corrcoef(ptco);