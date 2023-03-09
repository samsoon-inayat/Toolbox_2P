function [pva,cellNums,peakPos] = align_by_peaks(ptc,cellNums)

if nargin == 1
    [~,peakPos] = max(ptc,[],2);
    [~,cellNums] = sort(peakPos);
    pva = ptc(cellNums,:);
    peakPos = peakPos(cellNums);
    return;
end

if nargin == 2
    for ii = 1:length(ptc)
        tptc = ptc{ii};
        pva{ii} = tptc(cellNums,:);
    end
    peakPos = [];
end