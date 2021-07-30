function [avg,allVals] = getAveragesAndAllValues(distD)

avg = arrayfun(@(x) nanmean(x{1}),distD);

allVals = [];
for rr = 1:size(distD,1)
    for cc = 1:size(distD,2)
        allVals = [allVals;distD{rr,cc}(:)];
    end
end
