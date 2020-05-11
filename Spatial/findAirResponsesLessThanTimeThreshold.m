function cellList = findAirResponsesLessThanTimeThreshold(Rs,ccs,st,et)
ptc = findMeanRasters(Rs);
ptc = ptc(ccs,:);
ptc = normalizeSignal(ptc,2);
[~,peakPos] = max(ptc,[],2);
cols = size(Rs.rasters(:,:,1),2);
colsHalf = ceil(cols/2);
ts = Rs.cells(1).times - Rs.cells(1).times(colsHalf);
peakVals = find(ts>=st & ts<=et);
cellList = find(peakPos>=peakVals(1)&peakPos<=peakVals(end));
cellList = ccs(cellList);
