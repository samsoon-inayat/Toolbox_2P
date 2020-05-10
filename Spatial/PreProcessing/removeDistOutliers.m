function Disti = removeDistOutliers(Disti)
Dist = diff(Disti);
absoluteDeviation = abs(Dist - mean(Dist));
mad = median(absoluteDeviation);
sensitivityFactor = 50; % Whatever you want.
thresholdValue = sensitivityFactor * mad;
outlierIndexes = abs(absoluteDeviation) > thresholdValue;
inds = find(outlierIndexes)+1;
Disti(inds) = NaN; Disti = fillmissing(Disti,'spline');