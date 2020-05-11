function speed = removeSpeedOutliers(speed)
absoluteDeviation = abs(speed - mean(speed));
mad = median(absoluteDeviation);
sensitivityFactor = 50; % Whatever you want.
thresholdValue = sensitivityFactor * mad;
outlierIndexes = abs(absoluteDeviation) > thresholdValue;
inds = find(outlierIndexes);
speed(inds) = NaN; speed = fillmissing(speed,'spline');