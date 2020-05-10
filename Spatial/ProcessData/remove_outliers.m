function inds = remove_outliers(time_series,sensitivityFactor)
absoluteDeviation = abs(time_series - mean(time_series));
mad = median(absoluteDeviation);
% sensitivityFactor = 50; % Whatever you want.
thresholdValue = sensitivityFactor * mad;
outlierIndexes = abs(absoluteDeviation) > thresholdValue;
inds = find(outlierIndexes);
% time_series(inds) = [];
% time_series(inds) = NaN; time_series = fillmissing(time_series,'spline');