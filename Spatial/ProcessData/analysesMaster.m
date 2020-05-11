function analysesMaster
%%
aei = evalin('base','ei');

%% Cellular responses ... distance and temporal
dist_time_response_analyses

%% Correlation plots
for ii = 1:length(aei)
    ei = aei(1);
    cellCorrelation(ei,10+ii);
end

%% Speed Plots
tei = ei([1 2]);
getSpeedResponsePlot(tei);
