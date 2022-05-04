function Rs = findRasterProperties_1 (thispFolder,contextNumber,stimMarker,Rsi,thisRasterType,trials,ow)

if sum(ow) == -3
    return;
end

allContexts = contextDefinitions;
% 
% allContexts.rasterFunctions = {'find_MI';'getFits_myGaussFit';'fractalDim'};
% allContexts.varNames = {'info_metrics';'gauss_fit_on_mean';'fractal_dim'};

Rs = Rsi{1};
Rs_functions = {'find_MI_1_MC'};
Rs_varNames = {'info_metrics_MC'};

for ii = 1:length(Rs_functions)
    if ow(ii) == -1
        continue;
    end
    functionName = Rs_functions{ii};
    varName = Rs_varNames{ii};
    cmdTxt = sprintf('Rs.%s = %s(thispFolder,contextNumber,stimMarker,Rsi,thisRasterType,trials,ow(ii));',varName,functionName);
    eval(cmdTxt);
end