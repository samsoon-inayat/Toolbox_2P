function Rs = findRasterProperties (thispFolder,contextNumber,stimMarker,Rs,thisRasterType,trials,ow)

if sum(ow) == -3
    return;
end

allContexts = contextDefinitions;
% 
% allContexts.rasterFunctions = {'find_MI';'getFits_myGaussFit';'fractalDim'};
% allContexts.varNames = {'info_metrics';'gauss_fit_on_mean';'fractal_dim'};


Rs_functions = allContexts.rasterFunctions(1:3);
Rs_varNames = allContexts.varNames(1:3);

for ii = 1:length(Rs_functions)
    if ow(ii) == -1
        continue;
    end
    functionName = Rs_functions{ii};
    varName = Rs_varNames{ii};
    cmdTxt = sprintf('Rs.%s = %s(thispFolder,contextNumber,stimMarker,Rs,thisRasterType,trials,ow(ii));',varName,functionName);
    eval(cmdTxt);
end