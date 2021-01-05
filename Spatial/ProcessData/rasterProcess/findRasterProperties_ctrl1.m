function Rs = findRasterProperties_ctrl1 (thispFolder,contextNumber,stimMarker,Rs,thisRasterType,trials,ow)

if sum(ow) == -3
    return;
end

allContexts = contextDefinitions;
% 
% allContexts.rasterFunctions = {'find_MI';'getFits_myGaussFit';'fractalDim'};
% allContexts.varNames = {'info_metrics';'gauss_fit_on_mean';'fractal_dim'};


Rs_functions = allContexts.rasterFunctions(7:9);
Rs_varNames = allContexts.varNames(7:9);

for ii = 1:length(Rs_functions)
    if ow(ii) == -1
        continue;
    end
    functionName = Rs_functions{ii};
    varName = Rs_varNames{ii};
    varNameV = varName(1:(length(varName)-2));
    cmdTxt = sprintf('Rs.%s = %s(thispFolder,contextNumber,stimMarker,Rs,thisRasterType,trials,ow(ii));',varNameV,functionName);
    eval(cmdTxt);
end