function Rs = findRasterProperties (thispFolder,contextNumber,stimMarker,Rs,thisRasterType,trials,ow)

if sum(ow) == -3
    return;
end

allContexts = contextDefinitions;


Rs_functions = allContexts.rasterFunctions;
Rs_varNames = allContexts.varNames;

for ii = 1:length(Rs_functions)
    if ow(ii) == -1
        continue;
    end
    functionName = Rs_functions{ii};
    varName = Rs_varNames{ii};
    cmdTxt = sprintf('Rs.%s = %s(thispFolder,contextNumber,stimMarker,Rs,thisRasterType,trials,ow(ii));',varName,functionName);
    eval(cmdTxt);
end