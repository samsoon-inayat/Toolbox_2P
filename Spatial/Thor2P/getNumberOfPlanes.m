function nPlanes = getNumberOfPlanes(folderName)

if isstruct(folderName)
    ei = folderName;
else
    ei = thorGetExperimentInfo(folderName);
end
if ei.zFastEnable == 1
    nPlanes = ei.zSteps;
    return;
else
    nPlanes = 1;
end