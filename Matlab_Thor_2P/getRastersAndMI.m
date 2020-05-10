%% finding rasters and mutual information
function ei = getRastersAndMI(ei,trials,owr,owmi)
    R1 = getRasters(ei,trials,owr);
    temp = find_mutual_information(ei,R1,owmi);
    R1.zMI = temp.zMI;
    R1.MI = temp.MI;
    ei.rasters = R1;
end