function corrs = findPearsonCorrelation(R1,R2)
ccs = 1:size(R1.rasters,3);

for cn = 1:length(ccs)
    meanR1 = applyGaussFilt(nanmean(R1.rasters(:,:,cn)),5);
    meanR2 = applyGaussFilt(nanmean(R2.rasters(:,:,cn)),5);
    [RHO,PVAL] = corrcoef(meanR1',meanR2');
    corrs(cn) = RHO(1,2);
end

