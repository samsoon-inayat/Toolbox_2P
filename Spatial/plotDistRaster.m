function plotDistRaster(ha,A,cn,spca,minMax,titleText,meanPlot)

axes(ha);

if spca
    A_raster = A.caRasters(:,:,cn);
else
    A_raster = A.rasters(:,:,cn);
end

dist = A.dist;

if ~exist('minMax','var')
    minMax = [min(A_raster(:)) max(A_raster(:))];
end

if isempty(minMax)
    minMax = [min(A_raster(:)) max(A_raster(:))];
end

mSig = nanmean(A_raster);
% ft = fittype(A.formula);
% xs = 1:size(A.rasters,2);
% coeff = A.coeff(cn,:);
% fitplot = ft(coeff(1),coeff(2),coeff(3),xs);
mDR = minMax(1);
MDR = minMax(2);

if ~meanPlot
    plot(A.dist,nanmean(A_raster));
else
    imagesc(A_raster,[mDR MDR]);colorbar;
    xdata = get(gca,'XTick');
    set(gca,'XTickLabel',floor(dist(xdata)),'YDir','normal');
    meanResp = nanmean(A_raster);
    hold on;
    plot(size(A_raster,1)*meanResp/max(meanResp),'linewidth',1,'color','w');
%     plot(size(A_raster,1)*fitplot/max(fitplot),'linewidth',1.5,'color','r');
end

%     titleText = sprintf('Cell %d -- Air FR',ccs(cc),SI,pSIs,PCSI,PC);
if ~exist('titleText','var')
    titleText = [];
end
title(titleText);
