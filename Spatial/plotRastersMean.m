function plotRastersMean(R1,R2,ccs,spca,pmean)
% ccs = 1:size(R1.rasters,3);

numberOfRows = 2;
numberOfCols = 2;
graphsOnOneFigure = numberOfRows * numberOfCols;
numberOfData = length(ccs);
numberOfGroups = ceil(numberOfData/graphsOnOneFigure);
numberOfFullGroups = floor(numberOfData/graphsOnOneFigure);
indices = NaN(1,(numberOfGroups*graphsOnOneFigure));
indices(1:numberOfData) = 1:numberOfData;
groupIndices = reshape(indices,numberOfRows,numberOfCols,numberOfGroups);
for gg = 1:numberOfGroups
    groupIndices(:,:,gg) = groupIndices(:,:,gg)';
end

ff = makeFigureRowsCols(105,[NaN 0.1 7 5],'RowsCols',[numberOfRows numberOfCols],...
    'spaceRowsCols',[0.1 0.05],'rightUpShifts',[0.05 0.07],'widthHeightAdjustment',...
    [-70 -115]);

gg = 1;
while 1
    for rr = 1:numberOfRows
        for cc = 1:numberOfCols
        cni = groupIndices(rr,cc,gg);
        cn = ccs(cni);
%         cng(rr) = cn;
        if isnan(cn)
            continue;
        end
        if spca
            A_P_raster = R1.caRasters(:,:,cn);
            AP_P_raster = R2.caRasters(:,:,cn);
        else
            A_P_raster = R1.rasters(:,:,cn);
            AP_P_raster = R2.rasters(:,:,cn);
        end
        minD = min([A_P_raster(:);AP_P_raster(:)]);
        maxD = max([A_P_raster(:);AP_P_raster(:)]);
        meanR1 = applyGaussFilt(nanmean(R1.rasters(:,:,cn)),5);
        meanR2 = applyGaussFilt(nanmean(R2.rasters(:,:,cn)),5);
        [RHO,PVAL] = corrcoef(meanR1',meanR2');
        ylims = [min([meanR1 meanR2]) max([meanR1 meanR2])];
        
        
            axes(ff.h_axes(rr,cc));
%             if cc == 1
%                 if pmean%mod(rr,2) == 0
                    thisRaster = R1.rasters(:,:,cn);
                    axes(ff.h_axes(rr,cc));cla
                    plot(R1.dist,meanR1,'color','k','linewidth',1.5);
                    xlim([min(R1.dist) max(R1.dist)]);
                    ylim(ylims);
                    hold on;
                    plot(R1.dist,meanR2,'color','m','linewidth',1.5);
                    xlim([min(R1.dist) max(R1.dist)]);
%                     title(sprintf('%.3f , %.3f',RHO(1,2),PVAL(1,2)));
                    title(sprintf('%.3f',RHO(1,2)));
%                 else
%                     plotDistRaster(ff.h_axes(rr,cc),R1,cn,[minD maxD],sprintf('%s - %.2f',num2str(ccs(cn)),R1.SI(cn)));
%                 end
%             end
%             if cc == 2
%                 if pmean%mod(rr,2) == 0
%                     thisRaster = R2.rasters(:,:,cn);
%                     axes(ff.h_axes(rr,cc));
%                     plot(meanR2);
%                     ylim(ylims);
%                 else
%                     plotDistRaster(ff.h_axes(rr,cc),R2,cn,[minD maxD],sprintf('%s - %.2f',num2str(ccs(cn)),R2.SI(cn)));
%                 end
%             end
%             if cc == 2
%                 plot(meanR2-meanR1);
%             end
        end
    end
%     showCells(1000,ei,ccs(cng));
%     save2pdf('temp.pdf',gcf,600);
    display('Press key');
    gg = keyboardInput(gg,[1 numberOfGroups],[1 6],'');
    if gg < 0
        break;
    end
end
