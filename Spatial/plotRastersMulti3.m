function plotRastersMulti3(R1,R2,R3,R4,ccs,spca,pmean)
% ccs = 1:size(R1.rasters,3);

numberOfRows = 4;
numberOfCols = 1;
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

ff = makeFigureRowsCols(105,[NaN 0.1 7 5],'RowsCols',[numberOfRows 3],...
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
                AP_PP_raster = R3.caRasters(:,:,cn);
                AP_PP_rasterA = R4.caRasters(:,:,cn);
            else
                A_P_raster = R1.rasters(:,:,cn);
                AP_P_raster = R2.rasters(:,:,cn);
                AP_PP_raster = R3.rasters(:,:,cn);
                AP_PP_rasterA = R4.rasters(:,:,cn);
            end
            minD = min([A_P_raster(:);AP_P_raster(:);AP_PP_raster(:);AP_PP_rasterA(:)]);
            maxD = max([A_P_raster(:);AP_P_raster(:);AP_PP_raster(:);AP_PP_rasterA(:)]);
            plotDistRaster(ff.h_axes(rr,1),R1,cn,0,[minD maxD],sprintf('%s - %.2f',num2str(ccs(cni)),R1.SI(cn)),pmean);
            plotDistRaster(ff.h_axes(rr,2),R2,cn,0,[minD maxD],sprintf('%s - %.2f',num2str(ccs(cni)),R2.SI(cn)),pmean);
            plotDistRaster(ff.h_axes(rr,3),R3,cn,0,[minD maxD],sprintf('%s - %.2f',num2str(ccs(cni)),R3.SI(cn)),pmean);
%             plotDistRaster(ff.h_axes(rr,4),R4,cn,0,[minD maxD],sprintf('%s - %.2f',num2str(ccs(cni)),R4.SI(cn)),pmean);

    %         meanR1 = applyGaussFilt(nanmean(R1.rasters(:,:,cn)),5);
    %         meanR2 = applyGaussFilt(nanmean(R2.rasters(:,:,cn)),5);
    %         [RHO,PVAL] = corrcoef(meanR1',meanR2');
    %         ylims = [min([meanR1 meanR2]) max([meanR1 meanR2])];
%             axes(ff.h_axes(rr,cc));
%             if cc == 1
%                 plotDistRaster(ff.h_axes(rr,cc),R1,cn,0,[minD maxD],sprintf('%s - %.2f',num2str(ccs(cni)),R1.SI(cn)),pmean);
%             end
%             if cc == 2
%                 plotDistRaster(ff.h_axes(rr,cc),R2,cn,0,[minD maxD],sprintf('%s - %.2f',num2str(ccs(cni)),R2.SI(cn)),pmean);
%             end
%             if cc == 3
%                 plotDistRaster(ff.h_axes(rr,cc),R3,cn,0,[minD maxD],sprintf('%s - %.2f',num2str(ccs(cni)),R3.SI(cn)),pmean);
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
