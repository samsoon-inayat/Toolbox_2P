function place_field_analysis_2_for_171179_Feb_13_2018_Recording_3

ei = evalin('base','ei{5}');
b = ei.b;

ccs = ei.areCells;

allA_P = getRasters(ei,1:11);
allAP_P = getRasters(ei,13:23);


plotMarkers(ei,ei.b.air_puff_r(1:23),[],101);

numberOfRows = 6;
numberOfCols = 2;
graphsOnOneFigure = numberOfRows * 1;
numberOfData = length(ccs);
numberOfGroups = ceil(numberOfData/graphsOnOneFigure);
numberOfFullGroups = floor(numberOfData/graphsOnOneFigure);
indices = NaN(1,(numberOfGroups*graphsOnOneFigure));
indices(1:numberOfData) = 1:numberOfData;
groupIndices = reshape(indices,numberOfRows,1,numberOfGroups);
for gg = 1:numberOfGroups
    groupIndices(:,:,gg) = groupIndices(:,:,gg)';
end

ff = makeFigureRowsCols(105,[NaN 0.1 11 9],'RowsCols',[numberOfRows numberOfCols],...
    'spaceRowsCols',[0.05 0.045],'rightUpShifts',[0.03 0.03],'widthHeightAdjustment',...
    [-40 -55]);

gg = 1;
while 1
    for rr = 1:numberOfRows
        cn = groupIndices(rr,1,gg);
        cng(rr) = cn;
        if isnan(cn)
            continue;
        end
        A_P = allA_P(cn); A_P_raster = A_P.raster;
        A_P_SI = A_P.SI;
        AP_P = allAP_P(cn); AP_P_raster = AP_P.raster;
        AP_P_SI = AP_P.SI;
%         P = allP(cn); P_raster = P.raster;
%         P_SI = P.SI;
        minD = min([A_P_raster(:);AP_P_raster(:)]);
        maxD = max([A_P_raster(:);AP_P_raster(:)]);
        for cc = 1:numberOfCols
            axes(ff.h_axes(rr,cc));
            if cc == 1
            	plotDistRaster(ff.h_axes(rr,cc),A_P,[minD maxD],sprintf('%s - %.2f',num2str(ccs(cn)),A_P_SI));
            end
            if cc == 2
                plotDistRaster(ff.h_axes(rr,cc),AP_P,[minD maxD],sprintf('%s - %.2f',num2str(ccs(cn)),AP_P_SI));
            end
%             if cc == 3
%                 plotDistRaster(ff.h_axes(rr,cc),P,[minD maxD],sprintf('%s - %.2f',num2str(ccs(cn)),P_SI));
%             end
%             if cc == 4
%                 Ar = getTimeRaster_simple(b,caSig,tsp,ap_a_t_onsets,ap_a_t_offsets,0);
%                 imagesc(Ar.sigRaster);colorbar;
%             end
        end
    end
    showCells(1000,ei,ccs(cng));
    display('Press key');
    gg = keyboardInput(gg,[1 numberOfGroups],[1 6],'');
    if gg < 0
        break;
    end
end
