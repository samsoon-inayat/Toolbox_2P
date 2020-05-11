function plotRastersMulti(Rs,ccs,spca,pmean,keyFlag)
% ccs = 1:size(R1.rasters,3);
if ~exist('keyFlag','var')
    keyFlag = 1;
end
if length(ccs) < 4
    numberOfRows = length(ccs);
else
    numberOfRows = 4;
end
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

ff = makeFigureRowsCols(1066,[NaN 0.1 7 5],'RowsCols',[numberOfRows length(Rs)],...
    'spaceRowsCols',[0.1 0.05],'rightUpShifts',[0.05 0.07],'widthHeightAdjustment',...
    [-70 -115]);

gg = 1;
while 1
    for rr = 1:numberOfRows
        for cc = 1:numberOfCols
            cni = groupIndices(rr,cc,gg);
            if isnan(cni)
                continue;
            end
            cn = ccs(cni);
    %         cng(rr) = cn;
            if isnan(cn)
                continue;
            end
            for rii = 1:length(Rs)
                if spca
                    raster{rii} = Rs{rii}.caRasters(:,:,cn);
                else
                    raster{rii} = Rs{rii}.rasters(:,:,cn);
                end
                if isfield(Rs{rii},'dist')
                    mins(rii) = min(raster{rii}(:));
                    maxs(rii) = max(raster{rii}(:));
                else
                    minsT(rii) = min(raster{rii}(:));
                    maxsT(rii) = max(raster{rii}(:));
                end
            end
            if exist('mins','var')
                minD = min(mins);
                maxD = mean(maxs);
            end
            if exist('minsT','var')
                minDT = min(minsT);
                maxDT = mean(maxsT);
            end
            for rii = 1:length(Rs)
                if isfield(Rs{rii},'dist')
                    plotDistRaster(ff.h_axes(rr,rii),Rs{rii},cn,spca,[minD maxD],sprintf('%s - %.2f - %.2f - %d - %d',num2str(ccs(cni)),Rs{rii}.SI(cn),Rs{rii}.rs(cn),...
                        round(Rs{rii}.centers(cn)),round(Rs{rii}.pws(cn))),pmean);
                else
                    plotTimeRaster(ff.h_axes(rr,rii),Rs{rii},cn,spca,[minDT maxDT],sprintf('%s - %.2f',num2str(ccs(cni)),Rs{rii}.SI(cn)),pmean);
                end
%                 if rr == 1
%                     if isfield(Rs{rii},'dist')
%                         text(Rs{rii}.dist(5),length(Rs{rii}.onsets)+5,Rs{rii}.name,'FontSize',12,'FontWeight','Bold');
%                     else
%                         text(Rs{rii}.cells(1).times(1),length(Rs{rii}.onsets)+5,Rs{rii}.name,'FontSize',12,'FontWeight','Bold');
%                     end
%                 end
                if rr == numberOfRows
                    if isfield(Rs{rii},'dist')
                        xlabel('Distance (cm)');
                    else
                        xlabel('Time (secs)');
                    end
                end

            end
            
%             plotDistRaster(ff.h_axes(rr,2),R2,cn,0,[minD maxD],sprintf('%s - %.2f - %.2f - %.2f',num2str(ccs(cni)),R2.SI(cn),R2.pcw(cni),R2.pcc(cni)),pmean);
%             plotDistRaster(ff.h_axes(rr,3),R3,cn,0,[minD maxD],sprintf('%s - %.2f - %.2f - %.2f',num2str(ccs(cni)),R3.SI(cn),R3.pcw(cni),R3.pcc(cni)),pmean);
%             plotDistRaster(ff.h_axes(rr,4),R4,cn,0,[minD maxD],sprintf('%s - %.2f - %.2f - %.2f',num2str(ccs(cni)),R4.SI(cn),R4.pcw(cni),R4.pcc(cni)),pmean);
        end
    end
%     showCells(1000,ei,ccs(cng));
%     save2pdf('temp.pdf',gcf,600);
    if keyFlag
        display('Press key');
        gg = keyboardInput(gg,[1 numberOfGroups],[1 6],'');
        if gg < 0
            break;
        end
    else
        break;
    end
end
