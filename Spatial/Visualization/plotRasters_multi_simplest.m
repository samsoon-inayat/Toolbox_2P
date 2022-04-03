function plotRasters_multi(rasters,ccs,amdata)
if isempty(ccs)
    ccs = 1:size(rasters{1},3);
end

numberOfRows = length(rasters);
numberOfCols = 4;
graphsOnOneFigure = 1 * numberOfCols;
numberOfData = length(ccs);
numberOfGroups = ceil(numberOfData/graphsOnOneFigure);
numberOfFullGroups = floor(numberOfData/graphsOnOneFigure);
indices = NaN(1,(numberOfGroups*graphsOnOneFigure));
indices(1:numberOfData) = 1:numberOfData;
groupIndices = reshape(indices,1,numberOfCols,numberOfGroups);
% for gg = 1:numberOfGroups
%     groupIndices(:,:,gg) = groupIndices(:,:,gg)';
% end
% 
% ff = makeFigureRowsCols(1005,[NaN 0.1 7 5],'RowsCols',[numberOfRows numberOfCols],...
%     'spaceRowsCols',[0.1 0.05],'rightUpShifts',[0.05 0.07],'widthHeightAdjustment',...
%     [-70 -115]);

ff = makeFigureRowsCols(1005,[0.5 0.5 19 9 ],'RowsCols',[numberOfRows numberOfCols],...
    'spaceRowsCols',[0.1 0.05],'rightUpShifts',[0.05 0.07],'widthHeightAdjustment',...
    [-70 -115]);

gg = 1;
while 1
    for rr = 1:numberOfRows
        for cc = 1:numberOfCols
            cni = groupIndices(1,cc,gg);
            if isnan(cni)
                cla(ff.h_axes(rr,cc));
                set(ff.h_axes(rr,cc),'visible','off');
                continue;
            end
            set(ff.h_axes(rr,cc),'visible','on');
            cn = ccs(cni);
            Rs = rasters{rr};
            thisRaster = Rs.sp_rasters1(:,:,cn);
            if sum(thisRaster(:)) == 0
                continue;
            end
%             mdata = Rs.mdata;
            mSig = nanmean(thisRaster);
            axes(ff.h_axes(rr,cc));
            mtr = 1;%mode(thisRaster(thisRaster(:)>nanmean(thisRaster(:))));
            try
                imagesc(thisRaster,[0 mtr*max(thisRaster(:))]);colorbar;hold on;
            catch
                imagesc(thisRaster);colorbar;hold on;
            end
%             plot(size(thisRaster,1)*mSig/max(mSig),'linewidth',0.5,'color','w');
            title(sprintf('%d - %s',cn,''));
%             title(sprintf('%d-SI(%.2f)-Center(%.1f)-MaxFR(%.1f)',cn,A.SI(cn),A.centers(cn),max(thisRaster(:))));
            box off;
            set(gca,'FontSize',10,'FontWeight','Normal','linewidth',0.75,'Ydir','normal','TickDir','out');
%             xticks = [1:round(size(thisRaster,2)/3):size(thisRaster,2)];
%             set(gca,'XTick',xticks,'XTickLabels',mdata.xs(xticks));
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
