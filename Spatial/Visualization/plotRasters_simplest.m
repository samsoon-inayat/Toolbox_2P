function plotRasters_simplest(R,ccs)
rasters = R.sp_rasters1;
if ~exist('ccs','var')
    ccs = find(R.resp.vals);
end

% if length(ccs) == 1
%     figure(1005);
%     clf;
%     thisRaster = rasters(:,:,1);
%     mSig = nanmean(thisRaster);
%     imagesc(thisRaster,[0 0.5*max(thisRaster(:))]);colorbar;hold on;
%     plot(size(thisRaster,1)*mSig/max(mSig),'linewidth',0.75,'color','m');
% %     title(sprintf('%d',cn));
%     box off;
%     set(gca,'FontSize',10,'FontWeight','Normal','linewidth',0.75,'Ydir','normal','TickDir','out');
%     return;
% end

numberOfRows = 4;
numberOfCols = 4;
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

ff = makeFigureRowsCols(1005,[NaN 0.1 7 5],'RowsCols',[numberOfRows numberOfCols],...
    'spaceRowsCols',[0.1 0.05],'rightUpShifts',[0.05 0.07],'widthHeightAdjustment',...
    [-70 -115]);

gg = 1;
while 1
    for rr = 1:numberOfRows
        for cc = 1:numberOfCols
            cni = groupIndices(rr,cc,gg);
            if isnan(cni)
                cla(ff.h_axes(rr,cc));
                set(ff.h_axes(rr,cc),'visible','off');
                continue;
            end
            set(ff.h_axes(rr,cc),'visible','on');
            cn = ccs(cni);
            thisRaster = rasters(:,:,cn);
            mSig = nanmean(thisRaster);
            axes(ff.h_axes(rr,cc));
            mtr = 1;%mode(thisRaster(thisRaster(:)>0));
%             try
            f_TR = (fft2(thisRaster));
            af_TR = abs(f_TR);
            mf_TR = angle(f_TR);
            haa = findHaFD(af_TR,1:10);
            ham = findHaFD(mf_TR,1:10);
            htR = R.fractal_dim.HaFD(cn);
            [~,CRR] = findPopulationVectorPlot(thisRaster',1:size(thisRaster,2));
            haCRR = findHaFD(CRR,1:size(CRR,1));
            try
                imagesc(CRR);
%                 imagesc(thisRaster,[0 mtr*max(thisRaster(:))]);
            catch
                continue;
            end
            colorbar;hold on;
            plot(size(thisRaster,1)*mSig/max(mSig),'linewidth',0.5,'color','w');
            title(sprintf('%d - %.3f - %.3f - %.3f - %.3f',cn,htR,haa,ham,haCRR));
%             title(sprintf('%d-SI(%.2f)-Center(%.1f)-MaxFR(%.1f)',cn,A.SI(cn),A.centers(cn),max(thisRaster(:))));
            box off;
            set(gca,'FontSize',10,'FontWeight','Normal','linewidth',0.75,'Ydir','normal','TickDir','out');
            if isfield(R.resp,'cis')
                cis = R.resp.cis;
                plot([cis(2,1) cis(2,1)],[0 11],'linewidth',1.5,'color','r');
                plot([cis(1,3) cis(1,3)],[0 11],'linewidth',1.5,'color','m');
            end
%             catch
%             end
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
