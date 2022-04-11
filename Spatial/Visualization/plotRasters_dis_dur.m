function plotRasters_simplest(Rs,ccs)
R = Rs{1};
rasters = R.sp_rasters1;
if ~exist('ccs','var')
    ccs = find(R.resp.vals);% & R.resp.FR_based);
else
    if isempty(ccs)
        ccs = 1:length(R.resp.vals);
    end
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

ff = makeFigureRowsCols(1005,[0.5 2.5 19 7 ],'RowsCols',[numberOfRows 4],...
    'spaceRowsCols',[0.1 0.05],'rightUpShifts',[0.05 0.07],'widthHeightAdjustment',...
    [-70 -115]);


gg = 1;
while 1
    for rr = 1:numberOfRows
        for cc = 1:4
            cni = groupIndices(rr,1,gg);
            if isnan(cni)
                cla(ff.h_axes(rr,cc));
                set(ff.h_axes(rr,cc),'visible','off');
                continue;
            end
            set(ff.h_axes(rr,cc),'visible','on');
            cn = ccs(cni);
            R = Rs{cc};
            thisRaster = R.sp_rasters1(:,:,cn);
%             thisRaster = rasters(:,:,cn);
            mSig = nanmean(thisRaster);
            axes(ff.h_axes(rr,cc));
            mtr = 1;%mode(thisRaster(thisRaster(:)>0));
%             try
            f_TR = (fft2(thisRaster));
            af_TR = abs(f_TR);
            mf_TR = angle(f_TR);
            haa = findHaFD(af_TR,1:size(thisRaster,1));
            ham = findHaFD(mf_TR,1:size(thisRaster,1));
%             htR = R.fractal_dim.HaFD(cn);
%             [~,CRR] = findPopulationVectorPlot(thisRaster',1:size(thisRaster,2));
%             haCRR = findHaFD(CRR,1:size(CRR,1));
%             thisRaster = normalizeSignal(thisRaster,2);
%             thisRaster = thisRaster > 0;
            try
%                 imagesc(thisRaster);
%                 imagesc(CRR);
                imagesc(thisRaster,[0 max(thisRaster(:))]);
            catch
                imagesc(thisRaster);
%                 continue;
            end
            colormap jet;
            colorbar;hold on;
%             plot(size(thisRaster,1)*mSig/max(mSig),'linewidth',0.5,'color','w');
            try
%                 title(sprintf('%d - %.3f - %.3f',cn,R.info_metrics.ShannonMI_Zsh(cn),R.gauss_fit_on_mean.worked(cn)));
                cellrsp = R.resp.valsC(cn) & R.resp.FR_based710(cn);
                title(sprintf('%d (%.3f) %d',cn,R.info_metrics.ShannonMI_Zsh(cn),cellrsp));
%                 err
            catch
                title(sprintf('%d - %d',cn,R.resp.vals(cn)));
            end
%             title(sprintf('%d - %.3f - %.3f - %.3f - %.3f',cn,htR,haa,ham,haCRR));
%             title(sprintf('%d-SI(%.2f)-Center(%.1f)-MaxFR(%.1f)',cn,A.SI(cn),A.centers(cn),max(thisRaster(:))));
            box off;
            set(gca,'FontSize',10,'FontWeight','Normal','linewidth',0.75,'Ydir','normal','TickDir','out');
            if isfield(R,'resp')
            if isfield(R.resp,'cis')
                cis = R.resp.cis;
                plot([cis(2,1) cis(2,1)],[0 size(thisRaster,1)+1],'linewidth',1.5,'color','r');
                if ~strfind(R.marker_name,'motion')
                    plot([cis(1,3) cis(1,3)],[0 11],'linewidth',1.5,'color','g');
                end
            end
            end
            A = R;
            xs = 1:length(A.xs);
            if isfield(R,'context_info')
%             if (~isempty(strfind(R.context_info,'air')) | ~isempty(strfind(R.context_info,'belt'))) & ~isempty(strfind(R.context_info,'D'))
                fitplot = gauss_fit(xs,A.gauss_fit_on_mean.coefficients_Rs_mean(cn,1:3),A.gauss_fit_on_mean.gauss1Formula);
                plot(size(thisRaster,1)*fitplot/max(fitplot),'linewidth',0.5,'color','m');
                
%             end
            end
            ylims = ylim;
            lastBin = R.lastBin;
            if ~isempty(strfind(R.context_info,'airT')) || ~isempty(strfind(R.context_info,'airID'))
                for ii = 1:10
                    plot([lastBin(ii) lastBin(ii)],[ii-0.5 ii+0.5],'w');
                end
            else
                
            end
                
%             plot([49.5,49.5],ylims);
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
