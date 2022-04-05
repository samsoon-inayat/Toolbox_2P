function ff = plot_time_rasters(R,ccsi,ff)
cellList = ccsi;
for cc = 1:length(ff.h_axes)
    cn = cellList(cc);
    thisRaster = R.sp_rasters1(:,:,cn);
    minRasters(cc) = min(thisRaster(:));
    maxRasters(cc) = max(thisRaster(:));
end
minmin = min(minRasters(:));
maxmax = max(maxRasters(:));
for rr = 2
    for cc = 1:length(ff.h_axes)
        cn = cellList(cc);
        axes(ff.h_axes(1,cc));
        thisRaster = R.sp_rasters1(:,:,cn);
        mSig = nanmean(thisRaster);
%         ft = fittype(A.formula);
        xs = R.xs;
%         coeff = A.coeff(:,cn);
        fitplot = gauss_fit(1:length(xs),R.gauss_fit_on_mean.coefficients_Rs_mean(cn,1:3),R.gauss_fit_on_mean.gauss1Formula);
        if strfind(R.marker_name,'motion')
            imagesc(thisRaster,[min(thisRaster(:)) max(thisRaster(:))]);hold on;
        else
            imagesc(thisRaster,[min(thisRaster(:)) max(thisRaster(:))]);hold on;
        end
        plot(10*normalizeSignal(nanmean(thisRaster)),'b','linewidth',0.2);
%         if rr ==1 
            plot(size(thisRaster,1)*fitplot/max(fitplot),'linewidth',0.2,'color','m');
%         end
        if isfield(R,'resp')
            if isfield(R.resp,'cis')
                cis = R.resp.cis;
                if isempty(strfind(R.marker_name,'IT'))
                    plot([cis(2,1) cis(2,1)]+1,[0 size(thisRaster,1)+1],'linewidth',0.2,'color','m');
%                     if isempty(strfind(R.marker_name,'motion'))
%                         plot([cis(1,3) cis(1,3)]+1,[0 size(thisRaster,1)+1],'linewidth',0.1,'color','c');
%                     end
                end
            end
        end
        box off;
        try
%         text(size(thisRaster,2)+size(thisRaster,2)/15,2,sprintf('zMI = %.2f',R.info_metrics.ShannonMI_Zsh(cn)),'FontSize',5,'color','k','rotation',90);
        catch
        end
%         text(size(thisRaster,2)+size(thisRaster,2)/15,2,sprintf('Max FR = %d',round(max(thisRaster(:)))),'FontSize',5,'color','k','rotation',90);
%         text(1,size(thisRaster,1)+1.5,sprintf('Cell %d, (MFR = %d)',cn,round(max(thisRaster(:)))),'FontSize',5,'color','k');
        if max(thisRaster(:)) < 1
            text(1,size(thisRaster,1)+1.5,sprintf('%.1f A.U.',max(thisRaster(:))),'FontSize',5,'color','k');
        else
            text(1,size(thisRaster,1)+1.5,sprintf('%d A.U.',round(max(thisRaster(:)))),'FontSize',5,'color','k');
        end
        cols = size(thisRaster,2);
        colsHalf = ceil(cols/2);
        text(colsHalf,size(thisRaster,1)+1.5,sprintf('zMI = %.2f ',R.info_metrics.ShannonMI_Zsh(cn)),'FontSize',5,'color','k');
        set(gca,'FontSize',6,'linewidth',0.25,'Ydir','normal','TickDir','out');
        colormap jet;
        cols = size(thisRaster,2);
        colsHalf = ceil(cols/2);
        xticks = [1 colsHalf cols];
        xs = R.xs;
        xs = round(xs,1); 
        xticks = unique(xticks);
        if ~strcmp(R.marker_name,'airIT')
        for xli = 1:length(xticks)
            if xli == 1 || xli == length(xticks)
                xticklabels{xli} = sprintf('%.0f',xs(xticks(xli)));
            else
                xticklabels{xli} = sprintf('%.1f',xs(xticks(xli)));
            end
        end
        else
            xticklabels = {'0','7.5','15'};
        end
        set(gca,'XTick',xticks,'XTickLabels',xticklabels);
        if ~isempty(strfind(R.context_info,'airT'))
            for bb = 1:length(R.lastBin)
                xvalbin = R.lastBin(bb)-1;
                plot([xvalbin xvalbin]+0.15,[bb-0.75 bb+0.75],'c','linewidth',0.25);
            end
        end
        hx = xlabel('Time (sec)');

        if cc >1
            set(gca,'YTick',[]);
        end
        if cc == 1
%             h = ylabel('Trial #');
        end
        
        if cc == length(ff.h_axes)
            hca = gca;
%             ff.hc = putColorBar(hca,[-0.07 0 -0.05 0],{'0','Max FR (AU)'},6,'northoutside',[0.15 0.3 0.05 0.3]);
        end
        
        ts = xs;
        format_axes(gca);
%         set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    end
end
