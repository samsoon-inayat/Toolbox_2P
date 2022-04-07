function ff = plot_dist_rasters(R,ccsi,ff)
cellList = ccsi;
ddtt = 1;
A = R;
for cc = 1:length(ccsi)
    cn = cellList(cc);
    thisRaster = A.sp_rasters1(:,:,cn);
    minRasters(cc) = min(thisRaster(:));
    maxRasters(cc) = max(thisRaster(:));
end
minmin = min(minRasters(:));
maxmax = max(maxRasters(:));
for rr = 1
    for cc = 1:length(ccsi)
        cn = cellList(cc);
        axes(ff.h_axes(rr,cc));
%         thisRaster = A.rasters(:,:,cn);
        thisRaster = A.sp_rasters1(:,:,cn);
        mSig = nanmean(thisRaster);
%         ft = fittype(A.formula);
        xs = A.xs;
%         coeff = A.coeff(:,cn);
        fitplot = gauss_fit(1:length(xs),A.gauss_fit_on_mean.coefficients_Rs_mean(cn,1:3),A.gauss_fit_on_mean.gauss1Formula);
        imagesc(thisRaster,[min(thisRaster(:)) max(thisRaster(:))]);hold on;
        if rr ==1 
            plot(size(thisRaster,1)*fitplot/max(fitplot),'linewidth',0.2,'color','m');
        end
        box off;
        plot(size(thisRaster,1)*normalizeSignal(nanmean(thisRaster)),'b','linewidth',0.2);
        if rr == 1
%         text(size(thisRaster,2)+size(thisRaster,2)/20,1,sprintf('zMI = %.2f, Rs = %.2f',A.info_metrics.ShannonMI_Zsh(cn),...
%             A.gauss_fit_on_mean.coefficients_Rs_mean(cn,4)),'FontSize',5,'color','k','rotation',90);
        else
%        text(size(thisRaster,2)+size(thisRaster,2)/15,2,sprintf('Max FR %d - zMI = %.2f',round(max(thisRaster(:))),A.SI(cn)),'FontSize',5,'color','k','rotation',90);
        end
        cols = size(thisRaster,2);
        colsHalf = ceil(cols/2);
        if rr == 1
%             text(1,size(thisRaster,1)+1.5,sprintf('Cell %d, (MFR = %.0f)',cn,round(max(thisRaster(:)))),'FontSize',5,'color','k');
            text(1,size(thisRaster,1)+1.5,sprintf('%d A.U.',round(max(thisRaster(:)))),'FontSize',5,'color','k');
        end
        text(colsHalf,size(thisRaster,1)+1.5,sprintf('zMI = %.2f ',R.info_metrics.ShannonMI_Zsh(cn)),'FontSize',5,'color','k');
        set(gca,'FontSize',6,'FontWeight','Normal','linewidth',0.5,'Ydir','normal','TickDir','out');
        colormap jet;
        if rr == 1
%             cols = size(thisRaster,2);
%             xticks = [1:floor(cols/2):size(thisRaster,2)];
%             set(gca,'XTick',xticks,'XTickLabels',{'0','75','150'});
            hx = xlabel('Distance (cm)');
        end
        
        
        xticks = [1 colsHalf cols];
        xs = R.xs;
        xs = round(xs,1); 
        xticks = unique(xticks);
        if ~strcmp(R.marker_name,'airD')
            for xli = 1:length(xticks)
                if xli == 1 || xli == length(xticks)
                    xticklabels{xli} = sprintf('%.0f',xs(xticks(xli)));
                else
                    xticklabels{xli} = sprintf('%.0f',xs(xticks(xli)));
                end
            end
        else
            xticklabels = {'0','75','150'};
        end
        set(gca,'XTick',xticks,'XTickLabels',xticklabels);
        
        if ~isempty(strfind(R.context_info,'airID'))
                for bb = 1:length(R.lastBin)
                    xvalbin = R.lastBin(bb)-1;
                    plot([xvalbin xvalbin]+0.15,[bb-0.75 bb+0.75],'c','linewidth',0.25);
                end
            end
        
        if cc >1
            set(gca,'YTick',[]);
        end
        if cc == 1
%             set(gca,'YTick',[1 5 10]);
%             h = ylabel('Trial #');
%             changePosition(h,[3 0 0]);
%             text
        end
        
        if cc == length(ff.h_axes) && rr == 1
            hca = gca;
%             hc = putColorBar(hca,[-0.05 0 -0.05 0],{'0','Max FR(AU)'},6,'northoutside',[0.15 0.32 0.05 0.32]);
        end
        format_axes(gca);
    end
end
