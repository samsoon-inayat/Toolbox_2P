function ff = sample_rasters(R,ccs,ff)
n = 0;
if isfield(R,'fromFrames')
    ff = plot_time_rasters(R,ccs,ff);
else
    ff = plot_dist_rasters(R,ccs,ff);
end


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
        imagesc(thisRaster,[min(thisRaster(:)) max(thisRaster(:))]);hold on;
%         if rr ==1 
%             plot(size(thisRaster,1)*fitplot/max(fitplot),'linewidth',1,'color','m');
%         end
        if isfield(R.resp,'cis')
            cis = R.resp.cis;
            plot([cis(2,1) cis(2,1)]+1,[0 size(thisRaster,1)+1],'linewidth',0.5,'color','r');
            plot([cis(1,3) cis(1,3)]+1,[0 size(thisRaster,1)+1],'linewidth',0.5,'color','m');
        end
        box off;
%         text(size(thisRaster,2)+size(thisRaster,2)/15,2,sprintf('Max FR %d - zMI = %.2f',round(max(thisRaster(:))),R.info_metrics.ShannonMI_Zsh(cn)),'FontSize',4,'color','k','rotation',90);
        text(size(thisRaster,2)+size(thisRaster,2)/15,2,sprintf('Max FR = %d',round(max(thisRaster(:)))),'FontSize',5,'color','k','rotation',90);
        text(size(thisRaster,1)/2,size(thisRaster,1)+1.25,sprintf('Cell %d',cn),'FontSize',6,'color','k');
        set(gca,'FontSize',6,'FontWeight','Bold','linewidth',0.5,'Ydir','normal','TickDir','out');
        colormap jet;
            xticks = [1:10:size(thisRaster,2)];
            xs = R.xs;
            if isfield(R.resp,'cis')
                xs = xs - xs(cis(1,2));
                xticks = [cis(1,:) size(thisRaster,2)];
            end
            xs = round(xs);
            set(gca,'XTick',xticks,'XTickLabels',xs(xticks));
            if strfind(R.context_info,'air') & strfind(R.context_info,'D')
                for bb = 1:length(R.lastBin)
                    xvalbin = R.lastBin(bb)-1;
                    plot([xvalbin xvalbin]+0.15,[bb-0.75 bb+0.75],'r','linewidth',1);
                end
            end
                hx = xlabel('Time (sec)');

        if cc >1
            set(gca,'YTick',[]);
        end
        if cc == 1
            set(gca,'YTick',[1 5 10]);
            h = ylabel('Trials');
        end
        
        if cc == length(ff.h_axes)
            hca = gca;
            hc = putColorBar(hca,[-0.07 0 -0.05 0],{'0','Max FR (AU)'},6,'northoutside',[0.15 0.3 0.05 0.3]);
        end
        cols = size(thisRaster,2);
        colsHalf = ceil(cols/2);
        ts = xs;
%         set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    end
end
colormap parula;

function ff = plot_dist_rasters(R,ccsi,ff)
cellList = ccsi;
ddtt = 1;
A = R;
for cc = 1:4
    cn = cellList(cc);
    thisRaster = A.sp_rasters1(:,:,cn);
    minRasters(cc) = min(thisRaster(:));
    maxRasters(cc) = max(thisRaster(:));
end
minmin = min(minRasters(:));
maxmax = max(maxRasters(:));
for rr = 1
    for cc = 1:4
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
            plot(size(thisRaster,1)*fitplot/max(fitplot),'linewidth',0.5,'color','m');
        end
        box off;
        if rr == 1
        text(size(thisRaster,2)+size(thisRaster,2)/20,1,sprintf('Max FR %d - zMI = %.2f - Rs = %.2f',round(max(thisRaster(:))),A.info_metrics.ShannonMI_Zsh(cn),...
            A.gauss_fit_on_mean.coefficients_Rs_mean(cn,4)),'FontSize',3,'color','k','rotation',90);
        else
       text(size(thisRaster,2)+size(thisRaster,2)/15,2,sprintf('Max FR %d - zMI = %.2f',round(max(thisRaster(:))),A.SI(cn)),'FontSize',4,'color','k','rotation',90);
        end
        if rr == 1
            text(size(thisRaster,1)/2,size(thisRaster,1)+1.25,sprintf('Cell %d',cn),'FontSize',6,'color','k');
        end
        set(gca,'FontSize',6,'FontWeight','Bold','linewidth',0.75,'Ydir','normal','TickDir','out');
        colormap jet;
        if rr == 1
            cols = size(thisRaster,2);
            xticks = [1:floor(cols/2):size(thisRaster,2)];
            set(gca,'XTick',xticks,'XTickLabels',A.xs(xticks));
%                 if cc == 2
                hx = xlabel('Position (cm)');
%                     changePosition(hx,[-50 0 0]);                    
%                 end
        else
            xticks = [1:10:size(thisRaster,2)];
            xs = 0:0.2:100;
            set(gca,'XTick',xticks,'XTickLabels',xs(xticks));
            for bb = 1:length(A.lastBin)
                xvalbin = A.lastBin(bb)-1;
                plot([xvalbin xvalbin]+0.15,[bb-0.75 bb+0.75],'r','linewidth',1.5);
            end
            
                hx = xlabel('Time (sec)');
%                     changePosition(hx,[-20 0 0]);                    
        end

        if cc >1
            set(gca,'YTick',[]);
        end
        if cc == 1
            set(gca,'YTick',[1 5 10]);
            h = ylabel('Trials');
%             changePosition(h,[3 0 0]);
%             text
        end
        
        if cc == 4 && rr == 1
            hca = gca;
%             hc = putColorBar(hca,[0 0 -0.05 0],{'0','Max FR'},6,'northoutside',[0.15 0.22 0.05 0.22]);
        end
    end
end
colormap parula;
