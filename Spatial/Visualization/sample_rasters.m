function ff = sample_rasters(R,ccs,ff)

if iscell(R)
    for cns = 1:length(ccs)
        tff = ff;
        for rr = 1:length(R)
            tff.h_axes = ff.h_axes(cns,rr);
            tff.axesPos = ff.axesPos(cns,rr);
            sample_rasters(R{rr},ccs(cns),tff);
        end
    end
    return;
end

n = 0;
if isfield(R,'fromFrames')
    ff = plot_time_rasters(R,ccs,ff);
else
    ff = plot_dist_rasters(R,ccs,ff);
end
cm = colormap(gray);
cm = flipud(cm(1:size(cm,1),:));
colormap(cm);



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
%         fitplot = gauss_fit(1:length(xs),R.gauss_fit_on_mean.coefficients_Rs_mean(cn,1:3),R.gauss_fit_on_mean.gauss1Formula);
        if strfind(R.marker_name,'motion')
            imagesc(thisRaster,[min(thisRaster(:)) max(thisRaster(:))]);hold on;
        else
            imagesc(thisRaster,[min(thisRaster(:)) max(thisRaster(:))]);hold on;
        end
%         if rr ==1 
%             plot(size(thisRaster,1)*fitplot/max(fitplot),'linewidth',1,'color','m');
%         end
        if isfield(R,'resp')
            if isfield(R.resp,'cis')
                cis = R.resp.cis;
                plot([cis(2,1) cis(2,1)]+1,[0 size(thisRaster,1)+1],'linewidth',0.1,'color','m');
                if isempty(strfind(R.marker_name,'motion'))
                    plot([cis(1,3) cis(1,3)]+1,[0 size(thisRaster,1)+1],'linewidth',0.1,'color','c');
                end
            end
        end
        box off;
        try
        text(size(thisRaster,2)+size(thisRaster,2)/15,2,sprintf('zMI = %.2f',R.info_metrics.ShannonMI_Zsh(cn)),'FontSize',3.5,'color','k','rotation',90);
        catch
        end
%         text(size(thisRaster,2)+size(thisRaster,2)/15,2,sprintf('Max FR = %d',round(max(thisRaster(:)))),'FontSize',5,'color','k','rotation',90);
        text(1,size(thisRaster,1)+1.5,sprintf('Cell %d, (MFR = %d)',cn,round(max(thisRaster(:)))),'FontSize',5,'color','k');
        set(gca,'FontSize',6,'linewidth',0.25,'Ydir','normal','TickDir','out');
        colormap jet;
            xticks = [1:10:size(thisRaster,2)];
            xs = R.xs;
            if isfield(R,'resp')
                if isfield(R.resp,'cis')
                    xs = xs - xs(cis(1,2));
                    xticks = [cis(1,:) size(thisRaster,2)];
                else
                    xticks = [1:round(size(thisRaster,2)/3):size(thisRaster,2)];
                end
            else
                xticks = [1:round(size(thisRaster,2)/3):size(thisRaster,2)];
            end
%             xticks
            if ~ismember(size(thisRaster,2),xticks)
                xticks = [xticks size(thisRaster,2)];
            end
            xs = round(xs,1); 
            xticks = unique(xticks);
            set(gca,'XTick',xticks,'XTickLabels',xs(xticks));
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
            ntrials = size(thisRaster,1);
            set(gca,'YTick',[1  ntrials]);
            h = ylabel('Trials');
        end
        
        if cc == length(ff.h_axes)
            hca = gca;
            ff.hc = putColorBar(hca,[-0.07 0 -0.05 0],{'0','Max FR (AU)'},6,'northoutside',[0.15 0.3 0.05 0.3]);
        end
        cols = size(thisRaster,2);
        colsHalf = ceil(cols/2);
        ts = xs;
%         set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    end
end


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
            plot(size(thisRaster,1)*fitplot/max(fitplot),'linewidth',0.25,'color','m');
        end
        box off;
        if rr == 1
        text(size(thisRaster,2)+size(thisRaster,2)/20,1,sprintf('zMI = %.2f, Rs = %.2f',A.info_metrics.ShannonMI_Zsh(cn),...
            A.gauss_fit_on_mean.coefficients_Rs_mean(cn,4)),'FontSize',3.5,'color','k','rotation',90);
        else
       text(size(thisRaster,2)+size(thisRaster,2)/15,2,sprintf('Max FR %d - zMI = %.2f',round(max(thisRaster(:))),A.SI(cn)),'FontSize',4,'color','k','rotation',90);
        end
        if rr == 1
            text(1,size(thisRaster,1)+1.5,sprintf('Cell %d, (MFR = %.0f)',cn,round(max(thisRaster(:)))),'FontSize',5,'color','k');
        end
        set(gca,'FontSize',6,'FontWeight','Normal','linewidth',0.5,'Ydir','normal','TickDir','out');
        colormap jet;
        if rr == 1
            cols = size(thisRaster,2);
            xticks = [1:floor(cols/2):size(thisRaster,2)];
            set(gca,'XTick',xticks,'XTickLabels',A.xs(xticks));
            hx = xlabel('Position (cm)');
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
        
        if cc == length(ff.h_axes) && rr == 1
            hca = gca;
%             hc = putColorBar(hca,[0 0 -0.05 0],{'0','Max FR'},6,'northoutside',[0.15 0.22 0.05 0.22]);
        end
    end
end
