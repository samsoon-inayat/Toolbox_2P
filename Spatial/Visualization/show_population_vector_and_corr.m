function ff = show_population_vector_and_corr(mData,ff,Rs,allP,allC,minmaxCorr,maxBin)

if ~isempty(allP)
    FS = mData.axes_font_size;
    for sii = 1:size(allP,2)
        P = allP{sii};
        R = Rs{sii};

        %% first row
        axes(ff.h_axes(1,sii));changePosition(gca,[0 0.05 -0.091 -0.1]);
        if ~isempty(maxBin)
            P = P(:,1:maxBin);
        end
        imagesc(P);
        box off;
        if isfield(R.resp,'cis')
                hold on;
                cis = R.resp.cis;
                plot([cis(2,1) cis(2,1)]+1,[0 size(P,1)+1],'linewidth',0.5,'color','r');
                plot([cis(1,3) cis(1,3)]+1,[0 size(P,1)+1],'linewidth',0.5,'color','m');
            end
        if sii == 1
            h = ylabel('Cell No.'); %   changePosition(h,[0 0 0]);
        end
    %     text(3,size(P,1)+round(size(P,1)/7),sprintf('%s',raster_labels{sii}),'FontSize',FS,'FontWeight','Normal');
        set(gca,'Ydir','Normal','linewidth',0.5,'FontSize',FS,'FontWeight','Bold','YTick',[1]);
        set(gca,'YTick',[1 size(P,1)]);
        set(gca,'XTick',[]);

        %% 2nd row
        axes(ff.h_axes(2,sii));
        dec = -0.09;
        changePosition(gca,[0.0 0.05 dec dec]);
        if ~isempty(maxBin)
            corrPlot = allC{sii}(1:maxBin,1:maxBin);
        else
            corrPlot = allC{sii};
        end
        imagesc(corrPlot,[-1 1]);
%         axis equal;
        minC(sii) = min(corrPlot(:));
        maxC(sii) = max(corrPlot(:));
        box off;
        set(gca,'Ydir','Normal','linewidth',0.5,'FontSize',FS,'FontWeight','Bold');
        if sii == 1
           if  isfield(R,'fromFrames')
            h = ylabel('Time (sec)');    changePosition(h,[0 0 0]);
           else
               h = ylabel('Position (cm)');    changePosition(h,[0 0 0]);
            end
        end
        cols = size(P,2);
        colsHalf = round(cols/2);
        ts = round(R.xs(1:cols));
        set(gca,'XTick',[]);
        if sii == 1
            set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1)-2 ts(colsHalf) ts(cols)+2]);
        else
            set(gca,'YTick',[]);
        end
        if isfield(R,'fromFrames')
            h = xlabel('Time (sec)');%    changePosition(h,[0 0 0]);
        else
            h = xlabel('Position (cm)');%    changePosition(h,[0 0 0]);
        end
        cols = size(P,2);
        colsHalf = ceil(cols/2);
        ts = floor(R.xs+R.bin_width/2);
        set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
        if sii == 1
            set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
        else
            set(gca,'YTick',[]);
    end
            set(gca,'Ydir','Normal','linewidth',0.5,'FontSize',FS,'FontWeight','Bold');
    end

    colormap parula
    mI = min(minC);
    maxs = [1 1 1];
    if isempty(minmaxCorr)
        minmaxCorr = [mI maxs(2)];
    end
    for ii = 1:size(allP,2)
        axes(ff.h_axes(1,ii)); caxis([0 maxs(1)]);
        axes(ff.h_axes(2,ii)); caxis(minmaxCorr);
    %     axes(ff.h_axes(3,ii)); caxis([mIa maxs(3)]);
    end

    hc = putColorBar(ff.h_axes(1,size(allP,2)),[0.0 0.03 0 -0.05],[0 maxs(1)],6,'eastoutside',[0.07 0.07 0.1 0.1]);
    hc = putColorBar(ff.h_axes(2,size(allP,2)),[0.0 0.03 0 -0.05],minmaxCorr,6,'eastoutside',[0.07 0.07 0.1 0.1]);
    % hc = putColorBar(ff.h_axes(3,4),[0.0 0.03 0 -0.05],[mIa maxs(3)],6,'eastoutside',[0.07 0.07 0.1 0.1]);
end

%%
if isempty(allP)
    FS = mData.axes_font_size;
    for sii = 1:size(allC,2)
        R = Rs{sii};
        if ~isempty(maxBin)
            corrPlot = allC{sii}(1:maxBin,1:maxBin);
        else
            corrPlot = allC{sii};
        end
        C = corrPlot;
        axes(ff.h_axes(1,sii));
        dec = -0.09;
        changePosition(gca,[0.0 0.05 dec dec]);
        imagesc(corrPlot,[-1 1]);
        minC(sii) = min(corrPlot(:));
        maxC(sii) = max(corrPlot(:));
        box off;
        set(gca,'Ydir','Normal','linewidth',0.5,'FontSize',FS,'FontWeight','Bold');
        if sii == 1
           if  isfield(R,'fromFrames')
            h = ylabel('Time (sec)');    changePosition(h,[0 0 0]);
           else
               h = ylabel('Position (cm)');    changePosition(h,[0 0 0]);
            end
        end
        cols = size(C,2);
        colsHalf = round(cols/2);
        ts = floor(R.xs(1:cols)+R.bin_width/2);
        set(gca,'XTick',[]);
        if sii == 1
            set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
        else
            set(gca,'YTick',[]);
        end
        if isfield(R,'fromFrames')
            h = xlabel('Time (sec)');%    changePosition(h,[0 0 0]);
        else
            h = xlabel('Position (cm)');%    changePosition(h,[0 0 0]);
        end
        cols = size(C,2);
        colsHalf = ceil(cols/2);
        ts = floor(R.xs(1:cols)+R.bin_width/2);
        set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
        if sii == 1
            set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
        else
            set(gca,'YTick',[]);
    end
            set(gca,'Ydir','Normal','linewidth',1,'FontSize',FS,'FontWeight','Bold');
    end

    colormap parula
    mI = min(minC);
    maxs = [1 1 1];
    if isempty(minmaxCorr)
        minmaxCorr = [mI maxs(2)];
    end
    for ii = 1:size(allC,2)
        axes(ff.h_axes(1,ii)); caxis(minmaxCorr);

    %     axes(ff.h_axes(3,ii)); caxis([mIa maxs(3)]);
    end
    hc = putColorBar(ff.h_axes(1,size(allC,2)),[0.0 0.03 0 -0.05],minmaxCorr,6,'eastoutside',[0.07 0.07 0.1 0.1]);
    % hc = putColorBar(ff.h_axes(3,4),[0.0 0.03 0 -0.05],[mIa maxs(3)],6,'eastoutside',[0.07 0.07 0.1 0.1]);
end


