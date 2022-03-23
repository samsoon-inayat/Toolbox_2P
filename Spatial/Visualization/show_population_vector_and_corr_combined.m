function ff = show_population_vector_and_corr(mData,ff,Rs,allP,allC,minmaxCorr,maxBin,cbf)

if ~exist('cbf','var')
    cbf = 1;
end

cp = 1;

if ~isempty(allP)
    FS = mData.axes_font_size;
    for sii = 1:size(allP,2)
        P = allP{sii};
        R = Rs{sii};

        %% first row
        axes(ff.h_axes(1,sii));
        if cp
            changePosition(gca,[0 0.05 -0.091 -0.1]);
        end
        if ~isempty(maxBin)
            P = P(:,1:maxBin);
        end
        imagesc(P);
        box off;
        if isfield(R.resp,'cis')
                hold on;
                cis = R.resp.cis;
                if isempty(strfind(R.marker_name,'airIT'))
                plot([cis(2,1) cis(2,1)]+1,[0 size(P,1)+1],'linewidth',0.25,'color','m');
                end
                if strcmp(R.marker_name,'air55T') 
                    plot([cis(1,3) cis(1,3)]+1,[0 size(P,1)+1],'linewidth',0.25,'color','c');
                end
        end
        if sii == 1
            h = ylabel('Cell No.'); %   changePosition(h,[0 0 0]);
        end
    %     text(3,size(P,1)+round(size(P,1)/7),sprintf('%s',raster_labels{sii}),'FontSize',FS,'FontWeight','Normal');
        set(gca,'Ydir','Normal','linewidth',0.5,'FontSize',FS,'FontWeight','Bold','YTick',[1]);
        set(gca,'YTick',[1 size(P,1)]);
        set(gca,'XTick',[]);
%         text(0,size(P,1)+5,{'Pop. Activity'},'FontSize',5);
%         ht = title('Pop. Activity');
%         set(ht,'FontSize',5,'FontWeight','Normal');

        %% 2nd row
        axes(ff.h_axes(2,sii));
        dec = -0.09;
        if cp
            changePosition(gca,[0.0 0.05 dec dec]);
        end
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
        set(gca,'Ydir','Normal','linewidth',0.5,'FontSize',FS,'FontWeight','Bold','TickDir','out');
        if sii == 1 || (sii == 4 && strcmp(R.marker_name,'airIT'))
           if  isfield(R,'fromFrames')
                h = ylabel('Time (sec)');    changePosition(h,[0 0 0]);
           else
               h = ylabel('Distance (cm)');    changePosition(h,[0 0 0]);
            end
        end
        cols = size(P,2);
        colsHalf = round(cols/2);
%         ts = round(R.xs(1:cols));
        set(gca,'XTick',[]);
        if isfield(R,'fromFrames')
            h = xlabel('Time (sec)');%    changePosition(h,[0 0 0]);
        else
%             h = xlabel('Distance (cm)');%    changePosition(h,[0 0 0]);
        end
        cols = size(P,2);
        colsHalf = ceil(cols/2);
        ts = floor(R.xs+R.bin_width/2);
        if isfield(R.resp,'cis')
            xs = R.xs;
            xs = xs - xs(cis(1,2));
            if ~strfind(R.marker_name,'motion')
                xticks = [cis(1,:) size(corrPlot,2)];
            else
                xticks = [cis(1,1:2) size(corrPlot,2)];
            end
            xs = round(xs);
            set(gca,'XTick',xticks,'XTickLabels',xs(xticks));
            set(gca,'YTick',xticks,'YTickLabel',xs(xticks));
            hold on;
            cis = R.resp.cis;
%             plot([cis(2,1) cis(2,1)]+1,[0 size(corrPlot,1)+1],'linewidth',0.5,'color','r');
%             plot([cis(1,3) cis(1,3)]+1,[0 size(corrPlot,1)+1],'linewidth',0.5,'color','m');
%             plot([0 size(corrPlot,2)+1],[cis(2,1) cis(2,1)]+1,'linewidth',0.5,'color','r');
%             plot([0 size(corrPlot,2)+1],[cis(1,3) cis(1,3)]+1,'linewidth',0.5,'color','m');
        else
            if max(ts) < 50
                tlabels = {'0','7.5','15'};
            else
                if strcmp(R.marker_name,'airD')
                    tlabels = {'0','75','150'};
                else
                    tlabels = [];
                end
            end
            set(gca,'XTick',[1 colsHalf cols],'XTickLabel',tlabels);
            set(gca,'YTick',[1 colsHalf cols],'YTickLabel',tlabels);
        end
%         ht = title('Pop. Correlation');
%         set(ht,'FontSize',5,'FontWeight','Normal');
    end

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
    if cbf
    hc = putColorBar(ff.h_axes(1,size(allP,2)),[0.0 0.03 0 -0.05],[0 maxs(1)],6,'eastoutside',[0.07 0.07 0.1 0.1]);
    hc = putColorBar(ff.h_axes(2,size(allP,2)),[0.0 0.03 0 -0.05],minmaxCorr,6,'eastoutside',[0.07 0.07 0.1 0.1]);
    % hc = putColorBar(ff.h_axes(3,4),[0.0 0.03 0 -0.05],[mIa maxs(3)],6,'eastoutside',[0.07 0.07 0.1 0.1]);
    end
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
        set(gca,'Ydir','Normal','linewidth',0.5,'FontSize',FS,'FontWeight','Bold','TickDir','out');
        if sii == 1 || (sii == 4 && strcmp(R.marker_name,'airIT'))
           if  isfield(R,'fromFrames')
            h = ylabel('Time (sec)');    changePosition(h,[0 0 0]);
           else
               h = ylabel('Distance (cm)');    changePosition(h,[0 0 0]);
            end
        end
        cols = size(C,2);
        colsHalf = round(cols/2);
        ts = floor(R.xs(1:cols)+R.bin_width/2);
        if isfield(R,'fromFrames')
            h = xlabel('Time (sec)');%    changePosition(h,[0 0 0]);
        else
%             h = xlabel('Distance (cm)');%    changePosition(h,[0 0 0]);
        end
        cols = size(C,2);
        colsHalf = ceil(cols/2);
        ts = floor(R.xs(1:cols)+R.bin_width/2);
        if isfield(R.resp,'cis')
            cis = R.resp.cis;
            xs = R.xs;
            xs = xs - xs(cis(1,2));
            if ~strfind(R.marker_name,'motion')
                xticks = [cis(1,:) size(corrPlot,2)];
            else
                xticks = [cis(1,1:2) size(corrPlot,2)];
            end
            xs = round(xs);
            
            set(gca,'XTick',xticks,'XTickLabels',xs(xticks));
            set(gca,'YTick',xticks,'YTickLabel',xs(xticks));

        else
             if max(ts) < 50
                tlabels = {'0','7.5','15'};
            else
                tlabels = {'0','75','150'};
            end
            set(gca,'XTick',[1 colsHalf cols],'XTickLabel',tlabels);
            set(gca,'YTick',[1 colsHalf cols],'YTickLabel',tlabels);
        end
            set(gca,'Ydir','Normal','linewidth',0.5,'FontSize',FS,'FontWeight','Bold');
            if sii == 1
            ht = title('Avg. Pop. Correlation');
        set(ht,'FontSize',5,'FontWeight','Normal');
            end
    end

    mI = min(minC);
    maxs = [1 1 1];
    if isempty(minmaxCorr)
        minmaxCorr = [mI maxs(2)];
    end
    for ii = 1:size(allC,2)
        axes(ff.h_axes(1,ii)); caxis(minmaxCorr);

    %     axes(ff.h_axes(3,ii)); caxis([mIa maxs(3)]);
    end
    if cbf
    hc = putColorBar(ff.h_axes(1,size(allC,2)),[0.0 0.03 0 -0.05],minmaxCorr,6,'eastoutside',[0.07 0.07 0.1 0.1]);
    % hc = putColorBar(ff.h_axes(3,4),[0.0 0.03 0 -0.05],[mIa maxs(3)],6,'eastoutside',[0.07 0.07 0.1 0.1]);
    end
end


cm = colormap(gray);
cm = flipud(cm(1:size(cm,1),:));
colormap(cm);