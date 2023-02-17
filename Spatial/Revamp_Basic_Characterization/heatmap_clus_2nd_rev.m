    %% for heat maps exc and inh only air
%     figdim = 3;
%     hf = get_figure(5,[9 2 figdim figdim]);
            
    all_cells_list = all_exc_inh;
    all_cells_list = all_gV;
    sh = 0;
    good_FR = circshift(all_cells_list,sh,2);
    txl = circshift(event_type,sh,2);
    seq_nums = circshift(1:size(all_cells_list,2),sh,2);
  
    [OIo,mOIo,semOIo,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(good_FR,0.5,0.05);
    mOI = mCI; semOI = semCI;


    mUni = nanmean(uni,3); mUni1 = tril(mUni,-1) + tril(mUni,-1)'; mUni2 = triu(mUni,1) + triu(mUni,1)'; mmUni = min(mUni(:)); MmUni = max(mUni(:));
    semUni = nanstd(uni,[],3)/sqrt(5); semUni1 = tril(semUni,-1) + tril(semUni,-1)'; semUni2 = triu(semUni,1) + triu(semUni,1)'; msemUni = min(semUni(:)); MsemUni = max(semUni(:));
%     figure(1000);clf;subplot 131;imagesc(mUni,[mmUni MmUni]);set(gca,'YDir','normal');subplot 132;imagesc(mUni1,[mmUni MmUni]);set(gca,'YDir','normal');subplot 133;imagesc(mUni2,[mmUni MmUni]);set(gca,'YDir','normal');
    
    
   ff = makeFigureRowsCols(107,[5 3 6.9 2.5],'RowsCols',[1 2],'spaceRowsCols',[0.2 0.2],'rightUpShifts',[0.12 0.19],'widthHeightAdjustment',[-10 -300]);
    set(gcf,'color','w');     ylims = [0 1];
    stp = 0.57; widths = 0.5*([2 2 2 2 0.4 0.4]+1.5); gap = 0.99*1.75; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
    gap1 = 0.32; widths1 = 1*0.75; height1 = widths1 - 0.05;
    
    mats = {mOIo,mUni}; semmats = {semOIo,semUni};
    titles = {'Jaccard Distance','Complementation'};
    for ii = 1:length(mats)
        mOI = mats{ii}; semOI = semmats{ii};
        sz = size(mOI,1);        oM = ones(size(mOI));
        mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN; semOI(mask1 == 1) = NaN;
        maxI = max([mOI(:);semOI(:)]);            minI = min([mOI(:);semOI(:)]);

        mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 

        imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
        imAlpha(mask1 == 1) = 0;
        axes(ff.h_axes(1,ii));
        im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
        set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
        format_axes(gca);
        set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
        set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
        for rr = 1:size(mCI,1)
            for cc = 1:size(mCI,1)
                if rr == cc
                    text(rr-0.25,cc,sprintf('%.0f',mCI(rr,cc)),'FontSize',5);
                end
            end
        end
    %     plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 150.5],'w','linewidth',0.1); 
    %     plot([0 150.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
        set(gca,'Ydir','normal');ytickangle(15);        box on
        format_axes(gca);
        [hc,hca] = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.07 0.05 0.07 0.05]);
        colormap jet
        htit = title('Mean'); changePosition(htit,[0 0.07 0]); set(htit,'FontSize',6);
        ht = set_axes_top_text_no_line(gcf,gca,titles{ii},[0 -0.01 0 0]);set(ht,'FontSize',6,'FontWeight','Bold')

%         save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean_%d.pdf',ntrials,sh),600);
        hai = ff.h_axes(1,ii);
        pos = get(hai,'Position'); units = get(hai,'Units');
        ha = axes('Position',pos,'Visible','on','Units',units); poshca =  get(hca,'Position');
        set(ha,'Position',[hai.Position(1)+hai.Position(3)+gap1 (hai.Position(2)+hai.Position(4)-height1) widths1 height1]);
%         axes(ff.h_axes(1,jj(ii)+ii+1));
        im1 = imagesc(semOI,[minI,maxI]);    im1.AlphaData = imAlpha;
        set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
        format_axes(gca);
        set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
        set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',[],'yticklabels',[],'Ydir','reverse'); xtickangle(75);
        changePosition(gca,[0.0 0 -0.05 0]);
        set(gca,'Ydir','normal');
        htit = title('SEM'); changePosition(htit,[0 -0.2 0]); set(htit,'FontSize',6);
        box on
        format_axes(gca);
        colormap jet
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('OI_Map_%d_%d.pdf',ntrials,sh),600);
    %% distribution of JD
    all_bvs = []; all_jdvals = [];
    for an = 1:5
        jdvals = squareform(OIo{an},'tovector');
        all_minjdv(an) = min(jdvals);
        all_jdvals(an,:) = jdvals;
        bins = 0.2:0.05:1;
        [bvs,xs] = hist(jdvals,bins); bvs = bvs/sum(bvs);
        all_bvs(an,:) = bvs;
    end
    mbvs = mean(all_bvs); sembvs = std(all_bvs)/sqrt(5);
    ff = makeFigureWindow__one_axes_only(1000,[7 5 2.5 1.5],[0.15 0.2 0.75 0.75]);
    plot(xs,mbvs);
    shadedErrorBar(xs,mbvs,sembvs);
    ylabel('Fraction');xlabel('Jaccard Distance');
    format_axes(gca); set(gca,'xlim',[min(bins) max(bins)]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('JD.pdf'),600);
   %% for clustering
   ff = makeFigureRowsCols(107,[5 3 6.9 2.5],'RowsCols',[2 3],'spaceRowsCols',[0.3 -0.02],'rightUpShifts',[1 0.18],'widthHeightAdjustment',[-10 -300]);
    set(gcf,'color','w');     ylims = [0 1];
    stp = 0.25; widths = [2 2 2 0.4 0.4 0.4]+0.75; gap = 0.18; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
    stp = 0.35; widths = 0.5*([2 2 2 0.4 0.4 0.4]+1.95); gap = 0.2; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
    %
    mCI = 1-mOIo;
    hms = {mCI,mUni1,mUni2};
    ahc_col_th = 0.7;
    
    hms1 = {100-mCI,mUni1,mUni2};
    for hi = 1:length(hms1)
        mOI1 = hms1{hi}; mOI1(mask1==1) = 0; Di = squareform(mOI1,'tovector');%pdist(mOI1,@naneucdist);
        tree = linkage(Di,'average'); [c,d] = cophenet(tree,Di); r = corr(Di',d','type','spearman');    leafOrder = optimalleaforder(tree,Di);
        hf = figure(100000000); %     leafOrder1 = leafOrder([1:3 10:12 4:9]);
    %     leafOrder1 = circshift(leafOrder,3);
        figure(hf);clf
        [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);ylims = ylim; xlims = xlim;    set(gca,'xtick',[],'ytick',[]);
        set(gcf,'units','inches'); set(gcf,'Position',[5 2 0.9 0.5]); close(hf);
        axes(ff.h_axes(1,hi));
        [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);
        set(H,'linewidth',0.5); set(gca,'xticklabels',txl(leafOrder));xtickangle(45); format_axes(gca);
        if hi == 1
            hx = ylabel({'Euc. Dist.'});changePosition(hx,[0 0 0]);
        end
    %     xlim([xlims(1)+0.5 xlims(2)-0.5]);
%         changePosition(gca,[0.0 0.0 0.07 0.05]); text(0.5,ylims(2)+0,sprintf('CC = %.2f',c),'FontSize',6);
        hcct = set_axes_top_text_no_line(ff.hf,gca,sprintf('CC = %.2f',c),[0 -0.07 0 0]); set(hcct,'FontWeight','Normal')
        LOs{1,hi} = leafOrder;
    end
    
    for hi = 1:length(hms)
        mOI1 = hms{hi}; mOI1(mask1==1) = NaN; Di = pdist(mOI1,@naneucdist);
        tree = linkage(Di,'average'); [c,d] = cophenet(tree,Di); r = corr(Di',d','type','spearman');    leafOrder = optimalleaforder(tree,Di);
        hf = figure(100000000); %     leafOrder1 = leafOrder([1:3 10:12 4:9]);
    %     leafOrder1 = circshift(leafOrder,3);
        figure(hf);clf
        [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);ylims = ylim; xlims = xlim;    set(gca,'xtick',[],'ytick',[]);
        set(gcf,'units','inches'); set(gcf,'Position',[5 2 0.9 0.5]); close(hf);
        axes(ff.h_axes(2,hi));
        [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);
        set(H,'linewidth',0.5); set(gca,'xticklabels',txl(leafOrder));xtickangle(45); format_axes(gca);
        if hi == 1
            hx = ylabel({'Euc. Dist.'});changePosition(hx,[0 0 0]);
        end
    %     xlim([xlims(1)+0.5 xlims(2)-0.5]);
%         changePosition(gca,[0.0 0.0 0.07 0.05]); text(0.5,ylims(2)+0,sprintf('CC = %.2f',c),'FontSize',6);
        hcct = set_axes_top_text_no_line(ff.hf,gca,sprintf('CC = %.2f',c),[0 0.01 0 0]); set(hcct,'FontWeight','Normal')
    %     set_axes_top_text(ff.hf,ff.h_axes(1),sprintf('Cophenetic Correlation = %.2f',c));
        LOs{2,hi} = leafOrder;
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('OI_Map_cluster_%d.pdf',sh),600);

    
    %% for heat maps all_gV
%     figdim = 3;
%     hf = get_figure(5,[9 2 figdim figdim]);
            
    all_cells_list = all_gV;
%     all_cells_list = all_gV_A;
    sh = 0;
    good_FR = circshift(all_cells_list,sh,2);
    txl = circshift(event_type,sh,2);
    seq_nums = circshift(1:size(all_cells_list,2),sh,2);
  
    [OIo,mOIo,semOIo,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(good_FR,0.5,0.05);
    mOI = mCI; semOI = semCI;


    mUni = nanmean(uni,3); mUni1 = tril(mUni,-1) + tril(mUni,-1)'; mUni2 = triu(mUni,1) + triu(mUni,1)'; mmUni = min(mUni(:)); MmUni = max(mUni(:));
    semUni = nanstd(uni,[],3)/sqrt(5); semUni1 = tril(semUni,-1) + tril(semUni,-1)'; semUni2 = triu(semUni,1) + triu(semUni,1)'; msemUni = min(semUni(:)); MsemUni = max(semUni(:));
%     figure(1000);clf;subplot 131;imagesc(mUni,[mmUni MmUni]);set(gca,'YDir','normal');subplot 132;imagesc(mUni1,[mmUni MmUni]);set(gca,'YDir','normal');subplot 133;imagesc(mUni2,[mmUni MmUni]);set(gca,'YDir','normal');
    
    
   ff = makeFigureRowsCols(107,[5 3 6.9 2.5],'RowsCols',[1 2],'spaceRowsCols',[0.2 0.2],'rightUpShifts',[0.12 0.19],'widthHeightAdjustment',[-10 -300]);
    set(gcf,'color','w');     ylims = [0 1];
    stp = 0.57; widths = 0.5*([2 2 2 2 0.4 0.4]+1.5); gap = 0.99*1.75; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
    gap1 = 0.32; widths1 = 1*0.75; height1 = widths1 - 0.05;
    
    mats = {mCI,mUni}; semmats = {semCI,semUni};
    titles = {'Conjunction','Complementation'};
    for ii = 1:length(mats)
        mOI = mats{ii}; semOI = semmats{ii};
        sz = size(mOI,1);        oM = ones(size(mOI));
        mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN; semOI(mask1 == 1) = NaN;
        maxI = max([mOI(:);semOI(:)]);            minI = min([mOI(:);semOI(:)]);

        mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 

        imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
        imAlpha(mask1 == 1) = 0;
        axes(ff.h_axes(1,ii));
        im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
        set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
        format_axes(gca);
        set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
        set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
        for rr = 1:size(mCI,1)
            for cc = 1:size(mCI,1)
                if rr == cc
                    text(rr-0.25,cc,sprintf('%.0f',mCI(rr,cc)),'FontSize',5);
                end
            end
        end
    %     plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 150.5],'w','linewidth',0.1); 
    %     plot([0 150.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
        set(gca,'Ydir','normal');ytickangle(15);        box on
        format_axes(gca);
        [hc,hca] = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.07 0.05 0.07 0.05]);
        colormap jet
        htit = title('Mean'); changePosition(htit,[0 0.07 0]); set(htit,'FontSize',6);
        ht = set_axes_top_text_no_line(gcf,gca,titles{ii},[0 -0.01 0 0]);set(ht,'FontSize',6,'FontWeight','Bold')

%         save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean_%d.pdf',ntrials,sh),600);
        hai = ff.h_axes(1,ii);
        pos = get(hai,'Position'); units = get(hai,'Units');
        ha = axes('Position',pos,'Visible','on','Units',units); poshca =  get(hca,'Position');
        set(ha,'Position',[hai.Position(1)+hai.Position(3)+gap1 (hai.Position(2)+hai.Position(4)-height1) widths1 height1]);
%         axes(ff.h_axes(1,jj(ii)+ii+1));
        im1 = imagesc(semOI,[minI,maxI]);    im1.AlphaData = imAlpha;
        set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
        format_axes(gca);
        set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
        set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',[],'yticklabels',[],'Ydir','reverse'); xtickangle(75);
        changePosition(gca,[0.0 0 -0.05 0]);
        set(gca,'Ydir','normal');
        htit = title('SEM'); changePosition(htit,[0 -0.2 0]); set(htit,'FontSize',6);
        box on
        format_axes(gca);
        colormap jet
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('OI_Map_%d_%d.pdf',ntrials,sh),600);
   %% for clustering
   ff = makeFigureRowsCols(107,[5 3 6.9 2.5],'RowsCols',[2 3],'spaceRowsCols',[0.3 -0.02],'rightUpShifts',[1 0.18],'widthHeightAdjustment',[-10 -300]);
    set(gcf,'color','w');     ylims = [0 1];
    stp = 0.25; widths = [2 2 2 0.4 0.4 0.4]+0.75; gap = 0.18; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
    stp = 0.35; widths = 0.5*([2 2 2 0.4 0.4 0.4]+1.95); gap = 0.2; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
    %
    hms = {mCI,mUni1,mUni2};
    ahc_col_th = 0.7;
    
    hms1 = {100-mCI,mUni1,mUni2};
    for hi = 1:length(hms1)
        mOI1 = hms1{hi}; mOI1(mask1==1) = 0; Di = squareform(mOI1,'tovector');%pdist(mOI1,@naneucdist);
        tree = linkage(mOI1,'average'); [c,d] = cophenet(tree,Di); r = corr(Di',d','type','spearman');    leafOrder = optimalleaforder(tree,Di);
        hf = figure(100000000); %     leafOrder1 = leafOrder([1:3 10:12 4:9]);
    %     leafOrder1 = circshift(leafOrder,3);
        figure(hf);clf
        [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);ylims = ylim; xlims = xlim;    set(gca,'xtick',[],'ytick',[]);
        set(gcf,'units','inches'); set(gcf,'Position',[5 2 0.9 0.5]); close(hf);
        axes(ff.h_axes(1,hi));
        [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);
        set(H,'linewidth',0.5); set(gca,'xticklabels',txl(leafOrder));xtickangle(45); format_axes(gca);
        if hi == 1
            hx = ylabel({'Euc. Dist.'});changePosition(hx,[0 0 0]);
        end
    %     xlim([xlims(1)+0.5 xlims(2)-0.5]);
%         changePosition(gca,[0.0 0.0 0.07 0.05]); text(0.5,ylims(2)+0,sprintf('CC = %.2f',c),'FontSize',6);
        hcct = set_axes_top_text_no_line(ff.hf,gca,sprintf('CC = %.2f',c),[0 -0.07 0 0]); set(hcct,'FontWeight','Normal')
    end
    
    for hi = 1:length(hms)
        mOI1 = hms{hi}; mOI1(mask1==1) = NaN; Di = pdist(mOI1,@naneucdist);
        tree = linkage(Di,'average'); [c,d] = cophenet(tree,Di); r = corr(Di',d','type','spearman');    leafOrder = optimalleaforder(tree,Di);
        hf = figure(100000000); %     leafOrder1 = leafOrder([1:3 10:12 4:9]);
    %     leafOrder1 = circshift(leafOrder,3);
        figure(hf);clf
        [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);ylims = ylim; xlims = xlim;    set(gca,'xtick',[],'ytick',[]);
        set(gcf,'units','inches'); set(gcf,'Position',[5 2 0.9 0.5]); close(hf);
        axes(ff.h_axes(2,hi));
        [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);
        set(H,'linewidth',0.5); set(gca,'xticklabels',txl(leafOrder));xtickangle(45); format_axes(gca);
        if hi == 1
            hx = ylabel({'Euc. Dist.'});changePosition(hx,[0 0 0]);
        end
    %     xlim([xlims(1)+0.5 xlims(2)-0.5]);
%         changePosition(gca,[0.0 0.0 0.07 0.05]); text(0.5,ylims(2)+0,sprintf('CC = %.2f',c),'FontSize',6);
        hcct = set_axes_top_text_no_line(ff.hf,gca,sprintf('CC = %.2f',c),[0 0.01 0 0]); set(hcct,'FontWeight','Normal')
    %     set_axes_top_text(ff.hf,ff.h_axes(1),sprintf('Cophenetic Correlation = %.2f',c));
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('OI_Map_cluster_%d.pdf',sh),600);
   
%% for clustering
   ff = makeFigureRowsCols(107,[5 3 6.9 2.5],'RowsCols',[2 3],'spaceRowsCols',[0.3 -0.02],'rightUpShifts',[1 0.18],'widthHeightAdjustment',[-10 -300]);
    set(gcf,'color','w');     ylims = [0 1];
    stp = 0.25; widths = [2 2 2 0.4 0.4 0.4]+0.75; gap = 0.18; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
    stp = 0.35; widths = 0.5*([2 2 2 0.4 0.4 0.4]+1.95); gap = 0.2; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
    %
    hms = {mOIo,mUni1,mUni2};
    ahc_col_th = 0.7;
    
    hms1 = {100-mCI,mUni1,mUni2};
    hms1 = {1-mOIo,mUni1,mUni2};
    mOIo1 = mOIo;mOIo1(mask1==1) = NaN;
    for hi = 1:length(hms1)
        mOI1 = hms1{hi}; mOI1(mask1==1) = 0; Di = squareform(mOI1,'tovector');%pdist(mOI1,@naneucdist);
        tree = linkage(Di,'average'); [c,d] = cophenet(tree,Di); r = corr(Di',d','type','spearman');    leafOrder = optimalleaforder(tree,Di);
        hf = figure(100000000); %     leafOrder1 = leafOrder([1:3 10:12 4:9]);
    %     leafOrder1 = circshift(leafOrder,3);
        figure(hf);clf
        [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);ylims = ylim; xlims = xlim;    set(gca,'xtick',[],'ytick',[]);
        set(gcf,'units','inches'); set(gcf,'Position',[5 2 0.9 0.5]); close(hf);
        axes(ff.h_axes(1,hi));
        [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);
        set(H,'linewidth',0.5); set(gca,'xticklabels',txl(leafOrder));xtickangle(45); format_axes(gca);
        if hi == 1
            hx = ylabel({'Euc. Dist.'});changePosition(hx,[0 0 0]);
        end
    %     xlim([xlims(1)+0.5 xlims(2)-0.5]);
%         changePosition(gca,[0.0 0.0 0.07 0.05]); text(0.5,ylims(2)+0,sprintf('CC = %.2f',c),'FontSize',6);
        hcct = set_axes_top_text_no_line(ff.hf,gca,sprintf('CC = %.2f',c),[0 -0.07 0 0]); set(hcct,'FontWeight','Normal')
        LOs{1,hi} = leafOrder;
    end
    
    for hi = 1:length(hms)
        mOI1 = hms{hi}; mOI1(mask1==1) = NaN; Di = pdist(mOI1,@naneucdist);
        tree = linkage(Di,'average'); [c,d] = cophenet(tree,Di); r = corr(Di',d','type','spearman');    leafOrder = optimalleaforder(tree,Di);
        hf = figure(100000000); %     leafOrder1 = leafOrder([1:3 10:12 4:9]);
    %     leafOrder1 = circshift(leafOrder,3);
        figure(hf);clf
        [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);ylims = ylim; xlims = xlim;    set(gca,'xtick',[],'ytick',[]);
        set(gcf,'units','inches'); set(gcf,'Position',[5 2 0.9 0.5]); close(hf);
        axes(ff.h_axes(2,hi));
        [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);
        set(H,'linewidth',0.5); set(gca,'xticklabels',txl(leafOrder));xtickangle(45); format_axes(gca);
        if hi == 1
            hx = ylabel({'Euc. Dist.'});changePosition(hx,[0 0 0]);
        end
    %     xlim([xlims(1)+0.5 xlims(2)-0.5]);
%         changePosition(gca,[0.0 0.0 0.07 0.05]); text(0.5,ylims(2)+0,sprintf('CC = %.2f',c),'FontSize',6);
        hcct = set_axes_top_text_no_line(ff.hf,gca,sprintf('CC = %.2f',c),[0 0.01 0 0]); set(hcct,'FontWeight','Normal')
    %     set_axes_top_text(ff.hf,ff.h_axes(1),sprintf('Cophenetic Correlation = %.2f',c));
        LOs{2,hi} = leafOrder;
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('OI_Map_cluster_%d.pdf',sh),600);

    
     %% dissimilarity matrics for heat maps exc and inh only air
%     figdim = 3;
%     hf = get_figure(5,[9 2 figdim figdim]);
%     event_type = {'AOn-Exc','AOn-Inh','AOff-Exc','AOff-Inh','Arb-Exc','Arb-Inh','AOn-Exc','AOn-Inh','AOff-Exc','AOff-Inh','Arb-Exc','Arb-Inh'};
    all_cells_list = all_exc_inh;
%     all_cells_list = all_gV;
    sh = 0;
    good_FR = circshift(all_cells_list,sh,2);
    txl = circshift(event_type,sh,2);
    seq_nums = circshift(1:size(all_cells_list,2),sh,2);
  
    [OIo,mOIo,semOIo,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(good_FR,0.5,0.05);
    mOI = mCI; semOI = semCI;


    mUni = nanmean(uni,3); mUni1 = tril(mUni,-1) + tril(mUni,-1)'; mUni2 = triu(mUni,1) + triu(mUni,1)'; mmUni = min(mUni(:)); MmUni = max(mUni(:));
    semUni = nanstd(uni,[],3)/sqrt(5); semUni1 = tril(semUni,-1) + tril(semUni,-1)'; semUni2 = triu(semUni,1) + triu(semUni,1)'; msemUni = min(semUni(:)); MsemUni = max(semUni(:));
%     figure(1000);clf;subplot 131;imagesc(mUni,[mmUni MmUni]);set(gca,'YDir','normal');subplot 132;imagesc(mUni1,[mmUni MmUni]);set(gca,'YDir','normal');subplot 133;imagesc(mUni2,[mmUni MmUni]);set(gca,'YDir','normal');
    
    
    ff = makeFigureRowsCols(107,[5 3 6.9 2.5],'RowsCols',[1 3],'spaceRowsCols',[0.2 0.2],'rightUpShifts',[0.12 0.19],'widthHeightAdjustment',[-10 -300]);
    set(gcf,'color','w');     ylims = [0 1];
    stp = 0.56; widths = 0.5*([2 2 2 2 0.4 0.4]+1.5); gap = 0.24*1.75; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
    gap1 = 0.32; widths1 = 1*0.75; height1 = widths1 - 0.05;
    
    mats = {mCI,mUni1,mUni2}; semmats = {semCI,semUni1,semUni2};
    titles = {'Conjuction','Complementation1','Complementation2'};
    for ii = 1:length(mats)
        mOIt = mats{ii}; semOI = semmats{ii};
        mOI = pdist(mOIt); mOI = squareform(mOI);
        sz = size(mOI,1);        oM = ones(size(mOI));
        mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN; semOI(mask1 == 1) = NaN;
        maxI = max([mOI(:);semOI(:)]);            minI = min([mOI(:);semOI(:)]);

        mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 

        imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
        imAlpha(mask1 == 1) = 0;
        axes(ff.h_axes(1,ii));
        im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
        set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
        format_axes(gca);
        set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
        if ii == 1
            set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
        else
            set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',[],'Ydir','reverse'); xtickangle(45);
        end
        for rr = 1:size(mCI,1)
            for cc = 1:size(mCI,1)
                if rr == cc
                    text(rr-0.25,cc,sprintf('%.0f',mCI(rr,cc)),'FontSize',5);
                end
            end
        end
    %     plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 150.5],'w','linewidth',0.1); 
    %     plot([0 150.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
        set(gca,'Ydir','normal');ytickangle(15);        box on
        format_axes(gca);
        [hc,hca] = putColorBar(gca,[0.04 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.07 0.05 0.07 0.05]);
        colormap jet
        htit = title(titles{ii}); changePosition(htit,[0 0.07 0]); set(htit,'FontSize',6);
        if ii == 1
            ht = set_axes_top_text_no_line(gcf,gca,'Dissimilarity Matrices',[-0.07 0.05 0 -0.07]);set(ht,'FontSize',7,'FontWeight','Bold','Color',[0.25 0 0.25])
        end

% %         save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean_%d.pdf',ntrials,sh),600);
%         hai = ff.h_axes(1,ii);
%         pos = get(hai,'Position'); units = get(hai,'Units');
%         ha = axes('Position',pos,'Visible','on','Units',units); poshca =  get(hca,'Position');
%         set(ha,'Position',[hai.Position(1)+hai.Position(3)+gap1 (hai.Position(2)+hai.Position(4)-height1) widths1 height1]);
% %         axes(ff.h_axes(1,jj(ii)+ii+1));
%         im1 = imagesc(semOI,[minI,maxI]);    im1.AlphaData = imAlpha;
%         set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
%         format_axes(gca);
%         set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
%         set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',[],'yticklabels',[],'Ydir','reverse'); xtickangle(25);
%         changePosition(gca,[0.0 0 -0.05 0]);
%         set(gca,'Ydir','normal');
%         htit = title('SEM'); set(htit,'FontSize',5); changePosition(htit,[0 -0.2 0]); 
%         box on
%         format_axes(gca);
        colormap jet
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('OI_Map_%d_%d.pdf',ntrials,sh),600);

    
    %% light Figure 4 dissimilarity matrics for heat maps exc and inh only air
%     figdim = 3;
%     hf = get_figure(5,[9 2 figdim figdim]);
%     event_type = {'AOn-Exc','AOn-Inh','AOff-Exc','AOff-Inh','Arb-Exc','Arb-Inh','AOn-Exc','AOn-Inh','AOff-Exc','AOff-Inh','Arb-Exc','Arb-Inh'};
%     all_cells_list = all_exc_inh;
    all_cells_list = all_gV;
    sh = 0;
    good_FR = circshift(all_cells_list,sh,2);
    txl = circshift(event_type,sh,2);
    seq_nums = circshift(1:size(all_cells_list,2),sh,2);
  
    [OIo,mOIo,semOIo,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(good_FR,0.5,0.05);
    mOI = mCI; semOI = semCI;


    mUni = nanmean(uni,3); mUni1 = tril(mUni,-1) + tril(mUni,-1)'; mUni2 = triu(mUni,1) + triu(mUni,1)'; mmUni = min(mUni(:)); MmUni = max(mUni(:));
    semUni = nanstd(uni,[],3)/sqrt(5); semUni1 = tril(semUni,-1) + tril(semUni,-1)'; semUni2 = triu(semUni,1) + triu(semUni,1)'; msemUni = min(semUni(:)); MsemUni = max(semUni(:));
%     figure(1000);clf;subplot 131;imagesc(mUni,[mmUni MmUni]);set(gca,'YDir','normal');subplot 132;imagesc(mUni1,[mmUni MmUni]);set(gca,'YDir','normal');subplot 133;imagesc(mUni2,[mmUni MmUni]);set(gca,'YDir','normal');
    
    
    ff = makeFigureRowsCols(107,[5 3 6.9 1.5],'RowsCols',[1 3],'spaceRowsCols',[0.2 0.2],'rightUpShifts',[0.12 0.19],'widthHeightAdjustment',[-10 -400]);
    set(gcf,'color','w');     ylims = [0 1];
    stp = 0.34; widths = 0.45*([2 2 2 2 0.4 0.4]); gap = 0.99*1.25; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
    gap1 = 0.2; widths1 = 1*0.45; height1 = widths1 - 0.05;
    
    mats = {mCI,mUni,mOIo}; semmats = {semCI,semUni,semOIo};
    titles = {'Conjuction','Complementation','Jaccard Distance'};
    for ii = 1:length(mats)
        mOI = mats{ii}; semOI = semmats{ii};
%         if ii == 4
%             mOI = pdist(mOI); mOI = squareform(mOI);
%         end
        sz = size(mOI,1);        oM = ones(size(mOI));
        mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN; semOI(mask1 == 1) = NaN;
        maxI = max([mOI(:);semOI(:)]);            minI = min([mOI(:);semOI(:)]);

        mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 

        imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
        imAlpha(mask1 == 1) = 0;
        axes(ff.h_axes(1,ii));
        im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
        set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
        format_axes(gca);
        set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
        set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
        for rr = 1:size(mCI,1)
            for cc = 1:size(mCI,1)
                if rr == cc
                    if ii < 3
                    text(rr-0.25,cc,sprintf('%.0f',mCI(rr,cc)),'FontSize',5);
                    end
                end
            end
        end
    %     plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 150.5],'w','linewidth',0.1); 
    %     plot([0 150.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
        set(gca,'Ydir','normal');ytickangle(15);        box on
        format_axes(gca);
        [hc,hca] = putColorBar(gca,[0.07 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.07 0.05 0.07 0.08]);
        colormap jet
        htit = title('Mean'); changePosition(htit,[0 0.07 0]); set(htit,'FontSize',5,'FontWeight','Normal');
%         if ii == 1
            ht = set_axes_top_text_no_line(gcf,gca,titles{ii},[-0.01 0.07 0.05 0]);set(ht,'FontSize',7,'FontWeight','Bold','Color',[0.25 0 0.25])
%         end

%         save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean_%d.pdf',ntrials,sh),600);
        hai = ff.h_axes(1,ii);
        pos = get(hai,'Position'); units = get(hai,'Units');
        ha = axes('Position',pos,'Visible','on','Units',units); poshca =  get(hca,'Position');
        set(ha,'Position',[hai.Position(1)+hai.Position(3)+gap1 (hai.Position(2)+hai.Position(4)-height1) widths1 height1]);
%         axes(ff.h_axes(1,jj(ii)+ii+1));
        im1 = imagesc(semOI,[minI,maxI]);    im1.AlphaData = imAlpha;
        set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
        format_axes(gca);
        set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
        set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',[],'yticklabels',[],'Ydir','reverse'); xtickangle(25);
        changePosition(gca,[0.0 0 -0.05 0]);
        set(gca,'Ydir','normal');
        htit = title('SEM'); set(htit,'FontSize',5,'FontWeight','Normal'); changePosition(htit,[0 -0.2 0]); 
        box on
%         format_axes(gca);
        colormap jet
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('OI_Map_%d_%d.pdf',ntrials,sh),600);
  
    %%
    for ii = 1:5
        jd(ii,:) = squareform(OIo{ii},'tovector');
    end
    
    avar = jd;
    [within,dvn,xlabels,withinD] = make_within_table({'pop'},[3]); withinD3 = withinD;
    dataT = make_between_table({avar},dvn);
    ra = RMA(dataT,within,{0.05,{'lsd','hsd','bonferroni'}});
    ra.ranova
    print_for_manuscript(ra)
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'pop','bonferroni'},[1 1 1]);
    xdata = make_xdata([size(good_FR,2)],[1]);
%     h(h==1) = 0;
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
    set(gcf,'color','w');    set(gcf,'Position',[10 3 1.6 1.25]);
    MY = 1.1; ysp = 0.09; mY = 0; % responsive cells
    stp = 0.43; widths = [0.4 0.4 0.4 0.4 0.4 0.4]+0.1; gap = 0.09;
    adjust_axes(ff,[mY MY],stp,widths,gap,{'Jaccard Distance'});
    % s = generate_shades(length(bins)-1);
    tcolors = [mData.colors(1:2);mData.dcolors(8)];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'B-NB','B-M','NB-M'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
%     set_axes_top_text(ff.hf,ff.h_axes(1),'Complementary',[-0.04 0 0 0]);
    
    save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs_3_factors.pdf',ntrials),600);