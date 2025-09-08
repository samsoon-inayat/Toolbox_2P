sel_inds = 1:18;
sel_inds = 19:36;
all_cellsnew_hmPC = all_cellsnew_hm(:,19:36);
all_cellsnew_hmPC = all_cellsnew_hm(:,1:18);
txl = txl_all(:,sel_inds);
txl = txl_all18;

    %%
    [OIo,mOIo,semOIo,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(cell_popsPC,0.5,0.05);
    mOI = mCI; semOI = semCI;
    % mCI = mOIo; semCI = semOIo;

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
        ht = set_axes_top_text_no_line(gcf,gca,titles{ii},[0 0.125 0 0]);set(ht,'FontSize',6,'FontWeight','Bold')

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
    save_pdf(ff.hf,mData.pdf_folder,sprintf('OI_Map.pdf'),600);
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
        % Di = pdist(mOI1,@naneucdist);
        tree = linkage(mOI1,'average'); [c,d] = cophenet(tree,Di); r = corr(Di',d','type','spearman');    leafOrder = optimalleaforder(tree,Di);
        hf = figure(100000000); %     leafOrder1 = leafOrder([1:3 10:12 4:9]);
    %     leafOrder1 = circshift(leafOrder,3);
        figure(hf);clf
        [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);ylims = ylim; xlims = xlim;    set(gca,'xtick',1:length(txl),'ytick',[]);
        set(gcf,'units','inches'); set(gcf,'Position',[5 2 0.9 0.5]); close(hf);
        axes(ff.h_axes(1,hi));
        [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);
        set(H,'linewidth',0.5); set(gca,'xtick',1:length(txl),'xticklabels',txl(leafOrder));xtickangle(45); format_axes(gca);
        if hi == 1
            hx = ylabel({'Euc. Dist.'});changePosition(hx,[0 0 0]);
        end
    %     xlim([xlims(1)+0.5 xlims(2)-0.5]);
%         changePosition(gca,[0.0 0.0 0.07 0.05]); text(0.5,ylims(2)+0,sprintf('CC = %.2f',c),'FontSize',6);
        hcct = set_axes_top_text_no_line(ff.hf,gca,sprintf('CC = %.2f',c),[0 0.05 0 0]); set(hcct,'FontWeight','Normal')
    end
    
    for hi = 1:length(hms)
        mOI1 = hms1{hi}; mOI1(mask1==1) = NaN; Di = pdist(mOI1,@naneucdist);
        tree = linkage(Di,'average'); [c,d] = cophenet(tree,Di); r = corr(Di',d','type','spearman');    leafOrder = optimalleaforder(tree,Di);
        hf = figure(100000000); %     leafOrder1 = leafOrder([1:3 10:12 4:9]);
    %     leafOrder1 = circshift(leafOrder,3);
        figure(hf);clf
        [H,T,TC] = dendrogram(tree,0,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);ylims = ylim; xlims = xlim;    set(gca,'xtick',1:length(txl),'ytick',[]);
        set(gcf,'units','inches'); set(gcf,'Position',[5 2 0.9 0.5]); close(hf);
        axes(ff.h_axes(2,hi));
        % figure(300);clf
        [H,T,TC] = dendrogram(tree,0,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);
        set(H,'linewidth',0.5); set(gca,'xtick',1:length(txl),'xticklabels',txl(leafOrder));xtickangle(45); format_axes(gca);
        if hi == 1
            hx = ylabel({'Euc. Dist.'});changePosition(hx,[0 0 0]);
        end
    %     xlim([xlims(1)+0.5 xlims(2)-0.5]);
%         changePosition(gca,[0.0 0.0 0.07 0.05]); text(0.5,ylims(2)+0,sprintf('CC = %.2f',c),'FontSize',6);
        hcct = set_axes_top_text_no_line(ff.hf,gca,sprintf('CC = %.2f',c),[0 0.05 0 0]); set(hcct,'FontWeight','Normal')
    %     set_axes_top_text(ff.hf,ff.h_axes(1),sprintf('Cophenetic Correlation = %.2f',c));
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);

    
    %% for heat maps all_gV
%     figdim = 3;
%     hf = get_figure(5,[9 2 figdim figdim]);
            
    all_cells_list = all_cellsnew;
%     all_cells_list = all_gV_A;
    sh = 0;
    good_FR = circshift(all_cells_list,sh,2);
    % txl = circshift(event_type,sh,2);
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
        hcct = set_axes_top_text_no_line(ff.hf,gca,sprintf('CC = %.2f',c),[0 0.07 0 0]); set(hcct,'FontWeight','Normal')
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
        hcct = set_axes_top_text_no_line(ff.hf,gca,sprintf('CC = %.2f',c),[0 0.07 0 0]); set(hcct,'FontWeight','Normal')
    %     set_axes_top_text(ff.hf,ff.h_axes(1),sprintf('Cophenetic Correlation = %.2f',c));
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('OI_Map_cluster_%d.pdf',sh),600);
   


    %%
    tx = string(txl_all(:)');
mi_cols = find(contains(tx,"-MI-"));
pc_cols = find(contains(tx,"-PC-"));

% Order within each lens (keeps your given order):
types = ["-T","-D","-S"];
configs = ["C3-A1","C3-A2","C4-A1","C4-A2","C5-A1","C5-A2"];

% 18 columns per lens in the order configs × types:
lens_order = [];
for c = configs
    for t = types
        lens_order(end+1) = find(contains(tx, c) & (contains(tx,"-MI-") | contains(tx,"-PC-")) & contains(tx,t), 1, 'first');
    end
end

% Split by lens using names (robust to order):
mi18  = find(contains(tx,"-MI-"));
pc18  = find(contains(tx,"-PC-"));
assert(numel(mi18)==18 && numel(pc18)==18,'Expect 18 MI and 18 PC columns');

%%
S_MI = pairwise_stats_from_X(all_cellsnew_hm1, mi18, 'Tail','right');
S_PC = pairwise_stats_from_X(all_cellsnew_hm1, pc18, 'Tail','right');
%%
% 1) Distance from mean Jaccard
J = S_MI.mean.J;                    % 18x18
D = 1 - J;                          % convert to distance
D(1:size(D,1)+1:end) = 0;           % zero diagonal
dvec = squareform(D);

% 2) Hierarchical clustering
Z   = linkage(dvec,'average');
ord = optimalleaforder(Z,dvec);

% 3) Plot dendrogram + heatmap
figure(1);
subplot(1,2,1);
dendrogram(Z,0,'Reorder',ord,'Orientation','top');
title('MI — clustering (1 - mean Jaccard)'); set(gca,'XTickLabel',[]);

subplot(1,2,2);
Jord = J(ord,ord);
imagesc(Jord,[0 0.3]); axis square; colormap(parula); colorbar
set(gca,'XTick',1:18,'YTick',1:18, ...
        'XTickLabel', tx(mi18(ord)), 'YTickLabel', tx(mi18(ord)));
xtickangle(45); title('MI — mean Jaccard heatmap');

%%
% Ensure these exist:
tx   = string(txl_all);
pc18 = find(contains(tx,"-PC-"));   % 18 PC columns

% --- Cluster PC on 1 - mean(Jaccard)
J_pc   = S_PC.mean.J;                        % 18x18
D_pc   = 1 - J_pc; 
D_pc(1:size(D_pc,1)+1:end) = 0;              % zero diagonal
dvec_pc= squareform(D_pc);
Z_pc   = linkage(dvec_pc,'average');
ord_pc = optimalleaforder(Z_pc,dvec_pc);

% --- FDR on meta-analysis p-values, get pooled OR
p_pc = S_PC.meta.p_two_sided;                % 18x18
OR_pc= S_PC.meta.OR;
q_pc = mafdr(p_pc(:),'BHFDR',true);
q_pc = reshape(q_pc,size(p_pc));

% --- Reorder everything by dendrogram order
JO   = J_pc(ord_pc,ord_pc);
qO   = q_pc(ord_pc,ord_pc);
ORO  = OR_pc(ord_pc,ord_pc);
labs = tx(pc18(ord_pc));

% --- Plot heatmap + significance overlays
figure(3);
subplot(1,2,1);
dendrogram(Z_pc,0,'Reorder',ord_pc,'Orientation','top');
title('PC — clustering (1 - mean Jaccard)'); set(gca,'XTickLabel',[]);

subplot(1,2,2);
imagesc(JO,[0 0.25]); axis square; colormap(parula); colorbar
set(gca,'XTick',1:18,'YTick',1:18,'XTickLabel',labs,'YTickLabel',labs);
xtickangle(45); title('PC — mean Jaccard with FDR-significant overlaps'); hold on
[ri,ci] = find(qO<=0.05 & ORO>1);  plot(ci,ri,'k.','MarkerSize',10);   % enriched
[rd,cd] = find(qO<=0.05 & ORO<1);  plot(cd,rd,'kx','MarkerSize',8,'LineWidth',1); % depleted
hold off

%%
% --- Choose p-values (two-sided meta) and effect size (pooled OR)
p_mi = S_MI.meta.p_two_sided;      % 18x18 matrix
OR_mi = S_MI.meta.OR;              % 18x18 pooled odds ratios

% --- FDR (Benjamini–Hochberg)
q_mi = mafdr(p_mi(:),'BHFDR',true);
q_mi = reshape(q_mi, size(p_mi));

% --- Reorder by the dendrogram order from Step 1
qO  = q_mi(ord, ord);
ORO = OR_mi(ord, ord);

% --- Plot (re-use Jord from Step 1 if you kept it)
figure(2);
imagesc(J(ord,ord), [0 0.25]); axis square; colormap(parula); colorbar
set(gca,'XTick',1:18,'YTick',1:18, ...
        'XTickLabel', tx(mi18(ord)), 'YTickLabel', tx(mi18(ord)));
xtickangle(45); title('MI — mean Jaccard with FDR-significant overlaps'); hold on

% Enriched (q<=0.05 & OR>1) -> dots; Depleted (q<=0.05 & OR<1) -> crosses
[ri,ci] = find(qO<=0.05 & ORO>1);
plot(ci,ri,'k.','MarkerSize',10);

[rd,cd] = find(qO<=0.05 & ORO<1);
plot(cd,rd,'kx','MarkerSize',8,'LineWidth',1);

hold off

%%
% --- Build aligned indices for MI and PC: configs × types (6 × 3 = 18)
tx = string(txl_all);
configs = ["C3-A1","C3-A2","C4-A1","C4-A2","C5-A1","C5-A2"];
types   = ["-T","-D","-S"];

mi18 = zeros(1,18); pc18 = zeros(1,18); k = 0;
for c = configs
  for t = types
    k = k+1;
    mi18(k) = find(contains(tx,c) & contains(tx,"-MI-") & contains(tx,t),1);
    pc18(k) = find(contains(tx,c) & contains(tx,"-PC-") & contains(tx,t),1);
  end
end
labs = tx(mi18);   % labels in this consistent order

% --- Per-animal Jaccard(APC vs MI) for each concept
A = size(all_cellsnew_hm1,1);      % #animals
jac = nan(A,18);

for a = 1:A
  for j = 1:18
    Aset = logical(all_cellsnew_hm1{a, pc18(j)});   % PC membership
    Bset = logical(all_cellsnew_hm1{a, mi18(j)});   % MI membership
    kAB  = sum(Aset & Bset);
    uAB  = sum(Aset | Bset);
    if uAB==0
        jac(a,j) = NaN;
    else
        jac(a,j) = kAB / uAB;
    end
  end
end

jac_mean = nanmean(jac,1);
jac_sem  = nanstd(jac,[],1) ./ sqrt(A);

% --- Plot mean±SEM per concept
figure(4); errorbar(1:18, jac_mean, jac_sem, 'o'); ylim([0 1]);
xlim([0.5 18.5]); grid on
set(gca,'XTick',1:18,'XTickLabel',labs); xtickangle(45);
ylabel('Jaccard (PC vs MI)'); title('Concordance per concept');


%
% --- Combine per-animal Fisher p (right-tailed) per concept
p_comb = nan(1,18);

for j = 1:18
  pvec = nan(A,1);
  for a = 1:A
    Aset = logical(all_cellsnew_hm1{a, pc18(j)});
    Bset = logical(all_cellsnew_hm1{a, mi18(j)});
    k  = sum(Aset & Bset);
    ao = sum(Aset & ~Bset);
    bo = sum(~Aset & Bset);
    nn = sum(~Aset & ~Bset);
    try
      [~, pvec(a)] = fishertest([k ao; bo nn], 'Tail','right');  % enrichment
    catch
      % fallback (rare): hypergeometric right tail
      N = numel(Aset); K = sum(Aset); M = sum(Bset);
      pvec(a) = 1 - hygecdf(k-1, N, K, M);
    end
  end
  ok = pvec>0 & pvec<1 & ~isnan(pvec);
  if any(ok)
    z = norminv(1 - pvec(ok));                     % right-tailed z
    Z = sum(z) / sqrt(sum(ok));                    % Stouffer
    p_comb(j) = 1 - normcdf(Z);                    % combined p (right)
  end
end

% --- FDR across 18 concepts and add stars
q = mafdr(p_comb,'BHFDR',true);
hold on
sig = find(q<=0.05);
plot(sig, jac_mean(sig) + 0.05, 'k*', 'MarkerSize', 8);   % mark significant
for s = sig
    text(s, min(jac_mean(s)+0.10, 0.98), sprintf('q=%.3f', q(s)), ...
         'HorizontalAlignment','center','FontSize',7);
end
hold off

%%
% --- Build aligned indices (same as Step 4). Skip if you already have mi18/pc18.
tx = string(txl_all);
configs = ["C3-A1","C3-A2","C4-A1","C4-A2","C5-A1","C5-A2"];
types   = ["-T","-D","-S"];
mi18 = zeros(1,18); pc18 = zeros(1,18); k = 0;
for c = configs
  for t = types
    k = k+1;
    mi18(k) = find(contains(tx,c) & contains(tx,"-MI-") & contains(tx,t),1);
    pc18(k) = find(contains(tx,c) & contains(tx,"-PC-") & contains(tx,t),1);
  end
end

A = size(all_cellsnew_hm1,1);  % # animals

% --- Per-animal cross-lens Jaccard and right-tailed Fisher p (for enrichment)
Jacc = nan(A,18,18);      % Jaccard per animal
p_an = nan(A,18,18);      % per-animal Fisher p (right tail)

for a = 1:A
  for i = 1:18          % MI (rows)
    Ai = logical(all_cellsnew_hm1{a, mi18(i)});
    for j = 1:18        % PC (cols)
      Bj = logical(all_cellsnew_hm1{a, pc18(j)});
      kAB = sum(Ai & Bj);
      ao  = sum(Ai & ~Bj);
      bo  = sum(~Ai & Bj);
      nn  = sum(~Ai & ~Bj);

      uAB = kAB + ao + bo;
      if uAB==0, Jacc(a,i,j) = NaN; else, Jacc(a,i,j) = kAB / uAB; end

      try
        [~, p_an(a,i,j)] = fishertest([kAB ao; bo nn], 'Tail','right');
      catch
        N = numel(Ai); K = sum(Ai); M = sum(Bj);
        p_an(a,i,j) = 1 - hygecdf(kAB-1, N, K, M); % fallback
      end
    end
  end
end

% --- Mean Jaccard across animals
Jmean = squeeze(nanmean(Jacc,1));            % 18x18

% --- Stouffer combine per-pair p across animals (right-tailed)
p_comb = nan(18,18);
for i = 1:18
  for j = 1:18
    pv = squeeze(p_an(:,i,j));
    ok = pv>0 & pv<1 & ~isnan(pv);
    if any(ok)
      z  = norminv(1 - pv(ok));              % right-tailed z
      Z  = sum(z) / sqrt(sum(ok));
      p_comb(i,j) = 1 - normcdf(Z);          % combined p (right)
    end
  end
end
q = mafdr(p_comb(:),'BHFDR',true); q = reshape(q,size(p_comb));

% --- Order rows/cols by your earlier dendrograms (Step 1 and Step 3)
if exist('ord_mi','var'), rord = ord_mi; else, rord = 1:18; end
if exist('ord_pc','var'), cord = ord_pc; else, cord = 1:18; end

Jplot = Jmean(rord, cord);
qplot = q(rord, cord);

% --- Labels
labs_row = tx(mi18(rord));    % MI
labs_col = tx(pc18(cord));    % PC

% --- Plot
figure(5);
imagesc(Jplot,[0 1]); axis square; colormap(parula); colorbar
set(gca,'XTick',1:18,'YTick',1:18,'XTickLabel',labs_col,'YTickLabel',labs_row);
xtickangle(45);
title('Cross-lens overlap: MI (rows) vs PC (cols) — mean Jaccard'); hold on

% Mark significant enriched overlaps (q<=0.05)
[ri,ci] = find(qplot<=0.05);
plot(ci,ri,'k.','MarkerSize',10);
hold off
%%
% MI → PC loop (patched)
kTop = 3; thrQ = 0.05;

fprintf('\n=== MI → PC top matches ===\n');
for i = 1:18
    s  = Jmean(i,:);                % MI row vs all PC
    qs = q(i,:);                    % corresponding q-values
    valid = ~isnan(s);
    [~,ix] = sort(s(valid), 'descend');
    cols = find(valid);
    cols = cols(ix);

    fprintf('\nMI: %s\n', labs_row(i));
    shown = 0;
    for jj = 1:numel(cols)          % <— robust scalar loop
        j   = cols(jj);             %     j is a scalar
        star = '';
        val  = qs(j);               % scalar q-value
        if ~isnan(val) && val <= thrQ
            star = '*';
        end
        fprintf('  %2d) %-18s  J=%.3f  q=%g %s\n', ...
                shown+1, labs_col(j), s(j), val, star);
        shown = shown + 1;
        if shown >= kTop, break; end
    end
end
