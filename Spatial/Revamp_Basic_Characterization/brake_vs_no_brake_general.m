%% for VENN Diagrams with two circles
while 1
    %%
    % for different VENN Diagrams change good_FRV = all_exc(:,[1 3]) or
    % all_exc(:,[2 4]) or all_exc(:,[5 6]) same thing for all_inh
%     good_FRV = [cell_list_op(all_gV(:,[1:3 7]),[],'or',1) cell_list_op(all_gV(:,[4:6 8 9]),[],'or',1)]; 
%     good_FRV = [cell_list_op(all_gV(:,[1:3]),[],'or',1) cell_list_op(all_gV(:,[4:6]),[],'or',1)]; 
    good_FRV = all_gV;%[cell_list_op(all_gV(:,[1:3]),[],'or',1) cell_list_op(all_gV(:,[4:6]),[],'or',1)]; 
    cell_any = descriptiveStatistics(find_percent(cell_list_op(good_FRV,[],'or',1)));
    tcolors = mData.colors;
    hf = get_figure(6,[10 7 1.5 1]);
%     good_FRV = good_FR_bnb; %good_FRV = all_inh(:,[5 6]);
    [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(good_FRV,0.5,0.05);
    
    intcells = squeeze(all_CI_mat(1,2,:));    bcells = squeeze(all_CI_mat(1,1,:)); nbcells = squeeze(all_CI_mat(2,2,:));
    bcellsu = squeeze(uni(1,2,:)); nbcellsu = squeeze(uni(2,1,:));
    aVar_P = [bcellsu intcells nbcellsu];
    [within,dvn,xlabels,withinD] = make_within_table({'PT'},[3]); withinD3 = withinD;
    dataT = make_between_table({aVar_P},dvn);
    ra_pt = RMA(dataT,within,{0.05,{'hsd','bonferroni'}}); ra_pt.ranova
    print_for_manuscript(ra_pt); 
    
    
    p_gFR_C = [];
    for rr = 1:size(good_FRV,2)
        for cc = 1:size(good_FRV,2)
            gFR_C = cell_list_op(good_FRV(:,rr),good_FRV(:,cc),'and',1);
            temp_gFR_C = 100*exec_fun_on_cell_mat(gFR_C,'sum')./exec_fun_on_cell_mat(gFR_C,'length');
            [m_p_gFR_C(rr,cc),sem_p_gFR_C(rr,cc)] = findMeanAndStandardError(temp_gFR_C);
            p_gFR_C(rr,cc) = mean(100*exec_fun_on_cell_mat(gFR_C,'sum')./exec_fun_on_cell_mat(gFR_C,'length'));
        end
    end
    
    gFR_C_all = cell_list_op(good_FRV(:,1),good_FRV(:,2),'and',1);%    gFR_C_all = cell_list_op(gFR_C_all,good_FRV(:,3),'and',1);
    p_gFR_C_all = mean(100*exec_fun_on_cell_mat(gFR_C_all,'sum')./exec_fun_on_cell_mat(gFR_C_all,'length'));
    [mp_gFR_C_all,semp_gFR_C_all] = findMeanAndStandardError(100*exec_fun_on_cell_mat(gFR_C_all,'sum')./exec_fun_on_cell_mat(gFR_C_all,'length'));
    AVenn = [p_gFR_C(1,1) p_gFR_C(2,2)]; IVenn = [p_gFR_C(1,2)];
    [HVenn] = venn(AVenn,IVenn,'ErrMinMode','None');
%     AVenn = [p_gFR_C(1,1) p_gFR_C(2,2) p_gFR_C(3,3)]; IVenn = [p_gFR_C(1,2) p_gFR_C(1,3) p_gFR_C(2,3) p_gFR_C_all];
    [HVenn] = venn(AVenn,IVenn,'ErrMinMode','None');
    format_axes(gca);
    axis equal; axis off;
    changePosition(gca,[0.0 -0.0 -0.05 -0.05]); %put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -5 0]});
    pmchar=char(177);
    set(HVenn(1),'FaceColor',tcolors{3},'FaceAlpha',0.75); set(HVenn(2),'FaceColor',tcolors{4},'FaceAlpha',0.75);
    ylims = ylim;
    btxt = (sprintf('%0.1f%c%0.1f%%',m_p_gFR_C(1,1),pmchar,sem_p_gFR_C(1,1)));
    nbtxt = (sprintf('%0.1f%c%0.1f%%',m_p_gFR_C(2,2),pmchar,sem_p_gFR_C(2,2)));
    inttxt = (sprintf('%0.1f%c%0.1f%%',m_p_gFR_C(1,2),pmchar,sem_p_gFR_C(1,2)));
    clc
    disp(btxt);
    disp(inttxt);
    disp(nbtxt);
    save_pdf(hf,mData.pdf_folder,sprintf('conjunctive_cells_venn_diagram_%s.pdf',''),600);
    
    %%
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 1],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.29],'widthHeightAdjustment',...
        [10 -410]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 1 1.25]);
    MY = 60; ysp = 4; sigLinesStartYFactor = 1.5; mY = 0; % responsive cells
    stp = 0.43; widths = [0.4 0.4 0.4 0.4 0.4 0.4]+0.1; gap = 0.09;
    set(ff.hf,'Units','inches');
    for aii = 1:length(ff.h_axes)
        sel_ax = ff.h_axes(1,aii);
        set(sel_ax,'Units','inches');
        axesPosd = get(sel_ax,'Position');
        if aii == 1
            rt = stp; 
%             ylabel(sel_ax,{'Estimated','Marginal Means'});
        else
            rt = rt + widths(aii-1) + gap; set(sel_ax,'yticklabels',[]);
        end
        set(sel_ax,'Position',[rt axesPosd(2) widths(aii) axesPosd(4)]);
        ylim(sel_ax,[mY MY]);
    end
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_pt,{'PT','bonferroni'},[1.5 1 1]); 
	xdata = make_xdata([3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'Comp1','Conj','Comp2'};set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 25 50]); xtickangle(45);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs_3_factors.pdf',ntrials),600);
    %%
    break;
end

%% for VENN Diagrams with three circles
while 1
    %%
    % for different VENN Diagrams change good_FRV = all_exc(:,[1 3]) or
    % all_exc(:,[2 4]) or all_exc(:,[5 6]) same thing for all_inh
%     good_FRV = [cell_list_op(all_gV(:,[1:2 5]),[],'or',1) cell_list_op(all_gV(:,[3:4 6]),[],'or',1)]; 
    good_FRV = all_gV;
    cell_any = descriptiveStatistics(find_percent(cell_list_op(good_FRV,[],'or',1)));
    tcolors = mData.colors;
    hf = get_figure(6,[10 7 1.5 1]);
%     good_FRV = good_FR_bnb; %good_FRV = all_inh(:,[5 6]);
    [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(good_FRV,0.5,0.05);
    p_gFR_C = [];
    for rr = 1:size(good_FRV,2)
        for cc = 1:size(good_FRV,2)
            gFR_C = cell_list_op(good_FRV(:,rr),good_FRV(:,cc),'and',1);
            temp_gFR_C = 100*exec_fun_on_cell_mat(gFR_C,'sum')./exec_fun_on_cell_mat(gFR_C,'length');
            [m_p_gFR_C(rr,cc),sem_p_gFR_C(rr,cc)] = findMeanAndStandardError(temp_gFR_C);
            p_gFR_C(rr,cc) = mean(100*exec_fun_on_cell_mat(gFR_C,'sum')./exec_fun_on_cell_mat(gFR_C,'length'));
        end
    end
    
    gFR_C_all = cell_list_op(good_FRV(:,1),good_FRV(:,2),'and',1);%    gFR_C_all = cell_list_op(gFR_C_all,good_FRV(:,3),'and',1);
    p_gFR_C_all = mean(100*exec_fun_on_cell_mat(gFR_C_all,'sum')./exec_fun_on_cell_mat(gFR_C_all,'length'));
    [mp_gFR_C_all,semp_gFR_C_all] = findMeanAndStandardError(100*exec_fun_on_cell_mat(gFR_C_all,'sum')./exec_fun_on_cell_mat(gFR_C_all,'length'));
%     AVenn = [p_gFR_C(1,1) p_gFR_C(2,2)]; IVenn = [p_gFR_C(1,2)];
%     [HVenn] = venn(AVenn,IVenn,'ErrMinMode','None');
    AVenn = [p_gFR_C(1,1) p_gFR_C(2,2) p_gFR_C(3,3)]; IVenn = [p_gFR_C(1,2) p_gFR_C(1,3) p_gFR_C(2,3) p_gFR_C_all];
    [HVenn] = venn(AVenn,IVenn,'ErrMinMode','None');
    format_axes(gca);
    axis equal; axis off;
    changePosition(gca,[0.0 -0.0 -0.05 -0.05]); %put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -5 0]});
    pmchar=char(177);
    set(HVenn(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); set(HVenn(2),'FaceColor',tcolors{2},'FaceAlpha',0.75);
    ylims = ylim;
    btxt = (sprintf('%0.1f%c%0.1f%%',m_p_gFR_C(1,1),pmchar,sem_p_gFR_C(1,1)));
    nbtxt = (sprintf('%0.1f%c%0.1f%%',m_p_gFR_C(2,2),pmchar,sem_p_gFR_C(2,2)));
    inttxt = (sprintf('%0.1f%c%0.1f%%',m_p_gFR_C(1,2),pmchar,sem_p_gFR_C(1,2)));
    clc
    disp(btxt);
    disp(inttxt);
    disp(nbtxt);
    save_pdf(hf,mData.pdf_folder,sprintf('conjunctive_cells_venn_diagram_%s.pdf',''),600);
    %%
    break;
end

%% for Heap Map
while 1
    %% for heat maps
    figdim = 3;
    hf = get_figure(5,[9 2 figdim figdim]);
            
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
    
%     mOI = mUni;
%     mOI = mUni1; semOI = semUni1;
%     mOI = mUni2; semOI = semUni2;
    
    sz = size(mOI,1);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN; semOI(mask1 == 1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 

    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;

    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(75);
    for rr = 1:size(mCI,1)
        for cc = 1:size(mCI,1)
            if rr == cc
                text(rr-0.25,cc,sprintf('%.0f',mCI(rr,cc)),'FontSize',5);
            end
        end
    end
%     plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 150.5],'w','linewidth',0.1); 
%     plot([0 150.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
    set(gca,'Ydir','normal');ytickangle(15);
    box on
    changePosition(gca,[0.0 0 0.0 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.07 0.05 0.07 0.05]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean_%d.pdf',ntrials,sh),600);
    
    figdim = 1;
    hf = get_figure(6,[5 2 figdim+0.1 figdim]);
    im1 = imagesc(semOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',[],'yticklabels',[],'Ydir','reverse'); xtickangle(75);
    changePosition(gca,[0.0 0 -0.05 0]);
    set(gca,'Ydir','normal');
    box on
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_sem_%d.pdf',ntrials,sh),600);
   %%
   ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 3],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.35],'widthHeightAdjustment',[10 -450]);
    set(gcf,'color','w');    set(gcf,'Position',[10 3 6.95 1.25]);    ylims = [0 1];
    stp = 0.25; widths = [2 2 2 0.4 0.4 0.4]+0.1; gap = 0.13; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
    stp = 0.25; widths = [2 2 2 0.4 0.4 0.4]+0.1; gap = 0.17; adjust_axes(ff,ylims,stp,widths,gap,{'Euclidean Distance'});
    %
    ahc_col_th = 0.7;
    mOI1 = mCI;
%     mOI1 = mUni1;
    mOI1(mask1==1) = NaN; 
    Di = pdist(mOI1,@naneucdist);
    tree = linkage(Di,'average');
    [c,d] = cophenet(tree,Di); r = corr(Di',d','type','spearman');
    leafOrder = optimalleaforder(tree,Di);
    hf = figure(100000000);
%     leafOrder1 = leafOrder([1:3 10:12 4:9]);
%     leafOrder1 = circshift(leafOrder,3);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);ylims = ylim; xlims = xlim;
    set(gca,'xtick',[],'ytick',[]);
    set(gcf,'units','inches'); set(gcf,'Position',[5 2 0.9 0.5])
    %
    close(hf);
%     
    
    axes(ff.h_axes(1,1));
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);
    set(H,'linewidth',0.5);
    set(gca,'xticklabels',txl(leafOrder));xtickangle(45);
    format_axes(gca);
    hx = ylabel({'Euc. Dist.'});changePosition(hx,[0 0 0]);
%     xlim([xlims(1)+0.5 xlims(2)-0.5]);
    changePosition(gca,[0.0 0.0 0.07 0.05]);
    text(0.5,ylims(2)+0,sprintf('CC = %.2f',c),'FontSize',6);
%     set_axes_top_text(ff.hf,ff.h_axes(1),sprintf('Cophenetic Correlation = %.2f',c));

%     mOI1 = mCI;
    mOI1 = mUni1;
    mOI1(mask1==1) = NaN; 
    Di = pdist(mOI1,@naneucdist);
    tree = linkage(Di,'average');
    [c,d] = cophenet(tree,Di); r = corr(Di',d','type','spearman');
    leafOrder = optimalleaforder(tree,Di);
    hf = figure(100000000);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);ylims = ylim; xlims = xlim;
    close(hf);
%     
    
    axes(ff.h_axes(1,2));
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);
    set(H,'linewidth',0.5);
    set(gca,'xticklabels',txl(leafOrder));xtickangle(45);
    format_axes(gca);
%     xlim([xlims(1)+0.5 xlims(2)-0.5]);
    changePosition(gca,[0.0 0.0 0.07 0.05]);
    text(0.5,ylims(2)+1,sprintf('CC = %.2f',c),'FontSize',6);

    
%     mOI1 = mCI;
    mOI1 = mUni2;
    mOI1(mask1==1) = NaN; 
    Di = pdist(mOI1,@naneucdist);
    tree = linkage(Di,'average');
    [c,d] = cophenet(tree,Di); r = corr(Di',d','type','spearman');
    leafOrder = optimalleaforder(tree,Di);
    hf = figure(100000000);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);ylims = ylim; xlims = xlim;
    close(hf);
%     
    
    axes(ff.h_axes(1,3));
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)),'Reorder',leafOrder);
    set(H,'linewidth',0.5);
    set(gca,'xticklabels',txl(leafOrder));xtickangle(45);
    format_axes(gca);
%     xlim([xlims(1)+0.5 xlims(2)-0.5]);
    changePosition(gca,[0.0 0.0 0.07 0.05]);
    text(0.5,ylims(2)+1,sprintf('CC = %.2f',c),'FontSize',6);
%     delete(ff.h_axes(1,3));
    save_pdf(ff.hf,mData.pdf_folder,sprintf('OI_Map_cluster_%d.pdf',sh),600);
    
    %%
    break;
end

%% Agglomerative hierarchical clustering
while 1
    mOI1 = mCI; mOI1(mask1==1) = NaN; DiC = pdist(mOI1,@naneucdist); 
    mOI1 = mUni1; mOI1(mask1==1) = NaN; DiC1 = pdist(mOI1,@naneucdist); 
    mOI1 = mUni2; mOI1(mask1==1) = NaN; DiC2 = pdist(mOI1,@naneucdist);
    Di = mean([DiC1;DiC2]);
%     tree = linkage(mOI1,'average','euclidean');
    tree = linkage(Di,'average');
    hf = figure(1000);clf; set(hf,'Units','Inches');
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold',ahc_col_th*max(tree(:,3)));
    [c,d] = cophenet(tree,Di); r = corr(Di',d','type','spearman');
    hf = gcf;
    set(hf,'Position',[7 3 3.5 1.75]);
    set(hf,'Position',[7 3 2.3 1.25]);
%     set(hf,'Position',[7 3 2.4 1.15]);
    set(H,'linewidth',0.5);
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel({'Eucledian Distance'});changePosition(hx,[0 0 0]);
    changePosition(gca,[0.0 0.0 0.07 0.05]); ylims = ylim;
    text(0.5,ylims(2),sprintf('CC = %.2f',c),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);
    %%
    break;
end


%% for bar graphs type_by_cond and cond and type ... all together (for both responsive cells and response fidelity)
while 1
   %%
   ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 3],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.08 0.13],'widthHeightAdjustment',...
        [10 -60]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 1.6 1]);
    
    MY = 90; ysp = 5;
    MY = 10; ysp = 1;
    axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_Type','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([2 2],[1 1.5]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.colors([3 7]),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    maxY = maxY + 1;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) MY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Exc','Inh'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.07 0.07 0.01 -0.3]); 
    put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    put_axes_labels(gca,{[],[0 0 0]},{{'Trials (%)'},[0 0 0]});
    
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,2));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    maxY = maxY + 1;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) MY]); format_axes(gca);
    xticks = xdata; xticklabels = {'B','NB'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[]); xtickangle(45)
    changePosition(gca,[0.17 0.07 -0.2 -0.3]); %put_axes_labels(gca,{[],[0 0 0]},{{'Responsive','Cells (%)'},[0 0 0]});
%     text(2,MY+3,'Pooled','FontSize',6);
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,3));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Type','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([2],[1 1.5]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.colors([3 7]),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) MY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Exc','Inh'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[]); xtickangle(45)
    changePosition(gca,[0.05 0.07 -0.2 -0.3]); %put_axes_labels(gca,{[],[0 0 0]},{{'Responsive','Cells (%)'},[0 0 0]});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all3.pdf',ntrials),600);
    
    %%
    break;
end


%% for bar graphs type_by_cond and cond and type ... all together (for zMI)
while 1
   %%
   ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 3],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.08 0.13],'widthHeightAdjustment',...
        [10 -60]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 1.6 1]);
    
    MY = 2; ysp = 0.3;
    axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_Type','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([2 2],[1 1.5]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.colors([3 7]),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    maxY = maxY + 1;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) MY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Exc','Inh'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.07 0.07 0.01 -0.3]); 
    put_axes_labels(gca,{[],[0 0 0]},{{'zMI'},[0 0 0]});
    
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,2));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    maxY = maxY + 1;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) MY]); format_axes(gca);
    xticks = xdata; xticklabels = {'B','NB'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[]); xtickangle(45)
    changePosition(gca,[0.17 0.07 -0.2 -0.3]); %put_axes_labels(gca,{[],[0 0 0]},{{'Responsive','Cells (%)'},[0 0 0]});
%     text(2,MY+3,'Pooled','FontSize',6);
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,3));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Type','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([2],[1 1.5]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.colors([3 7]),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) MY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Exc','Inh'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[]); xtickangle(45)
    changePosition(gca,[0.05 0.07 -0.2 -0.3]); %put_axes_labels(gca,{[],[0 0 0]},{{'Responsive','Cells (%)'},[0 0 0]});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all3.pdf',ntrials),600);
    
    %%
    break;
end


%% for bar graphs three factors (all possible)
while 1
   %%
   clc
   ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 7],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.06 0.23],'widthHeightAdjustment',...
        [10 -350]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 4.85 1]);
%     MY = 2.55; ysp = 0.2; %zMI
%     MY = 1; ysp = 0.04; mY = 0.4; %HaFD
    MY = 15; ysp = 1.5; mY = 0; % responsive cells
%     MY = 100; ysp = 6; mY = 0; % response fidelity
%     MY = 1; ysp = 0.06; mY = 0; % zMINaN
    stp = 0.35; widths = [1.35 0.25 0.4 0.2 0.75 0.45 0.75]; gap = 0.05;
    axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT_ET_Cond','hsd'},[1.5 1 1]);
    h(h==1) = 0; xdata = make_xdata([2 2 2 2 2 2],[1 1.5]);    tcolors = repmat(mData.colors([7 6]),1,6);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca); xticks = [1.5 4 6.5 9 11.5 14]; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)

%     rectangle(gca,'Position',[1 12 1 1],'FaceColor',tcolors{1},'EdgeColor',tcolors{1}); text(2.5,12.7,'Exc','FontSize',6);
%     rectangle(gca,'Position',[5.5 12 1 1],'FaceColor',tcolors{2},'EdgeColor',tcolors{2}); text(7,12.7,'Inh','FontSize',6);
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,2));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'B','NB'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,3));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'ET','hsd'},[1.5 1 1]);
	xdata = make_xdata([3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,4));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);    tcolors = repmat(mData.colors([7 6]),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);   xticks = xdata; xticklabels = {'Exc','Inh'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,5));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_ET','hsd'},[1.5 1 1]);
    if ~hs_flag(1) h(h==1) = 0; end
	xdata = make_xdata([3 3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,6));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_CT','hsd'},[1.5 1 1]);
    if ~hs_flag(2) h(h==1) = 0; end
	xdata = make_xdata([2 2],[1 1.5]);    tcolors = repmat(mData.colors([7 6]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = [1.5 4]; xticklabels = {'B','NB'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,7));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'ET_by_CT','hsd'},[1.5 1 1]);
    if ~hs_flag(3) h(h==1) = 0; end
	xdata = make_xdata([2 2 2],[1 1.5]);    tcolors = repmat(mData.colors([7 6]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = [1.5 4 6.5]; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    %+++++++++++++++++++++++++
%     MY = 3;
    set(ff.hf,'Units','inches');
    for aii = 1:7
        sel_ax = ff.h_axes(1,aii);
        set(sel_ax,'Units','inches');
        axesPosd = get(sel_ax,'Position');
        if aii == 1
            rt = stp; 
            ylabel(sel_ax,'Cells (%)');
%             ylabel(sel_ax,'zMI');
%             ylabel(sel_ax,'Trials (%)');
%             ylabel(sel_ax,'Hausdorff FD');
        else
            rt = rt + widths(aii-1) + gap; set(sel_ax,'ytick',[]);
        end
        set(sel_ax,'Position',[rt axesPosd(2) widths(aii) axesPosd(4)]);
        ylim(sel_ax,[mY MY]);
    end

    save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs_3_factors.pdf',ntrials),600);
    
    %%
    break;
end

%% for bar graphs three factors (selected)
while 1
   %%
   clc
   ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 5],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.06 0.23],'widthHeightAdjustment',...
        [10 -350]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 4.85 1.3]);
%     MY = 2.55; ysp = 0.2; mY = -0.5; %zMI
%     MY = 1; ysp = 0.04; mY = 0.4; %HaFD
    MY = 15; ysp = 1.5; mY = 0; % responsive cells
%     MY = 100; ysp = 6; mY = 0; % response fidelity
%     MY = 1; ysp = 0.06; mY = 0; % zMINaN
    stp = 0.35; widths = [1.85 0.35 0.6 0.35 0.9]; gap = 0.07;
    set(ff.hf,'Units','inches');
    for aii = 1:length(ff.h_axes)
        sel_ax = ff.h_axes(1,aii);
        set(sel_ax,'Units','inches');
        axesPosd = get(sel_ax,'Position');
        if aii == 1
            rt = stp; 
            ylabel(sel_ax,'Estimated Marginal Means');
        else
            rt = rt + widths(aii-1) + gap; set(sel_ax,'yticklabels',[]);
        end
        set(sel_ax,'Position',[rt axesPosd(2) widths(aii) axesPosd(4)]);
        ylim(sel_ax,[mY MY]);
    end
    %++++++++
    axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_ET_CT','hsd'},[1.5 1 1]);
    h(h==1) = 0; xdata = make_xdata([2 2 2 2 2 2],[1 1.5]);    tcolors = repmat(mData.colors([7 6]),6,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca); xticks = [1.5 4 6.5 9 11.5 14]; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)

    rectangle(gca,'Position',[1 12 1 1],'FaceColor',tcolors{1},'EdgeColor',tcolors{1}); text(2.5,12.7,'Exc','FontSize',6);
    rectangle(gca,'Position',[5.5 12 1 1],'FaceColor',tcolors{2},'EdgeColor',tcolors{2}); text(7,12.7,'Inh','FontSize',6);
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,2));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'B','NB'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,3));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'ET','hsd'},[1.5 1 1]);
	xdata = make_xdata([3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,4));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);    tcolors = repmat(mData.colors([7 6]),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);   xticks = xdata; xticklabels = {'Exc','Inh'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,5));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'ET_by_CT','hsd'},[1.5 1 1]);
    if ~hs_flag(3) h(h==1) = 0; end
	xdata = make_xdata([2 2 2],[1 1.5]);    tcolors = repmat(mData.colors([7 6]),4,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigLineWidth',0.25,'BaseValue',0.01,'xdata',xdata,'sigFontSize',6,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = [1.5 4 6.5]; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)

    save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs_3_factors.pdf',ntrials),600);
    
    %%
    break;
end


%% for bar graphs three factors (for Figure 2)
while 1
   %%
   clc
   ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 5],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.06 0.23],'widthHeightAdjustment',...
        [10 -350]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 4.85 1.3]);
%     MY = 2.55; ysp = 0.2; mY = -0.5; %zMI
%     MY = 1; ysp = 0.04; mY = 0.4; %HaFD
    MY = 15; ysp = 1.5; mY = 0; % responsive cells
%     MY = 100; ysp = 6; mY = 0; % response fidelity
%     MY = 1; ysp = 0.06; mY = 0; % zMINaN
    stp = 0.35; widths = [1.85 0.35 0.6 0.35 0.9]; gap = 0.07;
    set(ff.hf,'Units','inches');
    for aii = 1:length(ff.h_axes)
        sel_ax = ff.h_axes(1,aii);
        set(sel_ax,'Units','inches');
        axesPosd = get(sel_ax,'Position');
        if aii == 1
            rt = stp; 
            ylabel(sel_ax,'Estimated Marginal Means');
        else
            rt = rt + widths(aii-1) + gap; set(sel_ax,'yticklabels',[]);
        end
        set(sel_ax,'Position',[rt axesPosd(2) widths(aii) axesPosd(4)]);
        ylim(sel_ax,[mY MY]);
    end
    %++++++++
    axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT_ET_by_Cond','hsd'},[1.5 1 1]);
    h(h==1) = 0; xdata = make_xdata([2 2 2 2 2 2],[1 1.5]);    tcolors = repmat(mData.colors([1 2]),6,1);
    combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca); xticks = [1.5 4 6.5 9 11.5 14]; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)

    rectangle(gca,'Position',[1 12 1 1],'FaceColor',tcolors{1},'EdgeColor',tcolors{1}); text(2.5,12.7,'Brake','FontSize',6);
    rectangle(gca,'Position',[5.5 12 1 1],'FaceColor',tcolors{2},'EdgeColor',tcolors{2}); text(7,12.7,'No-Brake','FontSize',6);
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,2));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'B','NB'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,3));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'ET','hsd'},[1.5 1 1]);
	xdata = make_xdata([3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,4));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);    tcolors = repmat(mData.colors([7 6]),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);   xticks = xdata; xticklabels = {'Exc','Inh'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,5));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'ET_by_CT','hsd'},[1.5 1 1]);
    if ~hs_flag(3) h(h==1) = 0; end
	xdata = make_xdata([2 2 2],[1 1.5]);    tcolors = repmat(mData.colors([7 6]),4,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigLineWidth',0.25,'BaseValue',0.01,'xdata',xdata,'sigFontSize',6,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = [1.5 4 6.5]; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)

    save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs_3_factors.pdf',ntrials),600);
    
    %%
    break;
end

%% Figure 2 Percentage of cells (Short Version)
while 1
    clc
    title_prefix = 'C'
    titles = {'1-Exc','2-Exc','3-Inh','4-Inh','5-AOn','6-AOn','7-AOff'};
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 7],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',...
        [10 -410]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 4.85 1.3]);
    MY = 12; ysp = 1.5; mY = 0; % responsive cells
    stp = 0.35; widths = [0.4 0.57 0.4 0.57 0.4 0.4 0.4]+0.1; gap = 0.09;
    set(ff.hf,'Units','inches');
    for aii = 1:length(ff.h_axes)
        sel_ax = ff.h_axes(1,aii);
        set(sel_ax,'Units','inches');
        axesPosd = get(sel_ax,'Position');
        if aii == 1
            rt = stp; 
            ylabel(sel_ax,{'Estimated','Marginal Means'});
        else
            rt = rt + widths(aii-1) + gap; set(sel_ax,'yticklabels',[]);
        end
        set(sel_ax,'Position',[rt axesPosd(2) widths(aii) axesPosd(4)]);
        ylim(sel_ax,[mY MY]);
    end
    
    %+++++++++++++++++++++++
    axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_exc,{'Cond','hsd'},[1.5 1 1]);
%     inds = [1;(find(diff(combs(:,1))==1)+1)]; combs = combs(inds,:); p = p(inds,:); h = h(inds,:);
	xdata = make_xdata([2],[1 1.5]);    tcolors = repmat(mData.colors([1 2]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp-0.3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'B','NB'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,2,{'Pooled'});
    %++++++++++++++++++++++++
    axes(ff.h_axes(1,2));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_exc,{'ET','hsd'},[1.5 1 1]);
	xdata = make_xdata([3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,3,{'Pooled'});
    set_axes_top_text(ff.hf,ff.h_axes(1:2),'Exc Cells');
    %+++++++++++++++++++++++
    axes(ff.h_axes(1,3));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_inh,{'Cond','hsd'},[1.5 1 1]);
%     inds = [1;(find(diff(combs(:,1))==1)+1)]; combs = combs(inds,:); p = p(inds,:); h = h(inds,:);
	xdata = make_xdata([2],[1 1.5]);    tcolors = repmat(mData.colors([1 2]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp-0.3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'B','NB'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,2,{'Pooled'});
    %++++++++++++++++++++++++
    axes(ff.h_axes(1,4));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_inh,{'ET','hsd'},[1.5 1 1]);
	xdata = make_xdata([3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,3,{'Pooled'});
    set_axes_top_text(ff.hf,ff.h_axes(3:4),'Inh Cells');
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,5));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_e1,{'Cond','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);    tcolors = repmat(mData.colors([1 2]),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);   xticks = xdata; xticklabels = {'B','NB'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,2,{'Pooled'});
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,6));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_e1,{'CT','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);    tcolors = repmat(mData.colors([7 6]),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);   xticks = xdata; xticklabels = {'Exc','Inh'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,2,{'Pooled'});
    set_axes_top_text(ff.hf,ff.h_axes(5:6),'AOn');
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,7));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_e2,{'Cond','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);    tcolors = repmat(mData.colors([1 2]),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);   xticks = xdata; xticklabels = {'B','NB'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,2,{'Pooled'});
    set_axes_top_text(ff.hf,ff.h_axes(7),'AOff');
    
%     for aii = 1:length(ff.h_axes)
%         sel_ax = ff.h_axes(1,aii);axes(sel_ax);
%         xlims = get(sel_ax,'xlim'); ylims = get(sel_ax,'ylim'); xlimshalf = xlims(1) + diff(xlims)/2; title_txt = sprintf('%s%s',title_prefix,titles{aii});
%         ht = text(xlimshalf,ylims(2)+ysp,title_txt,'FontSize',6); ht_ex = get(ht,'Extent'); ht_pos = get(ht,'Position'); 
%         new_x = xlims(1) + (xlims(2) - xlims(1) - ht_ex(3))/2; set(ht,'Position',[new_x ht_pos(2:3)]);
%     end
    
    save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs_3_factors.pdf',ntrials),600);
    %%
    
    break;
end

%% Figure 2 Percentage of cells (Detailed Version)
while 1
    clc
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 10],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',...
        [10 -410]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 4.85 1.3]);
    MY = 12; ysp = 1.5; mY = 0; % responsive cells
    stp = 0.35; widths = [0.35 0.52 0.35 0.52 0.35 0.35 0.35 0.35 0.35 0.35 0.35 0.35]-0.02; gap = 0.09;
    set(ff.hf,'Units','inches');
    for aii = 1:length(ff.h_axes)
        sel_ax = ff.h_axes(1,aii);
        set(sel_ax,'Units','inches');
        axesPosd = get(sel_ax,'Position');
        if aii == 1
            rt = stp; 
            ylabel(sel_ax,{'Estimated','Marginal Means'});
        else
            rt = rt + widths(aii-1) + gap; set(sel_ax,'yticklabels',[]);
        end
        set(sel_ax,'Position',[rt axesPosd(2) widths(aii) axesPosd(4)]);
        ylim(sel_ax,[mY MY]);
    end
    
    %+++++++++++++++++++++++
    axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_exc,{'Cond','hsd'},[1.5 1 1]);
%     inds = [1;(find(diff(combs(:,1))==1)+1)]; combs = combs(inds,:); p = p(inds,:); h = h(inds,:);
	xdata = make_xdata([2],[1 1.5]);    tcolors = repmat(mData.colors([1 2]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp-0.3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'B','NB'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,2,{'Pooled'});
    %++++++++++++++++++++++++
    axes(ff.h_axes(1,2));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_exc,{'ET','hsd'},[1.5 1 1]);
	xdata = make_xdata([3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,3,{'Pooled'});
    set_axes_top_text(ff.hf,ff.h_axes(1:2),'Exc Cells');
    %+++++++++++++++++++++++
    axes(ff.h_axes(1,3));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_inh,{'Cond','hsd'},[1.5 1 1]);
%     inds = [1;(find(diff(combs(:,1))==1)+1)]; combs = combs(inds,:); p = p(inds,:); h = h(inds,:);
	xdata = make_xdata([2],[1 1.5]);    tcolors = repmat(mData.colors([1 2]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp-0.3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'B','NB'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,2,{'Pooled'});
    %++++++++++++++++++++++++
    axes(ff.h_axes(1,4));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_inh,{'ET','hsd'},[1.5 1 1]);
	xdata = make_xdata([3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,3,{'Pooled'});
    set_axes_top_text(ff.hf,ff.h_axes(3:4),'Inh Cells');
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,5));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_e1,{'Cond','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);    tcolors = repmat(mData.colors([1 2]),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);   xticks = xdata; xticklabels = {'B','NB'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,2,{'Pooled'});
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,6));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_e1,{'CT','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);    tcolors = repmat(mData.colors([7 6]),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);   xticks = xdata; xticklabels = {'Exc','Inh'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,2,{'Pooled'});
    set_axes_top_text(ff.hf,ff.h_axes(5:6),'AOn');
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,7));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_e2,{'Cond','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);    tcolors = repmat(mData.colors([1 2]),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);   xticks = xdata; xticklabels = {'B','NB'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,2,{'Pooled'});
    set_axes_top_text(ff.hf,ff.h_axes(7),'AOff');
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,8));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_e2,{'CT','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);    tcolors = repmat(mData.colors([7 6]),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);   xticks = xdata; xticklabels = {'Exc','Inh'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,2,{'Pooled'});
    set_axes_top_text(ff.hf,ff.h_axes(7:8),'AOff');
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,9));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_e3,{'Cond','hsd'},[1.5 1 1]); h(h==1) = 0;
	xdata = make_xdata([2],[1 1.5]);    tcolors = repmat(mData.colors([1 2]),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);   xticks = xdata; xticklabels = {'B','NB'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,2,{'Pooled'});
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,10));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_e3,{'CT','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);    tcolors = repmat(mData.colors([7 6]),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);   xticks = xdata; xticklabels = {'Exc','Inh'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,2,{'Pooled'});
    set_axes_top_text(ff.hf,ff.h_axes(9:10),'Arb');
    
%     for aii = 1:length(ff.h_axes)
%         sel_ax = ff.h_axes(1,aii);axes(sel_ax);
%         xlims = get(sel_ax,'xlim'); ylims = get(sel_ax,'ylim'); xlimshalf = xlims(1) + diff(xlims)/2; title_txt = sprintf('%s%s',title_prefix,titles{aii});
%         ht = text(xlimshalf,ylims(2)+ysp,title_txt,'FontSize',6); ht_ex = get(ht,'Extent'); ht_pos = get(ht,'Position'); 
%         new_x = xlims(1) + (xlims(2) - xlims(1) - ht_ex(3))/2; set(ht,'Position',[new_x ht_pos(2:3)]);
%     end
    
    save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs_3_factors.pdf',ntrials),600);
    %%
    
    break;
end

%% Figure 2 Percentage of cells responsiveness (another vertsion ET-by-CT)
while 1
    clc
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',...
        [10 -410]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.75 1.25]);
    MY = 11; ysp = 0.75; sigLinesStartYFactor = 1.5; mY = 0; % responsive cells
    MY = 100; ysp = 5; sigLinesStartYFactor = 1.5; mY = 0; % responsive cells
    stp = 0.37; widths = [0.4 0.4 0.4 0.4]+0.1; gap = 0.1;
    set(ff.hf,'Units','inches');
    for aii = 1:length(ff.h_axes)
        sel_ax = ff.h_axes(1,aii);
        set(sel_ax,'Units','inches');
        axesPosd = get(sel_ax,'Position');
        if aii == 1
            rt = stp; 
            ylabel(sel_ax,{'Estimated','Marginal Means'});
        else
            rt = rt + widths(aii-1) + gap; set(sel_ax,'yticklabels',[]);
        end
        set(sel_ax,'Position',[rt axesPosd(2) widths(aii) axesPosd(4)]);
        ylim(sel_ax,[mY MY]);
    end
    
    
    %++++++++++++++++++++++++
    axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_cond1_exc,{'ET','bonferroni'},[1.5 1 1]); 
	xdata = make_xdata([3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,3,{'Exc'});
    %+++++++++++++++++++++++
     axes(ff.h_axes(1,2));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_cond1_inh,{'ET','bonferroni'},[1.5 1 1]); 
	xdata = make_xdata([3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,3,{'Inh'});
    set_axes_top_text(ff.hf,ff.h_axes(1:2),'Brake');
    %++++++++++++++++++++++++
    axes(ff.h_axes(1,3));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_cond2_exc,{'ET','bonferroni'},[1.5 1 1]); 
	xdata = make_xdata([3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,3,{'Exc'});
    %+++++++++++++++++++++++
     axes(ff.h_axes(1,4));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_cond2_inh,{'ET','bonferroni'},[1.5 1 1]); 
	xdata = make_xdata([3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,3,{'Inh'});
    set_axes_top_text(ff.hf,ff.h_axes(3:4),'No-Brake');
    
    
    
   
    save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs_3_factors.pdf',ntrials),600);
    %%
    
    break;
end

%% Figure 2 Percentage of Trials response fidelity
while 1
    clc
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 2],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',...
        [10 -410]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.7 1.25]);
    MY = 100; ysp = 7; mY = 0; % responsive cells
    MY = 2.25; ysp = 0.15; sigLinesStartYFactor = 0.1; mY = -0.5; % zMI
    stp = 0.45; widths = [0.95 0.95 0.57 0.4 0.4 0.4]+0.1; gap = 0.09;
    set(ff.hf,'Units','inches');
    for aii = 1:length(ff.h_axes)
        sel_ax = ff.h_axes(1,aii);
        set(sel_ax,'Units','inches');
        axesPosd = get(sel_ax,'Position');
        if aii == 1
            rt = stp; 
            ylabel(sel_ax,{'Estimated','Marginal Means'});
        else
            rt = rt + widths(aii-1) + gap; set(sel_ax,'yticklabels',[]);
        end
        set(sel_ax,'Position',[rt axesPosd(2) widths(aii) axesPosd(4)]);
        ylim(sel_ax,[mY MY]);
        if aii > 3
            delete(ff.h_axes(1,aii));
        end
    end
    
    
    %++++++++++++++++++++++++
    axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_cond1,{'CT_by_ET','bonferroni'},[1.5 1 1]);
	xdata = make_xdata([3 3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,3,{'Exc','Inh'});
    set_axes_top_text(ff.hf,ff.h_axes(1),'Brake');

     %++++++++++++++++++++++++
    axes(ff.h_axes(1,2));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_cond2,{'CT_by_ET','bonferroni'},[1.5 1 1]);
	xdata = make_xdata([3 3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,3,{'Exc','Inh'});
    set_axes_top_text(ff.hf,ff.h_axes(2),'No-Brake');
 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs_3_factors.pdf',ntrials),600);
    %%
    
    break;
end

%% Figure 2 MI
while 1
    clc
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 1],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',...
        [10 -410]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 1.5 1.3]);
    MY = 3; ysp = 0.2; mY = -0.5; % responsive cells
    stp = 0.35; widths = [0.95 1.2 0.57 0.4 0.4 0.4]+0.1; gap = 0.09;
    set(ff.hf,'Units','inches');
    for aii = 1:length(ff.h_axes)
        sel_ax = ff.h_axes(1,aii);
        set(sel_ax,'Units','inches');
        axesPosd = get(sel_ax,'Position');
        if aii == 1
            rt = stp; 
            ylabel(sel_ax,{'Estimated','Marginal Means'});
        else
            rt = rt + widths(aii-1) + gap; set(sel_ax,'yticklabels',[]);
        end
        set(sel_ax,'Position',[rt axesPosd(2) widths(aii) axesPosd(4)]);
        ylim(sel_ax,[mY MY]);
    end

    %++++++++++++++++++++++++
    axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT_by_ET','hsd'},[1.5 1 1]);
	xdata = make_xdata([3 3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,3,{'Exc','Inh'});
%     set_axes_top_text(ff.hf,ff.h_axes(2),'Pooled across Brake and No-Brake');
 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs_3_factors.pdf',ntrials),600);
    %%
    
    break;
end

%% Complementation and Conjunction
while 1
    clc
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',...
        [10 -410]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 3.2 1.25]);
    MY = 15; ysp = 0.75; sigLinesStartYFactor = 1.5; mY = 0; % responsive cells
    stp = 0.37; widths = [0.4 0.4 0.4 0.4 0.4 0.4]+0.1; gap = 0.09;
    set(ff.hf,'Units','inches');
    for aii = 1:length(ff.h_axes)
        sel_ax = ff.h_axes(1,aii);
        set(sel_ax,'Units','inches');
        axesPosd = get(sel_ax,'Position');
        if aii == 1
            rt = stp; 
            ylabel(sel_ax,{'Estimated','Marginal Means'});
        else
            rt = rt + widths(aii-1) + gap; set(sel_ax,'yticklabels',[]);
        end
        set(sel_ax,'Position',[rt axesPosd(2) widths(aii) axesPosd(4)]);
        ylim(sel_ax,[mY MY]);
    end

    %++++++++++++++++++++++++
    axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_p1_exc,{'ET','bonferroni'},[1.5 1 1]); 
	xdata = make_xdata([3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,3,{'Exc'});
    %+++++++++++++++++++++++
     axes(ff.h_axes(1,2));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_p1_inh,{'ET','bonferroni'},[1.5 1 1]); 
	xdata = make_xdata([3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,3,{'Inh'});
    set_axes_top_text(ff.hf,ff.h_axes(1:2),'Comp1 (Brake)');
    
    %++++++++++++++++++++++++
    axes(ff.h_axes(1,3));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_p2,{'ET','hsd'},[1.5 1 1]); 
	xdata = make_xdata([3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,[],[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,3,{'Pooled'});
    set_axes_top_text(ff.hf,ff.h_axes(3),'Conjunctive',[-0.01 0 0.2 0]);
    %+++++++++++++++++++++++
     axes(ff.h_axes(1,4));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_p3,{'ET','hsd'},[1.5 1 1]); 
	xdata = make_xdata([3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,[],[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,3,{'Pooled'});
    set_axes_top_text(ff.hf,ff.h_axes(4),'Comp2 (No-Brake)',[-0.02 0 0.2 0]);
    
   
    
     save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs_3_factors.pdf',ntrials),600);
    %%
    break;
end 


%% all exc and inh (all_gV variable) percentage of cells, response fidelity, and zMI
while 1
    clc
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 2],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',...
        [10 -410]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 1.5 1.25]);
    MY = 15; ysp = 0.75; sigLinesStartYFactor = 1.5; mY = 0; % responsive cells
%     MY = 70; ysp = 5; sigLinesStartYFactor = 1.5; mY = 0; % responsive fidelity
%     MY = 1.5; ysp = 0.15; sigLinesStartYFactor = 0.1; mY = -0.25; % responsive fidelity
    stp = 0.43; widths = [0.3 0.4 0.4 0.4 0.4 0.4]+0.1; gap = 0.09;
    set(ff.hf,'Units','inches');
    for aii = 1:length(ff.h_axes)
        sel_ax = ff.h_axes(1,aii);
        set(sel_ax,'Units','inches');
        axesPosd = get(sel_ax,'Position');
        if aii == 1
            rt = stp; 
            ylabel(sel_ax,{'Estimated','Marginal Means'});
        else
            rt = rt + widths(aii-1) + gap; set(sel_ax,'yticklabels',[]);
        end
        set(sel_ax,'Position',[rt axesPosd(2) widths(aii) axesPosd(4)]);
        ylim(sel_ax,[mY MY]);
    end

    %++++++++++++++++++++++++
    axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_cond_et,{'Cond','bonferroni'},[1.5 1 1]); 
	xdata = make_xdata([2],[1 1.5]);    tcolors = repmat(mData.colors([1 2]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'B','NB'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    %+++++++++++++++++++++++
     axes(ff.h_axes(1,2));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_cond_et,{'ET','bonferroni'},[1.5 1 1]); 
	xdata = make_xdata([4],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_axes_top_text(ff.hf,ff.h_axes(1:2),'All Responsive Cells');
 
     save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs_3_factors.pdf',ntrials),600);
    %%
    break;
end 


%% all exc and inh (all_gV variable) complementation conjunction
while 1
    clc
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 2],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',...
        [10 -410]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 1.6 1.25]);
    MY = 30; ysp = 2; sigLinesStartYFactor = 1.5; mY = 0; % responsive cells
%     MY = 70; ysp = 5; sigLinesStartYFactor = 1.5; mY = 0; % responsive fidelity
%     MY = 1.5; ysp = 0.15; sigLinesStartYFactor = 0.1; mY = -0.25; % responsive fidelity
    stp = 0.43; widths = [0.4 0.4 0.4 0.4 0.4 0.4]+0.1; gap = 0.09;
    set(ff.hf,'Units','inches');
    for aii = 1:length(ff.h_axes)
        sel_ax = ff.h_axes(1,aii);
        set(sel_ax,'Units','inches');
        axesPosd = get(sel_ax,'Position');
        if aii == 1
            rt = stp; 
            ylabel(sel_ax,{'Estimated','Marginal Means'});
        else
            rt = rt + widths(aii-1) + gap; set(sel_ax,'yticklabels',[]);
        end
        set(sel_ax,'Position',[rt axesPosd(2) widths(aii) axesPosd(4)]);
        ylim(sel_ax,[mY MY]);
    end
% 
%     %++++++++++++++++++++++++
    axes(ff.h_axes(1,2));
%     delete(gca);
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_et_pt,{'ET','bonferroni'},[1.5 1 1]); 
	xdata = make_xdata([3],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    %+++++++++++++++++++++++
     axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_et_pt,{'PT','bonferroni'},[1.5 1 1]); 
	xdata = make_xdata([3],[1 1.5]);    tcolors = repmat(mData.colors([3 5 10]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'Comp1','Conj','Comp2'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
%     ht = title('All Responsive Cells'); set(ht,'FontSize',6,'FontWeight','Normal')
    set_axes_top_text(ff.hf,ff.h_axes(1:2),'All Responsive Cells',[0 0 0 0]);
 
     save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs_3_factors.pdf',ntrials),600);
    %%
    break;
end 

