%% for VENN Diagrams with two circles
while 1
    %%
    % for different VENN Diagrams change good_FRV = all_exc(:,[1 3]) or
    % all_exc(:,[2 4]) or all_exc(:,[5 6]) same thing for all_inh
    good_FRV = [cell_list_op(all_gV(:,[1:2 5]),[],'or',1) cell_list_op(all_gV(:,[3:4 6]),[],'or',1)]; 
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
    AVenn = [p_gFR_C(1,1) p_gFR_C(2,2)]; IVenn = [p_gFR_C(1,2)];
    [HVenn] = venn(AVenn,IVenn,'ErrMinMode','None');
%     AVenn = [p_gFR_C(1,1) p_gFR_C(2,2) p_gFR_C(3,3)]; IVenn = [p_gFR_C(1,2) p_gFR_C(1,3) p_gFR_C(2,3) p_gFR_C_all];
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

%% for VENN Diagrams with three circles
while 1
    %%
    % for different VENN Diagrams change good_FRV = all_exc(:,[1 3]) or
    % all_exc(:,[2 4]) or all_exc(:,[5 6]) same thing for all_inh
    good_FRV = [cell_list_op(all_gV(:,[1:2 5]),[],'or',1) cell_list_op(all_gV(:,[3:4 6]),[],'or',1)]; 
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
    hf = get_figure(5,[5 2 figdim figdim]);
    
    good_FR = [all_exc(:,1),all_inh(:,1),all_exc(:,2),all_inh(:,2),all_exc(:,3),all_inh(:,3),all_exc(:,4),all_inh(:,4)]; % for comparison of air onset with air offset across brake and no-brake
    txl = {'AOn-Exc','Inh'}; txl = repmat(txl,1,4);
    
    good_FR = [all_exc(:,1),all_inh(:,1),all_exc(:,2),all_inh(:,2),all_exc(:,5),all_inh(:,5),all_exc(:,3),all_inh(:,3),all_exc(:,4),all_inh(:,4),all_exc(:,6),all_inh(:,6)]; % for comparison of air onset offset and arb
%     txl = {'AOn-Exc','AOn-Inh','AOff-Exc','AOff-Inh','Arb-Exc','Arb-Inh','AOn-Exc','AOn-Inh','AOff-Exc','AOff-Inh','Arb-Exc','Arb-Inh'};
    txl = {'B-AOn-Exc','B-AOn-Inh','B-AOff-Exc','B-AOff-Inh','B-Arb-Exc','B-Arb-Inh','NB-AOn-Exc','NB-AOn-Inh','NB-AOff-Exc','NB-AOff-Inh','NB-Arb-Exc','NB-Arb-Inh'};
%     txl = {'Exc','Inh'}; txl = repmat(txl,1,6);
    
%     good_FR = all_exc_inh;
%     txl = xtl;
%     
%     good_FR = all_gV;
%     txl = event_type;
    
   
    [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(good_FR,0.5,0.05);
    mOI = mCI; semOI = semCI;
    
    mUni = nanmean(uni,3); mUni1 = tril(mUni,-1) + tril(mUni,-1)'; mUni2 = triu(mUni,1) + triu(mUni,1)'; mmUni = min(mUni(:)); MmUni = max(mUni(:));
%     figure(1000);clf;subplot 131;imagesc(mUni,[mmUni MmUni]);set(gca,'YDir','normal');subplot 132;imagesc(mUni1,[mmUni MmUni]);set(gca,'YDir','normal');subplot 133;imagesc(mUni2,[mmUni MmUni]);set(gca,'YDir','normal');
    
    mOI = mUni2;
    
    sz = size(mOI,1);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
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
    plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 150.5],'w','linewidth',0.1); 
    plot([0 150.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
    set(gca,'Ydir','normal');ytickangle(15);
    box on
    changePosition(gca,[0.0 0 0.0 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.08 0.05 0.05 0.05]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',ntrials),600);
    %%
    break;
end

%% Agglomerative hierarchical clustering
while 1
    mOI1 = mOI;
%     mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1,@naneucdist);
%     tree = linkage(mOI1,'average','euclidean');
    tree = linkage(Di,'average');
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 3.5 1.75]);
    set(hf,'Position',[7 3 2.3 1.5]);
    set(H,'linewidth',0.5);
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel({'Eucledian Distance'});changePosition(hx,[0 0 0]);
    changePosition(gca,[0.0 0.0 0.07 0.05]);
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
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
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
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
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
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
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
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
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
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
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
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
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


%% for bar graphs three factors
while 1
   %%
   clc
   hs_flag = logical([0 0 1]);
   ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 7],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.05 0.23],'widthHeightAdjustment',...
        [10 -350]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 4.85 1]);
    MY = 2; ysp = 0.1;
    MY = 2.2; ysp = 0.001; mY = 1.8;
    axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_ET_CT','hsd'},[1.5 1 1]);
    h(h==1) = 0; xdata = make_xdata([2 2 2 2 2 2],[1 1.5]);    tcolors = repmat(mData.colors([7 3]),1,6);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca); xticks = xdata; xticklabels = {'Exc','Inh'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,2));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'B','NB'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,3));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'ET','hsd'},[1.5 1 1]);
	xdata = make_xdata([3],[1 1.5]);    tcolors = repmat(mData.colors([5 6 9]),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,4));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);    tcolors = repmat(mData.colors([3 7]),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);   xticks = xdata; xticklabels = {'Exc','Inh'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,5));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_ET','hsd'},[1.5 1 1]);
    if ~hs_flag(1) h(h==1) = 0; end
	xdata = make_xdata([3 3],[1 1.5]);    tcolors = repmat(mData.colors([5 6 9]),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,6));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_CT','hsd'},[1.5 1 1]);
    if ~hs_flag(2) h(h==1) = 0; end
	xdata = make_xdata([2 2],[1 1.5]);    tcolors = repmat(mData.colors([3 7]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'Exc','Inh'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    axes(ff.h_axes(1,7));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'ET_by_CT','hsd'},[1.5 1 1]);
    if ~hs_flag(3) h(h==1) = 0; end
	xdata = make_xdata([2 2 2],[1 1.5]);    tcolors = repmat(mData.colors([3 7]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'Exc','Inh'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    %+++++++++++++++++++++++++
%     MY = 3;
    set(ff.hf,'Units','inches');
    stp = 0.24; widths = [1.5 0.2 0.3 0.2 0.6 0.5 0.6]; gap = 0.09;
    for aii = 1:7
        sel_ax = ff.h_axes(1,aii);
        set(sel_ax,'Units','inches');
        axesPosd = get(sel_ax,'Position');
        if aii == 1
            rt = stp; ylabel(sel_ax,'Cells (%)');
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



%% for bar graphs CT_by_PT
while 1
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT_by_PT','hsd'},[1.5 1 1]);
    h(h==1) = 0;
	xdata = make_xdata([3 3],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = [mData.colors(1) mData.colors(10) mData.colors(2) mData.colors(1) mData.colors(10) mData.colors(2)];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 9.5]); format_axes(gca);
    xticks = xdata; xticklabels = {'B','Con','NB'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.4 -0.1]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('bar_graph.pdf',ntrials),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'PT','hsd'},[1.5 1 1]);
	xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors([1 10 2]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 9.5]); format_axes(gca);
    xticks = xdata; xticklabels = {'B','Con','NB'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.6 -0.1]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('bar_graph.pdf',ntrials),600);
    %%
    break;
end



%% for bar graphs ET_by_PT
while 1
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'ET_by_PT','hsd'},[1.5 1 1]);
    h(h==1) = 0;
	xdata = make_xdata([3 3 3],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = [mData.colors(1) mData.colors(10) mData.colors(2) mData.colors(1) mData.colors(10) mData.colors(2) mData.colors(1) mData.colors(10) mData.colors(2)];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 30]); format_axes(gca);
    xticks = xdata; xticklabels = {'B','Con','NB'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.2 -0.1]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('bar_graph.pdf',ntrials),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'PT','hsd'},[1.5 1 1]);
	xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors([1 10 2]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',-0.05);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 30]); format_axes(gca);
    xticks = xdata; xticklabels = {'B','Con','NB'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.55 -0.1]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('bar_graph.pdf',ntrials),600);
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'ET','hsd'},[1.5 1 1]);
	xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors([1 10 2]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 30]); format_axes(gca);
    xticks = xdata; xticklabels = {'AOn','AOff','Arb'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.6 -0.1]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('bar_graph.pdf',ntrials),600);
    %%
    break;
end



%% for bar graphs ET_by_PT
while 1
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_ET_CT','hsd'},[1.5 1 1]);
    h(h==1) = 0;
	xdata = make_xdata([2 2 2 2 2 2],[1 1.5]);
    hf = get_figure(5,[8 7 3.5 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = [mData.colors(1) mData.colors(10) mData.colors(2) mData.colors(1) mData.colors(10) mData.colors(2) mData.colors(1) mData.colors(10) mData.colors(2)];
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 30]); format_axes(gca);
    xticks = xdata; xticklabels = {'Exc','Inh'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.01 0.01 -0.1 -0.1]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('bar_graph.pdf',ntrials),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'ET_by_CT','hsd'},[1.5 1 1]);
	xdata = make_xdata([2 2 2],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',-0.05);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 30]); format_axes(gca);
    xticks = xdata; xticklabels = {'Exc','Inh'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.15 -0.1]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('bar_graph.pdf',ntrials),600);
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'ET','hsd'},[1.5 1 1]);
	xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors([1 10 2]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 30]); format_axes(gca);
    xticks = xdata; xticklabels = {'AOn','AOff','Arb'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.6 -0.1]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('bar_graph.pdf',ntrials),600);
    %%
    break;
end
