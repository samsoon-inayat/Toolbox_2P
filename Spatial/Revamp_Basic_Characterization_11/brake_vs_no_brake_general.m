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
    figdim = 4;
    hf = get_figure(5,[5 2 figdim figdim]);
    
    good_FR = [all_exc(:,1),all_inh(:,1),all_exc(:,2),all_inh(:,2),all_exc(:,3),all_inh(:,3),all_exc(:,4),all_inh(:,4)]; % for comparison of air onset with air offset across brake and no-brake
    txl = {'Exc','Inh'}; txl = repmat(txl,1,4);
    
    good_FR = [all_exc(:,1),all_inh(:,1),all_exc(:,2),all_inh(:,2),all_exc(:,5),all_inh(:,5),all_exc(:,3),all_inh(:,3),all_exc(:,4),all_inh(:,4),all_exc(:,6),all_inh(:,6)]; % for comparison of air onset offset and arb
    txl = {'Exc','Inh'}; txl = repmat(txl,1,6);
    
    good_FR = all_exc_inh;
    txl = xtl;
    
    good_FR = all_gV;
    txl = event_type;
    
   
    [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(good_FR,0.5,0.05);
    mOI = mCI; semOI = semCI;
    
    mUni = nanmean(uni,3); mUni1 = tril(mUni,-1) + tril(mUni,-1)'; mUni2 = triu(mUni,1) + triu(mUni,1)'; mmUni = min(mUni(:)); MmUni = max(mUni(:));
%     figure(1000);clf;subplot 131;imagesc(mUni,[mmUni MmUni]);set(gca,'YDir','normal');subplot 132;imagesc(mUni1,[mmUni MmUni]);set(gca,'YDir','normal');subplot 133;imagesc(mUni2,[mmUni MmUni]);set(gca,'YDir','normal');
    
%     mOI = mUni2;
    
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
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    for rr = 1:size(mCI,1)
        for cc = 1:size(mCI,1)
            if rr == cc
                text(rr-0.25,cc,sprintf('%.0f',mCI(rr,cc)),'FontSize',5);
            end
        end
    end
    plot([(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],[0 150.5],'w','linewidth',0.1); 
    plot([0 150.5],[(10.5+((ii-1)*10)) (10.5+((ii-1)*10))],'w','linewidth',0.1); 
    set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.0 0 0.0 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.1 0.05]);
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
    set(hf,'Position',[7 3 3.5 3]);
%     set(hf,'Position',[7 3 2 1]);
    set(H,'linewidth',0.5);
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel({'Eucledian','Distance'});changePosition(hx,[0 -0.3 0]);
    changePosition(gca,[0 0.0 0.07 -0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);
    %%
    break;
end


%% for bar graphs type_by_cond
while 1
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Type_by_Cond','hsd'},[1.5 1 1]);
	xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-Exc','NB-Exc','B-Inh','NB-Inh'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.5 -0.04]); 
%     put_axes_labels(gca,{[],[0 0 0]},{{'Trials (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end
