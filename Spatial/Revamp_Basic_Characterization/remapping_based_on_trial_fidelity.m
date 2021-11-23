function trial_to_trial_corr


%% find spatial across condition corr
while 1
    %%
    si = [Ar_t_D ArL_t_D Ars_t_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    props{1} = get_props_Rs(Rs,[10,40]);
    props{2} = get_props_Rs(Rs,[40,70]);
    props{3} = get_props_Rs(Rs,[70,100]);
    for pri = 1:length(props)
        t_props = props{pri};
        resp_OR = cell_list_op(t_props.good_FR,[],'or'); remap{pri} = find_population_vector_corr_remap(Rs,mR,resp_OR);
    end
    disp('Done');
    n = 0;
    %%
    break;
end
    
%% spatial across conditions corr SP
while 1
    meancorr_trials = [];
    for pri = 1:length(props)
        meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(remap{pri}.adj_SP_corr_diag,'nanmean')];
    end
    [within,dvn,xlabels] = make_within_table({'pop','Cond'},[3,2]);
    dataT = make_between_table({meancorr_trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'pop','hsd'},[1 1 1]);
    xdata = make_xdata([3],[1]);
    hf = get_figure(6,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'TF1','TF2','TF3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.5 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Spatial','Correlation(%)'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('spatial_across_conditions_corr.pdf'),600);
    %%
    break;
end


%% spatial across conditions corr PV
while 1
    meancorr_trials = [];
    for pri = 1:length(props)
        meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(remap{pri}.adj_PV_corr_diag,'nanmean')];
    end
    meancorr_trialsC = meancorr_trials(:,[1 3 5 2 4 6]);
    [within,dvn,xlabels] = make_within_table({'Cond','pop'},[2,3]);
    dataT = make_between_table({meancorr_trialsC},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_pop','hsd'},[1 1 1]);
    xdata = make_xdata([3 3],[1 2]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'TF1','TF2','TF3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.5 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Population','Correlation(%)'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('PV_across_conditions_corr.pdf'),600);
    %%
    break;
end


%% spatial across conditions corr SP _ cell ordering difference
while 1
    meancorr_trials = [];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outR.adj_PV_corr_diag,'nanmean')];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outR.adj_PV_corr_diag_ind,'nanmean')];
    [within,dvn,xlabels] = make_within_table({'C_ord','Cond'},[2,2]);
    dataT = make_between_table({meancorr_trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'C_ord','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1 2]);
    hf = get_figure(5,[8 7 1.25 1]);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Same Order','Peak Order'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.5 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Spatial','Correlation(%)'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('spatial_across_conditions_corr_PV_ind.pdf'),600);
    %%
    break;
end


%% spatial across conditions corr RR
while 1
    meancorr_trials = [];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outT.adj_RR_SP,'nanmean')];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outU.adj_RR_SP,'nanmean')];
    [within,dvn,xlabels] = make_within_table({'CT','Cond'},[2,2]);
    dataT = make_between_table({meancorr_trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.025,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0.8 0.93]); format_axes(gca);
    xticks = xdata; xticklabels = {'gT','gU'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.5 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Rate Remapping','Score'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('spatial_across_conditions_corr_RR.pdf'),600);
    %%
    break;
end

%% Overlap Indices ImageSC spatial
while 1
    ntrials = 50;
    si = [Ab_T Abs_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    propsB = get_props_Rs(Rs,ntrials);
    exc = propsB.good_FR_and_exc;     sup = propsB.good_FR_and_inh;    com = propsB.good_FR_and_untuned;
    respE_OR = cell_list_op(exc,[],'or');     respS_OR = cell_list_op(sup,[],'or');     respC_OR = cell_list_op(com,[],'or');
    
    si_nb = [Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
    props_nb = get_props_Rs(o.Rs(:,si_nb),ntrials);
    resp_nb = props_nb.good_FR;
    
    si = [Ar_t_D ArL_t_D Ars_t_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    propsAt = get_props_Rs(Rs,ntrials);
    gauss = propsAt.good_FR_and_Gauss; n_gauss = propsAt.good_FR_and_notGauss;
    respG_OR = cell_list_op(gauss,[],'or'); respNG_OR = cell_list_op(n_gauss,[],'or');
    
    resp = [respE_OR(:,1) respS_OR(:,1) respC_OR(:,1) gauss n_gauss];
    resp = [gauss n_gauss resp_nb(:,4:6)];
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = [{'Exc','Sup','Com'} {'gT','gT','gT','gU','gU','gU'}]; 
    txl = [{'gT3','gT4','gT5','gU3','gU4','gU5'} rasterNamesTxt(si_nb(4:6))]; 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(6,[8 7 1.5 1.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    text(-0.3,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_tu_spatial.pdf'),600);
    %%
    break;
end

%% spatial agglomerative hierarchical clustering
while 1
    mOI1 = mOI;
    mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1);
    tree = linkage(Di);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','right','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 1.25 1.25]);
    set(H,'linewidth',1);
    set(gca,'yticklabels',txl(TC));ytickangle(30);
    format_axes(gca);
    hx = xlabel('Eucledian Distance');%changePosition(hx,[-0.051 0 0]);
    changePosition(gca,[0 0.0 -0.05 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster_tu_spatial.pdf'),600);
    %%
    break;
end


