%% start with this
while 1
    siS = [Ar_t_D ArL_t_D Ars_t_D Ar_t_D ArL_t_D Ars_t_D]; si = siS; Rs = o.Rs(:,si); mR = o.mR(:,si);
    props1S = get_props_Rs(Rs,50);
    siD = [Ar_t_D ArL_t_D Ars_t_D]; si = siD; RsD = o.Rs(:,si); mRD = o.mR(:,si);
    siT = [Ar_t_T ArL_t_T Ars_t_T]; si = siT; RsT = o.Rs(:,si); mRT = o.mR(:,si);
    propsD = get_props_Rs(RsD,50); propsT = get_props_Rs(RsT,50);
    
    dzMI = prop_op(propsD,propsT,0);
%     gFR_D_g_T = cell_list_op(props1S.good_FR,dzMI.resp_D_g_T,'and');
%     gFR_T_g_D = cell_list_op(props1S.good_FR,dzMI.resp_T_g_D,'and');
    mean_diff = exec_fun_on_cell_mat(dzMI.diff_D_T,'nanmean');
    mDiff = mean(mean_diff(:)); sDiff = std(mean_diff(:));
    dzMI = prop_op(propsD,propsT,mDiff+sDiff);
    gFR_D_g_T = dzMI.resp_D_g_T;
    gFR_T_g_D = dzMI.resp_T_g_D;
    gFR_D_g_T = cell_list_op(dzMI.resp_D_g_T,propsD.good_FR,'and');
    gFR_T_g_D = cell_list_op(dzMI.resp_T_g_D,propsD.good_FR,'and');
    gFR_Comp = cell_list_op(dzMI.resp_complex,propsD.good_FR,'and');
%     [dzMI.resp_D_g_T_perc;dzMI.resp_T_g_D_perc]
    gauss = propsD.good_FR_and_Gauss_loose; n_gauss = propsD.good_FR_and_notGauss_loose;
    break;
end

%% compare responsivity in distance and time rasters
while 1
    resp = [propsD.good_FR propsT.good_FR];
    perc_resp = 100*exec_fun_on_cell_mat(resp,'sum')./exec_fun_on_cell_mat(resp,'length');
    [within,dvn,xlabels] = make_within_table({'DiTi','Cond'},[2,3]);
    dataT = make_between_table({perc_resp},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DiTi_TI_Cond','hsd'},[1 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3 3],[1 2]); 
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.colors(1:6),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.25,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(7:end));
    ylims = ylim;
    format_axes(gca);
    maxY1 = maxY;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Ar-t-D','ArL-t-D','Ar*-t-D','Ar-i-D','ArL-i-D','Ar*-i-D','Ar-t-T','ArL-t-T','Ar*-t-T','Ar-i-T','ArL-i-T','Ar*-i-T'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.07 0.0 -0. -0.03]); put_axes_labels(gca,{[],[0 0 0]},{{'Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('responsive_cells_zMID_zMIT.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DiTi_by_Cond','hsd'},[1 1 1]);
    xdata = make_xdata([3 3],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.dcolors(1:3),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.25,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(3:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Trials','I-Trials'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 50]); xtickangle(45);
    changePosition(gca,[0.15 0.0 -0.25 -0.03]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('responsive_cells_zMID_zMIT_pooled_DiTi_TI.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DiTi','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.dcolors(3),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(2:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY1]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dist','Time'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.5 -0.08]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('responsive_cells_zMID_zMIT_pooled_DiTi.pdf'),600);
    
    %%
    break;
end

%% compare the difference between zMID and zMIT across trials and intertrials of No-Brake Conditions
while 1
    mean_diff = exec_fun_on_cell_mat(dzMI.diff_D_T,'nanmean');
    [within,dvn,xlabels] = make_within_table({'TI','Cond'},[2,3]);
    dataT = make_between_table({mean_diff},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_Cond','hsd'},[1 1 1]);
    xdata = make_xdata([3 3],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(1:6);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.25,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Ar','ArL','Ar*'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.09 0.0 -0.5 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{'zMI Difference'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('mean_diff_zMID_zMIT.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'TI','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.dcolors(1:2),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
%     make_bars_hollow(hbs(2:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Trials','I-Trials'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.09 0.0 -0.5 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{'zMI Difference'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('mean_diff_zMID_zMIT_pooled.pdf'),600);
    %%
    break;
end

%% population vector and correlation spatial temporal
while 1
    cell_type = 'time_enc';
%     cell_type = 'time_enc';
    titles = {'Ar-t-D','ArL-t-D','Ar*-t-D','Ar-t-D','ArL-t-D','Ar*-t-D'};
    an = 4;
%     resp = cell_list_op(props1S.good_FR,dzMI.resp_D_g_T,'and');
%     resp = cell_list_op(propsT.good_FR,dzMI.resp_T_g_D,'and');
    resp = [dzMI.resp_D_g_T dzMI.resp_T_g_D];
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 6],...
        'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.055 0.1],'widthHeightAdjustment',...
        [25 -60]);    set(gcf,'color','w');    set(gcf,'Position',[10 5 5 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('population_vector_corr_dist_time_%s.pdf',cell_type),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 6],...
        'spaceRowsCols',[0 -0.024],'rightUpShifts',[0.055 0.25],'widthHeightAdjustment',...
        [20 -310]);    set(gcf,'color','w');    set(gcf,'Position',[10 8 5 0.85]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('average_population_vector_corr_dist_time_%s.pdf',cell_type),600);
    %%
    break;
end

%% Comparison of zMIs across the conditions for distance and time zMIs and trials and intertrials
while 1
    %%
    all_zMID = exec_fun_on_cell_mat(propsD.zMI,'nanmean');
    all_zMIT = exec_fun_on_cell_mat(propsT.zMI,'nanmean');
    [within,dvn,xlabels] = make_within_table({'DiTi','Cond'},[2,3]);
    dataT = make_between_table({all_zMID,all_zMIT},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DiTi_by_Cond','hsd'},[1 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3 3],[1 2]); 
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.colors(1:6),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.25,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(7:end));
    ylims = ylim;
    format_axes(gca);
    maxY1 = maxY;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Ar-t-D','ArL-t-D','Ar*-t-D','Ar-i-D','ArL-i-D','Ar*-i-D','Ar-t-T','ArL-t-T','Ar*-t-T','Ar-i-T','ArL-i-T','Ar*-i-T'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.03 0.0 0.03 -0.03]); put_axes_labels(gca,{[],[0 0 0]},{{'zMI'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('zMIs_D_T.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DiTi_by_Cond','hsd'},[1 1 1]);
    xdata = make_xdata([3 3],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.dcolors(1:3),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     maxY = maxY;
    make_bars_hollow(hbs(3:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Trials','I-Trials'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    changePosition(gca,[0.15 0.0 -0.25 -0.03]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('zMIs_DiTI_TI_pooled.pdf'),600);
    %%
    break;
end


%% Comparison of the difference of distance and time zMIs of responsive cells during trials versus intertrials
while 1
    %%
    resp = [propsD.good_FR(:,1:3) propsT.good_FR(:,4:6)];
    respDo = propsD.good_FR(:,1:3); respTo = propsT.good_FR(:,4:6);
    respD = cell_list_op(respDo,[],'or'); respT = cell_list_op(respTo,[],'or');
    respDf = cell_list_op(respD(:,1),respT(:,1),'sep');
    respTf = cell_list_op(respT(:,1),respD(:,1),'sep');
%     dzMI = prop_op(propsD,propsT,0);
    dist_diff_zMI1 = [];
    for rr = 1:size(dzMI.diff_D_T,1)
        t_rD = respDf{rr,1}; t_rT = respTf{rr,1};
        dist_diff_zMI1{rr,1} = [dzMI.diff_D_T{rr,1};dzMI.diff_D_T{rr,3};dzMI.diff_D_T{rr,5}];
        dist_diff_zMI1{rr,2} = [dzMI.diff_D_T{rr,2};dzMI.diff_D_T{rr,4};dzMI.diff_D_T{rr,6}];
    end
%     dist_diff_zMI1 = dzMI.diff_D_T;
    [distDo,allVals] = getAveragesAndAllValues(dist_diff_zMI1);
    minBin = min(allVals);
    maxBin = max(allVals);
    incr = 1;
    tcolors = mData.colors;
    hf = get_figure(8,[5 7 2.25 1.5]);hold on;
    [ha,hb,~,bins] = plotAverageDistributions(dist_diff_zMI1,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'pdf_or_cdf','pdf');
    %%
    resp_sel = [repmat(respDf,1,3) repmat(respTf,1,3)];
    diff_zMI = exec_fun_on_cell_mat(dzMI.diff_D_T,'nanmean',resp_sel);
    [within,dvn,xlabels] = make_within_table({'TI','Cond'},[2,3]);
    dataT = make_between_table({diff_zMI},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_Cond','hsd'},[1 1 1]);
    ptab = 0;
    if ptab h(h==1) = 0; end
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    
    if ptab
    [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.5);
    else
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'yspacing',0.2);
    end
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    txl = rasterNamesTxt([siD siT]); 
    xticks = xdata; xticklabels = {'Trials','Inter-Trials'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 (maxY)]); xtickangle(45);
    if 0
    changePosition(gca,[0.08 0.01 0.0 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. Correlation',[-1 3 0]});
    ha = gca; ptable = extras.pvalsTable;
    display_p_table(ha,hbs,[0 -0.07 0 0.9],ptable);
    else
    changePosition(gca,[-0.01 -0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. zMI Difference',[0 0 0]});
    end
    %%
    break;
end


%% Comparison of the difference of distance and time zMIs of responsive cells 
while 1
    %%
    dzMI = prop_op(propsD.zMI,propsT.zMI,0.5,props1S.good_FR);
    dist_diff_zMI1 = [];
    for rr = 1:size(dzMI.diff_D_T,1)
        dist_diff_zMI1{rr,1} = [dzMI.diff_D_T{rr,1};dzMI.diff_D_T{rr,3};dzMI.diff_D_T{rr,5}];
        dist_diff_zMI1{rr,2} = [dzMI.diff_D_T{rr,2};dzMI.diff_D_T{rr,4};dzMI.diff_D_T{rr,6}];
    end
%     dist_diff_zMI1 = dzMI.diff_D_T;
    [distDo,allVals] = getAveragesAndAllValues(dist_diff_zMI1);
    minBin = min(allVals);
    maxBin = max(allVals);
    incr = 1;
    tcolors = mData.colors;
    hf = get_figure(8,[5 7 2.25 1.5]);hold on;
    [ha,hb,~,bins] = plotAverageDistributions(dist_diff_zMI1,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'pdf_or_cdf','pdf');
    %%
    diff_zMI = exec_fun_on_cell_mat(dzMI.diff_D_T,'nanmean');
    [within,dvn,xlabels] = make_within_table({'T','Cond'},[2,3]);
    dataT = make_between_table({diff_zMI},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_T','hsd'},[1 1 1]);
    ptab = 0;
    if ptab h(h==1) = 0; end
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    
    if ptab
    [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.5);
    else
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'yspacing',0.2);
    end
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    txl = rasterNamesTxt([siD siT]); 
    xticks = xdata; xticklabels = {'Trials','Inter-Trials'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 (maxY)]); xtickangle(45);
    if 0
    changePosition(gca,[0.08 0.01 0.0 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. Correlation',[-1 3 0]});
    ha = gca; ptable = extras.pvalsTable;
    display_p_table(ha,hbs,[0 -0.07 0 0.9],ptable);
    else
    changePosition(gca,[-0.01 -0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. zMI Difference',[0 0 0]});
    end
    %%
    break;
end

%% Comparison of the percentages with zMID > zMIT and zMIT > zMID (differential spatial and temporal encoding during trials and inter-trials
while 1
    %%
    perc_D = 100*exec_fun_on_cell_mat(gFR_D_g_T,'sum')./exec_fun_on_cell_mat(gFR_D_g_T,'length');
    perc_T = 100*exec_fun_on_cell_mat(gFR_T_g_D,'sum')./exec_fun_on_cell_mat(gFR_T_g_D,'length');
    
    [within,dvn,xlabels] = make_within_table({'DiTi','Cond'},[2,3]);
    dataT = make_between_table({perc_D,perc_T},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DiTi','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.dcolors(1:2),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
%     make_bars_hollow(hbs(2:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dist','Time'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.09 0.0 -0.5 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('mean_percCells_pooled_dzMI_based.pdf'),600);
   %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DiTi_by_Cond','hsd'},[1 1 1]);
    xdata = make_xdata([3 3],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.dcolors(1:3),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.25,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(3:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Trials','I-Trials'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    changePosition(gca,[0.15 0.0 -0.25 -0.03]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('diff_ST_all_pooled.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DiTi_TI_Cond','hsd'},[1 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3 3],[1 2]);
    hf = get_figure(5,[8 7 3 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;

    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'yspacing',10);

    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    txl = rasterNamesTxt([siD siT]); 
    xticks = xdata; xticklabels = xlabels%{'DistE','TimeE','DistE','TimeE'};
    xticks = xdata; xticklabels = {'Ar-t','ArL-t','Ar*-t','Ar-i','ArL-i','Ar*-i','Ar-t','ArL-t','Ar*-t','Ar-i','ArL-i','Ar*-i'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
%     set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 (100)]); xtickangle(45);
    changePosition(gca,[-0.01 -0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
%     save_pdf(hf,mData.pdf_folder,sprintf('diff_ST_all.pdf'),600);
    %%
    break;
end

%% Comparison of the peak locations for zMID > zMIT and zMIT > zMID from distance based rasters
while 1
    %%
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            DgT = dzMI.resp_D_g_T{rr,cc};
            TgD = dzMI.resp_T_g_D{rr,cc};
            tpeak_loc = propsD.peak_locations{rr,cc};
            peak_locD(rr,cc) = nanmean(tpeak_loc(DgT));
            peak_locT(rr,cc) = nanmean(tpeak_loc(TgD));
        end
    end
    
    %%
    [within,dvn,xlabels] = make_within_table({'DiTi','TI','Cond'},[2,2,3]);
    dataT = make_between_table({peak_locD,peak_locT},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    
   %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DiTi_by_TI','hsd'},[1 1 1]);
    xdata = make_xdata([2 2],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.dcolors(1:2),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     maxY = maxY;
    make_bars_hollow(hbs(3:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Trials','I-Trials'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    changePosition(gca,[0.15 0.0 -0.25 -0.03]); put_axes_labels(gca,{[],[0 0 0]},{{'Peak Locations (cm)'},[0 -20 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Peak_Locations_DiTI_TI_pooled_dist.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DiTi_TI_Cond','hsd'},[1 1 1]);
    xdata = make_xdata([3 3 3 3],[1 2]);
    hf = get_figure(5,[8 7 6.99 1.5]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;

    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'yspacing',10);

    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    txl = rasterNamesTxt([siD siT]); 
    xticks = xdata; xticklabels = xlabels%{'DistE','TimeE','DistE','TimeE'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 (100)]); xtickangle(45);
    changePosition(gca,[-0.01 -0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Peak Locations (cm)',[0 0 0]});
%     save_pdf(hf,mData.pdf_folder,sprintf('diff_ST_all.pdf'),600);
    %%
    break;
end

%% time Comparison of the peak locations for zMID > zMIT and zMIT > zMID from time based rasters
while 1
    %%
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            DgT = dzMI.resp_D_g_T{rr,cc};
            TgD = dzMI.resp_T_g_D{rr,cc};
            tpeak_loc = propsT.peak_locations{rr,cc};
            peak_locD(rr,cc) = nanmean(tpeak_loc(DgT));
            peak_locT(rr,cc) = nanmean(tpeak_loc(TgD));
        end
    end
    
    %%
    [within,dvn,xlabels] = make_within_table({'DiTi','TI','Cond'},[2,2,3]);
    dataT = make_between_table({peak_locD,peak_locT},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    
   %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DiTi_by_TI','hsd'},[1 1 1]);
    xdata = make_xdata([2 2],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.dcolors(1:2),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     maxY = maxY;
    make_bars_hollow(hbs(3:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Trials','I-Trials'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    changePosition(gca,[0.15 0.0 -0.25 -0.03]); put_axes_labels(gca,{[],[0 0 0]},{{'Peak Locations (sec)'},[0 -3 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Peak_Locations_DiTI_TI_pooled_time.pdf'),600);
    
    %%
    break;
end

%% Overlap Indices ImageSC zMID>zMIT
while 1
    ntrials = 50;
    si = siS;
    resp = [dzMI.resp_D_g_T props1S.good_FR];% resp_speed];
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = [{'D-Ar-t-D'}    {'D-ArL-t-D'}    {'D-Ar*-t-D'}    {'D-Ar-i-D'}    {'D-ArL-i-D'}    {'D-Ar*-i-D'} rasterNamesTxt(si)]; 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.5 1.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    text(0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean_DT.pdf',ntrials),600);
    %%
    break;
end

%% agglomerative hierarchical clustering zMID>zMIT
while 1
    mOI1 = mOI;
    mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1);
    tree = linkage(Di);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','right','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 1.25 2]);
    set(H,'linewidth',1);
    set(gca,'yticklabels',txl(TC));ytickangle(30);
    format_axes(gca);
    hx = xlabel('Eucledian Distance');%changePosition(hx,[-0.051 0 0]);
    changePosition(gca,[0 0.0 0.05 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster_DT.pdf'),600);
    %%
    break;
end


%% Overlap Indices ImageSC zMIT>zMID
while 1
    ntrials = 50;
    si = siS;
    resp = [dzMI.resp_T_g_D props1S.good_FR];% resp_speed];
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = [{'T-Ar-t-T'}    {'T-ArL-t-T'}    {'T-Ar*-t-T'}    {'T-Ar-i-T'}    {'T-ArL-i-T'}    {'T-Ar*-i-T'} rasterNamesTxt(si)]; 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.5 1.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    text(0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean_DT_T.pdf',ntrials),600);
    %%
    break;
end

%% agglomerative hierarchical clustering zMIT>zMID
while 1
    mOI1 = mOI;
    mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1);
    tree = linkage(Di);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','right','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 1.25 2]);
    set(H,'linewidth',1);
    set(gca,'yticklabels',txl(TC));ytickangle(30);
    format_axes(gca);
    hx = xlabel('Eucledian Distance');%changePosition(hx,[-0.051 0 0]);
    changePosition(gca,[0 0.0 0.05 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster_DT_T.pdf'),600);
    %%
    break;
end


%% compare responsivity in gauss n_gauss D>T and T>D
while 1
    resp = [gauss n_gauss gFR_D_g_T gFR_T_g_D gFR_Comp];
    perc_resp = 100*exec_fun_on_cell_mat(resp,'sum')./exec_fun_on_cell_mat(resp,'length');
    [within,dvn,xlabels] = make_within_table({'DiTi','Cond'},[5,3]);
    dataT = make_between_table({perc_resp},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    
    
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DiTi','hsd'},[1 1 1]);
    xdata = make_xdata([2,3],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.dcolors(1:5),4,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',7,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
%     make_bars_hollow(hbs(2:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'gT','gU','Dist','Time','Comp'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    changePosition(gca,[0.06 0.05 -0.1 -0.1]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('percent_DiTi_trials.pdf'),600);
    
    %%
    break;
end

%% compare zMI in gauss n_gauss D>T and T>D
while 1
    all_zMI_C = exec_fun_on_cell_mat(propsD.zMI,'nanmean',gFR_Comp);
    all_zMI_D = exec_fun_on_cell_mat(propsD.zMI,'nanmean',gFR_D_g_T);
    all_zMI_T = exec_fun_on_cell_mat(propsD.zMI,'nanmean',gFR_T_g_D);
    all_zMI_G = exec_fun_on_cell_mat(propsD.zMI,'nanmean',gauss);
    all_zMI_nG = exec_fun_on_cell_mat(propsD.zMI,'nanmean',n_gauss);
    [within,dvn,xlabels] = make_within_table({'DiTi','Cond'},[5,3]);
    dataT = make_between_table({all_zMI_G,all_zMI_nG,all_zMI_D,all_zMI_T,all_zMI_C},dvn);
    ra = RMA(dataT,within);
    ra.ranova
   
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DiTi','hsd'},[1 1 1]);
    xdata = make_xdata([2 3],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.dcolors(1:5),4,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
%     make_bars_hollow(hbs(2:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'gT','gU','Dist','Time','Comp'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.06 0.05 -0.1 -0.1]); put_axes_labels(gca,{[],[0 0 0]},{{'zMI'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('zMIs_DiTi_trials.pdf'),600);

%%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DiTi_by_Cond','hsd'},[1 1 1]);
    xdata = make_xdata([3 3 3 3 3],[1 2]); 
    hf = get_figure(5,[8 7 3.25 3]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.dcolors(1:5),4,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
%     make_bars_hollow(hbs(2:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Ar','ArL','Ar*'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.06 0.05 -0.1 -0.1]); put_axes_labels(gca,{[],[0 0 0]},{{'zMI'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('zMIs_DiTi_by_Cond_trials.pdf'),600);
    
    %%
    break;
end


%% Overlap Indices ImageSC gauss n gauss d>t t>d
while 1
    ntrials = 50;
    si = siS;
    resp = [gFR_D_g_T gFR_T_g_D gFR_Comp gauss n_gauss];
%     resp = [dzMI.resp_T_g_D props1S.good_FR];% resp_speed];
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = [{'Dist-Ar'}    {'Dist-ArL'}    {'Dist-Ar*'}    {'Time-Ar'}    {'Time-ArL'}    {'Time-Ar*'} {'Comp-Ar'}    {'Comp-ArL'}    {'Comp-Ar*'} {'gT-Ar'}    {'gT-ArL'}    {'gT-Ar*'}    {'gU-Ar'}    {'gU-ArL'}    {'gU-Ar*'} ]; 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 2.5 2.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    text(0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean_DT_T.pdf',ntrials),600);
    %%
    break;
end

%% agglomerative hierarchical clustering 
while 1
    mOI1 = mOI;
    mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1);
    tree = linkage(Di);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 5 2.25 1.25]);
    set(H,'linewidth',1);
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel('Eucledian Distance');changePosition(hx,[0 -0.05 0]);
    changePosition(gca,[0.03 0.0 0.05 0.01]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster_DT_T.pdf'),600);
    %%
    break;
end

