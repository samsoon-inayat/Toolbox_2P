function trial_to_trial_and_cond_corr

%% find spatial across condition corr
while 1
    %%
    si = [Ar_t_D ArL_t_D Ars_t_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    gauss = props1.vals; %n_gauss = props1.good_FR_and_notGauss;
    respE_OR = cell_list_op(gauss,[],'or'); outT = find_population_vector_corr_remap(Rs,mR,respE_OR);
%     respC_OR = cell_list_op(n_gauss,[],'or'); outU = find_population_vector_corr_remap(Rs,mR,respC_OR);
%     respR_OR = cell_list_op(props1.good_FR,[],'or'); outR = find_population_vector_corr_remap(Rs,mR,respR_OR);
    disp('Done');
    n = 0;
    %%
    break;
end
    
%% spatial across conditions corr
while 1
    meancorr_trials = [];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outT.adj_SP_corr_diag,'nanmean')];
    [within,dvn,xlabels] = make_within_table({'Cond'},[2]);
    dataT = make_between_table({meancorr_trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1]);
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
    xticks = xdata; xticklabels = {'3-4','4-5'};
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
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outT.adj_PV_corr_diag,'nanmean')];
    [within,dvn,xlabels] = make_within_table({'Cond'},[2]);
    dataT = make_between_table({meancorr_trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1]);
    hf = get_figure(5,[8 7 1.25 1]);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3-4','4-5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.5 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Pop. Vec. ','Correlation(%)'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('spatial_across_conditions_corr_PV.pdf'),600);
    %%
    break;
end

%% spatial across conditions corr RR
while 1
    meancorr_trials = [];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outT.adj_RR_SP,'nanmean')];
    [within,dvn,xlabels] = make_within_table({'Cond'},[2]);
    dataT = make_between_table({meancorr_trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
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
    xticks = xdata; xticklabels = {'3-4','4-5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.5 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Rate Remapping','Score'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('spatial_across_conditions_corr_RR.pdf'),600);
    %%
    break;
end

