function trial_to_trial_and_cond_corr

%% find spatial across condition corr
while 1
    %%
    si = [Ar_t_D ArL_t_D Ars_t_D];
    si = [Ar_i_T ArL_i_T Ars_i_T];
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
    
%% spatial and PV across conditions corr
while 1
    meancorr_trials = [];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outT.adj_SP_corr_diag,'nanmean')];
    meancorr = [];
    for an = 1:5
        ttt = outT.all_PV_corr_diag{an};
%         ttt = outT.all_SP_corr_diag{an};
%         ttt = outT.all_RR_SP{an};
        tttt = exec_fun_on_cell_mat(ttt,'nanmean');
        meancorr(an,:) = [tttt(1,2),tttt(2,3),tttt(1,3)];
    end
    
    [within,dvn,xlabels] = make_within_table({'Cond'},[3]);
    dataT = make_between_table({meancorr},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY ;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3-4','4-5','3-5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.08 0.01 -0.4 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Correlation',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('spatial_across_conditions_corr.pdf'),600);
    %%
    break;
end


%% RR score across conditions corr
while 1
%     meancorr_trials = [];
%     meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outT.adj_SP_corr_diag,'nanmean')];
    meancorr = [];
    for an = 1:5
        ttt = outT.all_RR_SP{an};
        tttt = exec_fun_on_cell_mat(ttt,'nanmean');
        meancorr(an,:) = [tttt(1,2),tttt(2,3),tttt(1,3)];
    end
    
    [within,dvn,xlabels] = make_within_table({'Cond'},[3]);
    dataT = make_between_table({meancorr},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
    xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    maxY = maxY ;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0.7 maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3-4','4-5','3-5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0.7 0.8 0.9 1]); xtickangle(45)
    changePosition(gca,[0.08 0.01 -0.4 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Score',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('spatial_across_conditions_corr.pdf'),600);
    %%
    break;
end

