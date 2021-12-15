function trials_resp_based_analysis

%% spatial
si_nb_T = [Ar_t_D ArL_t_D Ars_t_D]; si_nb_IT = [Ar_i_T ArL_i_T Ars_i_T];
si_b_A = [Ab_T,Abs_T];
si = si_nb_T;
si = si_b_A;
si = [Lb_T ArL_L_T Lbs_T Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
si = [Ar_i_T ArL_i_T Ars_i_T];
si = [Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
ntrials = {[0,20],[20,40],[40,60],[60,80],[80,100]};
for ii = 1:length(ntrials)
    props{ii} = get_props_Rs(o.Rs(:,si),ntrials{ii});
end
disp('done')
    %% For percentage of cells over belt
    perc = [];
    for ii = 1:length(ntrials)
        resp = props{ii}.good_FR; tperc = 100*exec_fun_on_cell_mat(resp,'sum')./exec_fun_on_cell_mat(resp,'length');
        perc = [perc tperc];
    end
    perc1 = [];
    for ii = 1:length(si)
        perc1 = [perc1 perc(:,[ii:length(si):(length(ntrials)*length(si))])];
    end
    [within,dvn,xlabels] = make_within_table({'Cond','CT'},[length(si),length(ntrials)]);
    dataT = make_between_table({perc1},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_CT','hsd'},[1 1 1]);
    xdata = make_xdata([repmat(length(ntrials),1,length(si))],[1 2]); 
%     h(h==1) = 0;
    hf = get_figure(5,[8 7 5.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.dcolors(1:length(ntrials)),length(si),1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
%     make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'R0','R20','R40','R60','R80'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
%     changePosition(gca,[0.05 0.0 -0.5 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('percent_cells_trials_resp.pdf'),600);
    
%% For zMIs of cells over belt
    zMI10 = exec_fun_on_cell_mat(props_nb_T10.zMI,'nanmean',props_nb_T10.good_FR);
    zMI40 = exec_fun_on_cell_mat(props_nb_T40.zMI,'nanmean',props_nb_T40.good_FR);
    zMI70 = exec_fun_on_cell_mat(props_nb_T70.zMI,'nanmean',props_nb_T70.good_FR);

    [within,dvn,xlabels] = make_within_table({'CT','Cond'},[3,3]);
    dataT = make_between_table({zMI10,zMI40,zMI70},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([4],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(1:4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.85,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
%     make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 6]); format_axes(gca);
    xticks = xdata; xticklabels = {'R10','R40','R70'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.5 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_cells_trials_resp.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Bin','hsd'},[1 1 1]);
    xdata = make_xdata([3],[1 2]); 
    hf = get_figure(6,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(3:5);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 6]); format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.4 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_cells_on_belt_pooled_Bin.pdf'),600);


