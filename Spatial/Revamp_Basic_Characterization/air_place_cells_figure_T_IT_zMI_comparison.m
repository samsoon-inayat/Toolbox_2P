selected_property = 'good_FR';
cmdTxt = sprintf('good_FR = props1.%s;',selected_property);

%% start with this
while 1
    siS = si_seq(setdiff(1:11,[1 11 9 2 10])); si = siS; Rs = o.Rs(:,si); mRs = o.mR(:,si);
    props1S = get_props_Rs(Rs,50);
    siD = si_no_brake_dist; si = siD; RsD = o.Rs(:,si); mRsD = o.mR(:,si);
    siT = si_no_brake_time; si = siT; RsT = o.Rs(:,si); mRsT = o.mR(:,si);
    
    propsD = get_props_Rs(RsD,50); propsT = get_props_Rs(RsT,50);
    
    dzMI = prop_op(propsD.zMI,propsT.zMI,0.3);
    gFR_D_g_T = cell_list_op(props1S.good_FR,dzMI.resp_D_g_T,'and');
    gFR_T_g_D = cell_list_op(props1S.good_FR,dzMI.resp_T_g_D,'and');
%     [dzMI.resp_D_g_T_perc;dzMI.resp_T_g_D_perc]
    break;
end

%% Comparison of zMIs across the conditions for distance and time zMIs and trials and intertrials
while 1
    %%
    siS = si_seq(setdiff(1:11,[1 11 9 2 10])); si = siS; Rs = o.Rs(:,si); mRs = o.mR(:,si);
    props1 = get_props_Rs(Rs,50);
    siD = si_no_brake_dist; si = siD; RsD = o.Rs(:,si); mRsD = o.mR(:,si);
    siT = si_no_brake_time; si = siT; RsT = o.Rs(:,si); mRsT = o.mR(:,si);
    
    propsD = get_props_Rs(RsD,50); propsT = get_props_Rs(RsT,50);
    
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            good_FR = props1.good_FR{rr,cc};
            tzMI = propsD.zMI{rr,cc}; zMID = tzMI(good_FR);
            tzMI = propsT.zMI{rr,cc}; zMIT = tzMI(good_FR);
            all_zMID(rr,cc) = nanmean(zMID); all_zMIT(rr,cc) = nanmean(zMIT);
        end
    end
    
    [within,dvn,xlabels] = make_within_table({'DT','Cond','T'},[2,3,2]);
    dataT = make_between_table({all_zMID,all_zMIT},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DT_Cond_T','hsd'},[1 1 1]);
    ptab = 0;
    if ptab h(h==1) = 0; end
    hf = get_figure(5,[8 7 6.99 1.5]);
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
    xticks = xdata; xticklabels = txl;
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 (maxY)]); xtickangle(45);
    if 0
    changePosition(gca,[0.08 0.01 0.0 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. Correlation',[-1 3 0]});
    ha = gca; ptable = extras.pvalsTable;
    display_p_table(ha,hbs,[0 -0.07 0 0.9],ptable);
    else
    changePosition(gca,[-0.01 -0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. zMI',[0 0 0]});
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
    [within,dvn,xlabels] = make_within_table({'Cond','T'},[3,2]);
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
    siS = si_seq(setdiff(1:11,[1 11 9 2 10])); si = siS; Rs = o.Rs(:,si); mRs = o.mR(:,si);
    props1 = get_props_Rs(Rs,50);
    siD = si_no_brake_dist; si = siD; RsD = o.Rs(:,si); mRsD = o.mR(:,si);
    siT = si_no_brake_time; si = siT; RsT = o.Rs(:,si); mRsT = o.mR(:,si);
    
    propsD = get_props_Rs(RsD,50); propsT = get_props_Rs(RsT,50);
    
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            good_FR = props1.good_FR{rr,cc};
            tzMI = propsD.zMI{rr,cc}; zMID = tzMI(good_FR);
            tzMI = propsT.zMI{rr,cc}; zMIT = tzMI(good_FR);
            dzmi = zMID - zMIT;
            perc_D(rr,cc) = 100*sum(dzmi > 0.5)/length(good_FR);
            perc_T(rr,cc) = 100*sum(dzmi < -0.5)/length(good_FR);
        end
    end
    
    [within,dvn,xlabels] = make_within_table({'DT','Cond','T'},[2,3,2]);
    dataT = make_between_table({perc_D,perc_T},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DT_by_T','hsd'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(6:7); tcolors = repmat(tcolors,2,1);

    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'yspacing',15);

    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    txl = rasterNamesTxt([siD siT]);
    xticks = xdata; xticklabels = {'D-E','T-E','D-E','T-E'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    changePosition(gca,[0.1 -0.02 -0.05 0.01]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('diff_ST_all_pooled.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DT_Cond_T','hsd'},[1 1 1]);
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
    changePosition(gca,[-0.01 -0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
%     save_pdf(hf,mData.pdf_folder,sprintf('diff_ST_all.pdf'),600);
    %%
    break;
end

%% Comparison of the peak locations for zMID > zMIT and zMIT > zMID from distance based rasters
while 1
    %%
    siS = si_seq(setdiff(1:11,[1 11 9 2 10])); si = siS; Rs = o.Rs(:,si); mRs = o.mR(:,si);
    props1 = get_props_Rs(Rs,50);
    siD = si_no_brake_dist; si = siD; RsD = o.Rs(:,si); mRsD = o.mR(:,si);
    siT = si_no_brake_time; si = siT; RsT = o.Rs(:,si); mRsT = o.mR(:,si);
    
    propsD = get_props_Rs(RsD,50); propsT = get_props_Rs(RsT,50);
    
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            good_FR = props1.good_FR{rr,cc};
            tzMI = propsD.zMI{rr,cc}; zMID = tzMI(good_FR);
            tzMI = propsT.zMI{rr,cc}; zMIT = tzMI(good_FR);
            tpeak_loc = propsD.peak_locations{rr,cc}(good_FR);
            peak_locD(rr,cc) = nanmean(tpeak_loc(zMID > zMIT));
%             tpeak_loc = propsT.peak_locations{rr,cc}(good_FR);
            peak_locT(rr,cc) = nanmean(tpeak_loc(zMIT > zMID));
        end
    end
    
    
    [within,dvn,xlabels] = make_within_table({'DT','Cond','T'},[2,3,2]);
    dataT = make_between_table({peak_locD,peak_locT},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DT_by_T','hsd'},[1 1 1]);
    hf = get_figure(5,[8 7 6.99 1.5]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;

    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'yspacing',10);

    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    txl = rasterNamesTxt([siD siT]); 
    xticks = xdata; xticklabels = {'D-E','T-E','D-E','T-E'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 (100)]); xtickangle(45);
    changePosition(gca,[-0.01 -0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. Peak Location (cm)',[0 0 0]});

    %%
    break;
end

%% Comparison of the peak locations for zMID > zMIT and zMIT > zMID from time based rasters
while 1
    %%
    siS = si_seq(setdiff(1:11,[1 11 9 2 10])); si = siS; Rs = o.Rs(:,si); mRs = o.mR(:,si);
    props1 = get_props_Rs(Rs,50);
    siD = si_no_brake_dist; si = siD; RsD = o.Rs(:,si); mRsD = o.mR(:,si);
    siT = si_no_brake_time; si = siT; RsT = o.Rs(:,si); mRsT = o.mR(:,si);
    
    propsD = get_props_Rs(RsD,50); propsT = get_props_Rs(RsT,50);
    
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            good_FR = props1.good_FR{rr,cc};
            tzMI = propsD.zMI{rr,cc}; zMID = tzMI(good_FR);
            tzMI = propsT.zMI{rr,cc}; zMIT = tzMI(good_FR);
            tpeak_loc = propsT.peak_locations{rr,cc}(good_FR);
            peak_locD(rr,cc) = nanmean(tpeak_loc(zMID > zMIT));
%             tpeak_loc = propsT.peak_locations{rr,cc}(good_FR);
            peak_locT(rr,cc) = nanmean(tpeak_loc(zMIT > zMID));
        end
    end
    
    
    [within,dvn,xlabels] = make_within_table({'DT','Cond','T'},[2,3,2]);
    dataT = make_between_table({peak_locD,peak_locT},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DT_by_T','hsd'},[1 1 1]);
    hf = get_figure(6,[8 2 6.99 1.5]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;

    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'yspacing',10);

    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    txl = rasterNamesTxt([siD siT]); 
    xticks = xdata; xticklabels = {'D-E','T-E','D-E','T-E'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 (100)]); xtickangle(45);
    changePosition(gca,[-0.01 -0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. Peak Location (sec)',[0 0 0]});

    %%
    break;
end

%% compare centers dist
while 1
    %%
    siS = si_seq(setdiff(1:11,[1 11 9 2 10])); si = siS; Rs = o.Rs(:,si); mRs = o.mR(:,si);
    props1 = get_props_Rs(Rs,50);
    siD = si_no_brake_dist; si = siD; RsD = o.Rs(:,si); mRsD = o.mR(:,si);
    siT = si_no_brake_time; si = siT; RsT = o.Rs(:,si); mRsT = o.mR(:,si);
    
    propsD = get_props_Rs(RsD,50); propsT = get_props_Rs(RsT,50);
    
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            good_FR = props1.good_FR{rr,cc};
            tzMI = propsD.zMI{rr,cc}; zMID = tzMI(good_FR);
            tzMI = propsT.zMI{rr,cc}; zMIT = tzMI(good_FR);
            tpeak_loc = propsD.MFR{rr,cc}(good_FR);
            peak_locD(rr,cc) = nanmean(tpeak_loc(zMID > zMIT));
%             tpeak_loc = propsT.peak_locations{rr,cc}(good_FR);
            peak_locT(rr,cc) = nanmean(tpeak_loc(zMIT > zMID));
        end
    end
    
    
    [within,dvn,xlabels] = make_within_table({'DT','Cond','T'},[2,3,2]);
    dataT = make_between_table({peak_locD,peak_locT},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'T','hsd'},[1 1 1]);
    hf = get_figure(5,[8 7 6.99 1.5]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;

    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'yspacing',10);

    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    txl = rasterNamesTxt([siD siT]); 
    xticks = xdata; xticklabels = {'Trials','Inter-Trials','Trials','Inter-Trials'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 (100)]); xtickangle(45);
    changePosition(gca,[-0.01 -0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. Peak Location (cm)',[0 0 0]});

    %%
    break;
end

%% compare PWs dist
while 1
    % for trials
    respT1 = dzMIo.resp_D_g_T; respT2 = dzMIo.resp_T_g_D; 
    mZMIsD = (exec_fun_on_cell_mat(reduce_Rs(oDo.props.PWs,respT1),'nanmean'));
    mZMIsDI = (exec_fun_on_cell_mat(reduce_Rs(oDo.props.PWs,respT2),'nanmean'));
    
    respT1 = dzMIoI.resp_D_g_T; respT2 = dzMIoI.resp_T_g_D; 
    mZMIsD1 = (exec_fun_on_cell_mat(reduce_Rs(oIDo.props.PWs,respT1),'nanmean'));
    mZMIsDI1 = (exec_fun_on_cell_mat(reduce_Rs(oIDo.props.PWs,respT2),'nanmean'));
    
    [within,dvn,xlabels] = make_within_table({'TI','P','Cond'},[2,2,3]);
    dataT = make_between_table({mZMIsD,mZMIsDI,mZMIsD1,mZMIsDI1},dvn);
    ra = RMA(dataT,within,0.05);

    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_Cond','bonferroni'},[1 1 1]);
    hf = get_figure(5,[8 7 3.25 5]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = xlabels; set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.35 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Peak Location (cm)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_Comparison_dist_time_peak_locations.pdf'),600);
    ra.ranova
    ra.mauchly
    %%
    break;
end
%% compare MFR dist
while 1
    % for trials
    respT1 = dzMIo.resp_D_g_T; respT2 = dzMIo.resp_T_g_D; 
    mZMIsD = (exec_fun_on_cell_mat(reduce_Rs(oDo.props.MFR,respT1),'nanmean'));
    mZMIsDI = (exec_fun_on_cell_mat(reduce_Rs(oDo.props.MFR,respT2),'nanmean'));
    
    respT1 = dzMIoI.resp_D_g_T; respT2 = dzMIoI.resp_T_g_D; 
    mZMIsD1 = (exec_fun_on_cell_mat(reduce_Rs(oIDo.props.MFR,respT1),'nanmean'));
    mZMIsDI1 = (exec_fun_on_cell_mat(reduce_Rs(oIDo.props.MFR,respT2),'nanmean'));
    
    [within,dvn,xlabels] = make_within_table({'TI','P','Cond'},[2,2,3]);
    dataT = make_between_table({mZMIsD,mZMIsDI,mZMIsD1,mZMIsDI1},dvn);
    ra = RMA(dataT,within,0.05);

    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_Cond','bonferroni'},[1 1 1]);
    hf = get_figure(5,[8 7 3.25 5]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = xlabels; set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.35 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Peak Location (cm)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_Comparison_dist_time_peak_locations.pdf'),600);
    ra.ranova
    ra.mauchly
    %%
    break;
end

%% population vector and correlation  temporal
while 1
    titles = {'Ar-i-T','ArL-i-T','Ar*-i-T'};
    inds = 2:2:6;
    dzMI = prop_op(propsD.zMI,propsT.zMI,0.5,props1S.good_Gauss);
    an = 2;
    Rs1 = RsT(:,inds); mR1 = mRsT(:,inds);
    resp1 = dzMI.resp_D_g_T;
%     resp1 = dzMI.resp_T_g_D;
    resp = resp1(:,inds);

    ff = makeFigureRowsCols(108,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.1],'widthHeightAdjustment',...
        [-35 -60]);    set(gcf,'color','w');    set(gcf,'Position',[7 8 3.45 2]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs1,mR1,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs1(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_temporal_%s.pdf',selected_property),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(109,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.2],'widthHeightAdjustment',...
        [-35 -240]);    set(gcf,'color','w');    set(gcf,'Position',[7 5 3.45 1]);3
    ff = show_population_vector_and_corr(mData,ff,Rs1(an,:),[],aCRc,[],[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_temporal_%s.pdf',selected_property),600);
%     %%
%     break;
% end
% 
% 
% 
% %% population vector and correlation spatial 
% while 1
    titles = {'Ar-t-D','ArL-t-D','Ar*-t-D'};
    inds = 1:2:6;
%     dzMI = prop_op(propsD.zMI,propsT.zMI,0.5,props1S.good_FR);
    an = 4;
    Rs1 = RsD(:,inds); mR1 = mRsD(:,inds);
    resp = resp1(:,inds);
    ff = makeFigureRowsCols(102,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.1],'widthHeightAdjustment',...
        [-35 -60]);    set(gcf,'color','w');    set(gcf,'Position',[3 8 3.45 2]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs1,mR1,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs1(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_spatial_%s.pdf',selected_property),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(103,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.2],'widthHeightAdjustment',...
        [-35 -240]);    set(gcf,'color','w');    set(gcf,'Position',[3 5 3.45 1]);
    ff = show_population_vector_and_corr(mData,ff,Rs1(an,:),[],aCRc,[],[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_spatial_%s.pdf',selected_property),600);
    %%
    break;
end
