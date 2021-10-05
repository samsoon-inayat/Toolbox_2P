selected_property = 'good_FR';
cmdTxt = sprintf('good_FR = props1.%s;',selected_property);

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
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DT_by_T','hsd'},[1 1 1]);
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
%             ttPL_D = propsD.peak_locations{rr,cc}; ttPL_T = propsT.peak_locations{rr,cc};
            diff_zMI(rr,cc) = nanmean(zMID - zMIT);
            dist_diff_zMI{rr,cc} = zMID-zMIT;
        end
        dist_diff_zMI1{rr,1} = [dist_diff_zMI{rr,1};dist_diff_zMI{rr,3};dist_diff_zMI{rr,5}];
        dist_diff_zMI1{rr,2} = [dist_diff_zMI{rr,2};dist_diff_zMI{rr,4};dist_diff_zMI{rr,6}];
    end
    [distDo,allVals] = getAveragesAndAllValues(dist_diff_zMI1);
    minBin = min(allVals);
    maxBin = max(allVals);
    incr = 1;
    tcolors = mData.colors;
    hf = get_figure(8,[5 7 2.25 1.5]);hold on;
    [ha,hb,~,bins] = plotAverageDistributions(dist_diff_zMI1,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'pdf_or_cdf','pdf');
    
    %%
    [within,dvn,xlabels] = make_within_table({'Cond','T'},[3,2]);
    dataT = make_between_table({diff_zMI},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'T','hsd'},[1 1 1]);
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
            perc_zMID_g_zMIT(rr,cc) = 100*sum(zMID - zMIT > 0)/length(zMID);
            perc_zMIT_g_zMID(rr,cc) = 100*sum(zMIT - zMID > 0)/length(zMID);
        end
    end
    
    
    [within,dvn,xlabels] = make_within_table({'DT','Cond','T'},[2,3,2]);
    dataT = make_between_table({perc_zMID_g_zMIT,perc_zMIT_g_zMID},dvn);
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
    xticks = xdata; xticklabels = {'Trials','Inter-Trials','Trials','Inter-Trials'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 (100)]); xtickangle(45);
    changePosition(gca,[-0.01 -0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. zMI Difference',[0 0 0]});

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
    xticks = xdata; xticklabels = {'Trials','Inter-Trials','Trials','Inter-Trials'};
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
    xticks = xdata; xticklabels = {'Trials','Inter-Trials','Trials','Inter-Trials'};
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
