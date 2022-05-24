function air_place_cells_figure
%% Load Data
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei = evalin('base','ei'); 
    
    selContexts = [3 3 4 4 5 5]; rasterNames = {'airD','airID','airD','airID','airD','airID'};
    oD_seq = get_data(ei,selContexts,rasterNames);
    
    selContexts = [3 3 4 4 5 5]; rasterNames = {'airT','airIT','airT','airIT','airT','airIT'};
    oT_seq = get_data(ei,selContexts,rasterNames);
    
    selContexts = [3 4 5]; rasterNames = {'air77T','air77T','air77T'};
    o77T = get_data(ei,selContexts,rasterNames);

    selContexts = [3 4 5 3 4 5]; rasterNames = {'airT','airT','airT','airIT','airIT','airIT'};
    oT = get_data(ei,selContexts,rasterNames);

    selContexts = [3 4 5 3 4 5]; rasterNames = {'airD','airD','airD','airID','airID','airID'};
    oD = get_data(ei,selContexts,rasterNames);

    selContexts = [3 4 5 3 4 5]; rasterNames = {'airT','airT','airT','airID','airID','airID'};
    oTD = get_data(ei,selContexts,rasterNames);

    selContexts = [3 3 4 4 5 5]; rasterNames = {'airD','airIT','airD','airIT','airD','airIT'};
    oDT = get_data(ei,selContexts,rasterNames);

%     dzMI = prop_op(oD.props.zMI,oT.props.zMI,0.25);
    
    selContexts = [3 4 5]; rasterNames = {'airT','airT','airT'}; oTo = get_data(ei,selContexts,rasterNames);

    selContexts = [3 4 5]; rasterNames = {'airD','airD','airD'}; oDo = get_data(ei,selContexts,rasterNames);
    
    selContexts = [3 4 5]; rasterNames = {'airIT','airIT','airIT'}; oITo = get_data(ei,selContexts,rasterNames);

    selContexts = [3 4 5]; rasterNames = {'airID','airID','airID'}; oIDo = get_data(ei,selContexts,rasterNames);
    
%     dzMIo = prop_op(oDo.props.zMI,oTo.props.zMI,0.1);
%     dzMIoI = prop_op(oIDo.props.zMI,oITo.props.zMI,0.1);
%     dzMI_seq = prop_op(oD_seq.props.zMI,oT_seq.props.zMI,0.1);
%     
%     pop_DgT_T = dzMIo.resp_D_g_T; 
%     pop_TgD_T = dzMIo.resp_T_g_D; 
%     
%     pop_DgT_IT = dzMIoI.resp_D_g_T;
%     pop_TgD_IT = dzMIoI.resp_T_g_D;
%     
%     pop_DgT_and_good_Gauss_T = cell_list_op(pop_DgT_T,oDo.props.good_Gauss,'and');
%     pop_TgD_and_good_Gauss_T = cell_list_op(pop_TgD_T,oTo.props.good_Gauss,'and');
%     
%     pop_DgT_and_good_Gauss_IT = cell_list_op(pop_DgT_IT,oIDo.props.good_Gauss,'and');
%     pop_TgD_and_good_Gauss_IT = cell_list_op(pop_TgD_IT,oITo.props.good_Gauss,'and');
    
    
    break
end
n = 0;
%% find overlap between the two types of populations (zMID > zMIT and zMIT > zMID) across conditions and trials-inter/trials
while 1
    [OI,mOI,semOI] = get_overlap_index([pop_DgT_T pop_TgD_T pop_DgT_IT pop_TgD_IT]);
    an = 1;
    figure(100);clf;imagesc(mOI,[0 1]);colorbar

    [OI,mOI,semOI] = get_overlap_index([pop_DgT_and_good_Gauss_T pop_TgD_and_good_Gauss_T pop_DgT_and_good_Gauss_IT pop_TgD_and_good_Gauss_IT]);
    an = 1;
    figure(101);clf;imagesc(mOI,[0 1]);colorbar
    break;
end
%% plot the dzMI values sequentially and run anova on the mean
while 1
    more_dist = dzMI_seq.diff_D_T;
%     popD = cell_list_op(dzMI_seq.resp_D_g_T,oD_seq.props.good_Gauss,'and');
%     popT = cell_list_op(dzMI_seq.resp_T_g_D,oT_seq.props.good_Gauss,'and');
    popD = dzMI_seq.resp_D_g_T; popT = dzMI_seq.resp_T_g_D;
    for rr = 1:size(more_dist,1)
        vals = [];
        for cc = 1:size(more_dist,2)
            vals(:,cc) = more_dist{rr,cc};
            polarity(:,cc) = popD{rr,cc} - popT{rr,cc};
        end
        m_vals(rr,:) = nanmean(vals);
        D_counts = sum(polarity == 1,2);
        T_counts = sum(polarity == -1,2);
        dt_ind = D_counts
    end
    [within,dvn,xlabels] = make_within_table({'C','T'},[3,2]);
    dataT = make_between_table({m_vals},dvn);
    ra = RMA(dataT,within,0.05);
    ra.ranova
    ra.mauchly
%%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'C_by_T','bonferroni'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = xlabels; set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.09 0.02 -0.15 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'zMID-zMIT',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('dzMI_Trials_IT.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'T','bonferroni'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Trials','I-Trials'}; set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.09 0.02 -0.55 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'zMID-zMIT',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('dzMI_Trials_IT_pooled.pdf'),600);
    %%
    break;
end
%% compare raw percentages of cells where zMID > zMIT for both trials and inter-trials
while 1
    [within,dvn,xlabels] = make_within_table({'T','P','C'},[2,2,3]);
    TPC = [find_percent(pop_DgT_T) find_percent(pop_TgD_T)] ;
    ITPC = [find_percent(pop_DgT_IT) find_percent(pop_TgD_IT)] ;
    dataT = make_between_table({TPC,ITPC},dvn);
    ra = RMA(dataT,within,0.05);
%     writetable(dataT,'spss_three_way_anova.xls');
%%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'T_by_P','bonferroni'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Tr-DgT','Tr-TgD','ITr-DgT','Itr-TgD'}; set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.35 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('percentages_zMI_Comparison_dist_time.pdf'),600);
    ra.ranova
    ra.mauchly
    %%
    break;
end
%% compare percentages - good Gauss of cells where zMID > zMIT for both trials and inter-trials
while 1
    [within,dvn,xlabels] = make_within_table({'T','P','C'},[2,2,3]);
    TPC = [find_percent(pop_DgT_and_good_Gauss_T) find_percent(pop_TgD_and_good_Gauss_T)] ;
    ITPC = [find_percent(pop_DgT_and_good_Gauss_IT) find_percent(pop_TgD_and_good_Gauss_IT)] ;
    dataT = make_between_table({TPC,ITPC},dvn);
    ra = RMA(dataT,within,0.05);
%%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'T_by_P','bonferroni'},[1 1 1]);
    hf = get_figure(6,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Tr-DgT','Tr-TgD','ITr-DgT','Itr-TgD'}; set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.35 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('percentages_zMI_Comparison_dist_time_good_Gauss.pdf'),600);
    ra.ranova
    ra.mauchly
    %%
    break;
end
%% compare zMIsD of the populations with zMID > zMIT for both trials and inter-trials
while 1
    % for trials
    props = oDo.props; propsI = oIDo.props;
    respT = dzMIo.resp_D_g_T; respTI = dzMIoI.resp_D_g_T;
    mZMIsD = exec_fun_on_cell_mat(reduce_Rs(props.zMI,respT),'nanmean');
    
    mZMIsDI = exec_fun_on_cell_mat(reduce_Rs(propsI.zMI,respTI),'nanmean');
   
    [within,dvn,xlabels] = make_within_table({'T','C'},[2,3]);
    dataT = make_between_table({mZMIsD,mZMIsDI},dvn);
    ra = RMA(dataT,within,0.05);
    ra.ranova
    ra.mauchly
%%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'T_by_C','bonferroni'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = xlabels; set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('zMID_Comparison_dist_time_zMID_g_zMIT.pdf'),600);
    
     [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'T','bonferroni'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Tr','ITr'}; set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.55 -0.05]);% put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('zMID_Comparison_dist_time_zMID_g_zMIT_pooled.pdf'),600);
    %%
    break;
end

%% compare zMIsT of the populations with zMIT > zMID for both trials and inter-trials
while 1
    % for trials
    props = oTo.props; propsI = oITo.props;
    respT = dzMIo.resp_T_g_D; respTI = dzMIoI.resp_T_g_D;
    mZMIsD = exec_fun_on_cell_mat(reduce_Rs(props.zMI,respT),'nanmean');
    
    mZMIsDI = exec_fun_on_cell_mat(reduce_Rs(propsI.zMI,respTI),'nanmean');
   
    [within,dvn,xlabels] = make_within_table({'T','C'},[2,3]);
    dataT = make_between_table({mZMIsD,mZMIsDI},dvn);
    ra = RMA(dataT,within,0.05);
    ra.ranova
    ra.mauchly
%%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'T_by_C','bonferroni'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = xlabels; set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('zMIT_Comparison_dist_time_zMIT_g_zMID.pdf'),600);
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'T','bonferroni'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Tr','ITr'}; set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.55 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('zMIT_Comparison_dist_time_zMIT_g_zMID_pooled.pdf'),600);
    
    %%
    break;
end

%% compare peak locations in distance rasters of the populations with zMIT > zMID and zMID > zMIT for both trials and inter-trials and conditions
while 1
     % for trials
    respT = dzMIo.resp_D_g_T; respTI = dzMIoI.resp_D_g_T;
    respT1 = dzMIo.resp_T_g_D; respTI1 = dzMIoI.resp_T_g_D;
    
%     respT = pop_DgT_and_good_Gauss_T; respTI = pop_DgT_and_good_Gauss_IT;
%     respT1 = pop_TgD_and_good_Gauss_T; respTI1 = pop_TgD_and_good_Gauss_IT;
    
    props = oDo.props; propsI = oIDo.props;
    mZMIsD = exec_fun_on_cell_mat(reduce_Rs(props.peak_locations,respT),'nanmean');
    mZMIsD1 = exec_fun_on_cell_mat(reduce_Rs(props.peak_locations,respT1),'nanmean');
    
    mZMIsDI = exec_fun_on_cell_mat(reduce_Rs(propsI.peak_locations,respTI),'nanmean');
    mZMIsDI1 = exec_fun_on_cell_mat(reduce_Rs(propsI.peak_locations,respTI1),'nanmean');
   
    [within,dvn,xlabels] = make_within_table({'T','DT','C'},[2,2,3]);
    dataT = make_between_table({mZMIsD,mZMIsD1,mZMIsDI,mZMIsDI1},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    ra.mauchly
%%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'T_by_DT','hsd'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Tr-Dist','Tr-Time','ITr-Dist','ITr-Time'}; set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.1 0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'cm',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('peak_locations_Comparison_dist_time_zMID_g_zMIT.pdf'),600);
    %%
    break;
end

%% compare peak locations in time rasters of the populations with zMIT > zMID and zMID > zMIT for both trials and inter-trials and conditions
while 1
     % for trials
    respT = dzMIo.resp_D_g_T; respTI = dzMIoI.resp_D_g_T;
    respT1 = dzMIo.resp_T_g_D; respTI1 = dzMIoI.resp_T_g_D;
    props = oTo.props; propsI = oITo.props;
    mZMIsD = exec_fun_on_cell_mat(reduce_Rs(props.peak_locations,respT),'nanmean');
    mZMIsD1 = exec_fun_on_cell_mat(reduce_Rs(props.peak_locations,respT1),'nanmean');
    
    mZMIsDI = exec_fun_on_cell_mat(reduce_Rs(propsI.peak_locations,respTI),'nanmean');
    mZMIsDI1 = exec_fun_on_cell_mat(reduce_Rs(propsI.peak_locations,respTI1),'nanmean');
   
    [within,dvn,xlabels] = make_within_table({'T','DT','C'},[2,2,3]);
    dataT = make_between_table({mZMIsD,mZMIsD1,mZMIsDI,mZMIsDI1},dvn);
    ra = RMA(dataT,within);
    ra.ranova
%     ra.mauchly
%%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'T_by_DT','hsd'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Tr-Dist','Tr-Time','ITr-Dist','ITr-Time'}; set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.1 0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'sec',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('peak_locationsT_Comparison_dist_time_zMID_g_zMIT.pdf'),600);
    %%
    break;
end
%% compare centers dist
while 1
    % for trials
    respT1 = dzMIo.resp_D_g_T; respT2 = dzMIo.resp_T_g_D; 
    mZMIsD = (exec_fun_on_cell_mat(reduce_Rs(oDo.props.centers,respT1),'nanmean'));
    mZMIsDI = (exec_fun_on_cell_mat(reduce_Rs(oDo.props.centers,respT2),'nanmean'));
    
    respT1 = dzMIoI.resp_D_g_T; respT2 = dzMIoI.resp_T_g_D; 
    mZMIsD1 = (exec_fun_on_cell_mat(reduce_Rs(oIDo.props.centers,respT1),'nanmean'));
    mZMIsDI1 = (exec_fun_on_cell_mat(reduce_Rs(oIDo.props.centers,respT2),'nanmean'));
    
    [within,dvn,xlabels] = make_within_table({'TI','P','Cond'},[2,2,3]);
    dataT = make_between_table({mZMIsD,mZMIsDI,mZMIsD1,mZMIsDI1},dvn);
    ra = RMA(dataT,within,0.05);

    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_P','bonferroni'},[1 1 1]);
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

%% Population vectors and correlation
while 1
    %% population vector and correlation single animal
    Rs = oD_seq.Rs; mR = oD_seq.mR; resp = dzMI_seq.resp_D_g_T;%   resp = cell_list_op(dzMIo.resp_T_g_D,[],'or');
    an = 3; R = Rs{an,1};
    ff = makeFigureRowsCols(106,[1 0.5 6 1],'RowsCols',[2 6],...
        'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.07 0.1],'widthHeightAdjustment',...
        [0.01 -60]);
    set(gcf,'color','w');
    set(gcf,'Position',[1 7 6.99 2]);
%     [resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC,resp_exc_inh] = get_responsive_fraction(Rs);
%     resp = get_cell_list(resp_valsC,[]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[-0.1 1],[]);
    set_obj(ff,{'FontWeight','Normal','FontSize',6,'LineWidth',0.25});
    ht = get_obj(ff,'title'); hyl = get_obj(ff,'ylabel'); changePosition(hyl(1,1),[4 0 0]);
    set_obj(ht,'String',repmat({'Pop. Activity';'Pop. Correlation'},1,6));
    set_obj(ht,{'FontSize',5,'FontWeight','Normal'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('population_vector_corr_zMID_g_zMIT_D.pdf'),600);
    
    Rs = oT_seq.Rs; mR = oT_seq.mR; resp = dzMI_seq.resp_D_g_T;%   resp = cell_list_op(dzMIo.resp_T_g_D,[],'or');
    R = Rs{an,1};
    ff = makeFigureRowsCols(107,[1 0.5 6 1],'RowsCols',[2 6],...
        'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.07 0.1],'widthHeightAdjustment',...
        [0.01 -60]);
    set(gcf,'color','w');
    set(gcf,'Position',[1 3 6.99 2]);
%     [resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC,resp_exc_inh] = get_responsive_fraction(Rs);
%     resp = get_cell_list(resp_valsC,[]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[-0.1 1],[]);
    set_obj(ff,{'FontWeight','Normal','FontSize',6,'LineWidth',0.25});
    ht = get_obj(ff,'title'); hyl = get_obj(ff,'ylabel'); changePosition(hyl(1,1),[4 0 0]);
    set_obj(ht,'String',repmat({'Pop. Activity';'Pop. Correlation'},1,6));
    set_obj(ht,{'FontSize',5,'FontWeight','Normal'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('population_vector_corr_zMID_g_zMIT_T.pdf'),600);

    
    Rs = oD_seq.Rs; mR = oD_seq.mR; resp = dzMI_seq.resp_T_g_D;%   resp = cell_list_op(dzMIo.resp_T_g_D,[],'or');
     R = Rs{an,1};
    ff = makeFigureRowsCols(108,[1 0.5 6 1],'RowsCols',[2 6],...
        'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.07 0.1],'widthHeightAdjustment',...
        [0.01 -60]);
    set(gcf,'color','w');
    set(gcf,'Position',[10 7 6.99 2]);
%     [resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC,resp_exc_inh] = get_responsive_fraction(Rs);
%     resp = get_cell_list(resp_valsC,[]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[-0.1 1],[]);
    set_obj(ff,{'FontWeight','Normal','FontSize',6,'LineWidth',0.25});
    ht = get_obj(ff,'title'); hyl = get_obj(ff,'ylabel'); changePosition(hyl(1,1),[4 0 0]);
    set_obj(ht,'String',repmat({'Pop. Activity';'Pop. Correlation'},1,6));
    set_obj(ht,{'FontSize',5,'FontWeight','Normal'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('population_vector_corr_zMIT_g_zMID_D.pdf'),600);
    
    Rs = oT_seq.Rs; mR = oT_seq.mR; resp = dzMI_seq.resp_T_g_D;%   resp = cell_list_op(dzMIo.resp_T_g_D,[],'or');
    R = Rs{an,1};
    ff = makeFigureRowsCols(109,[1 0.5 6 1],'RowsCols',[2 6],...
        'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.07 0.1],'widthHeightAdjustment',...
        [0.01 -60]);
    set(gcf,'color','w');
    set(gcf,'Position',[10 3 6.99 2]);
%     [resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC,resp_exc_inh] = get_responsive_fraction(Rs);
%     resp = get_cell_list(resp_valsC,[]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[-0.1 1],[]);
    set_obj(ff,{'FontWeight','Normal','FontSize',6,'LineWidth',0.25});
    ht = get_obj(ff,'title'); hyl = get_obj(ff,'ylabel'); changePosition(hyl(1,1),[4 0 0]);
    set_obj(ht,'String',repmat({'Pop. Activity';'Pop. Correlation'},1,6));
    set_obj(ht,{'FontSize',5,'FontWeight','Normal'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('population_vector_corr_zMIT_g_zMID_T.pdf'),600);

%% average correlation of all animals
    ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 6],...
        'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.075 0.2],'widthHeightAdjustment',...
        [0.01 -250]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 5 6.99 1]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[-0.1 1],[]);
    set_obj(ff,{'FontWeight','Normal','FontSize',6,'LineWidth',0.25});
    ht = get_obj(ff,'title'); 
    set_obj(ht,'String',{'Avg. Pop. Correlation','Avg. Pop. Correlation','Avg. Pop. Correlation'});
    set_obj(ht,{'FontSize',5,'FontWeight','Normal'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_average_population_vector_corrD.pdf'),600);
%%
    break;
end
%% compare zMIs overall (for all cells)
while 1
   
    mZMIsD = exec_fun_on_cell_mat(oDo.props.zMI,'nanmean'); 
    mZMIsT = exec_fun_on_cell_mat(oTo.props.zMI,'nanmean');
    mZMIsD1 = exec_fun_on_cell_mat(oIDo.props.zMI,'nanmean'); 
    mZMIsT1 = exec_fun_on_cell_mat(oITo.props.zMI,'nanmean');
    
    [within,dvn,xlabels] = make_within_table({'T','DT','Cond'},[2,2,3]);
    dataT = make_between_table({mZMIsD,mZMIsD1,mZMIsT,mZMIsT1},dvn);
    ra = RMA(dataT,within,0.05);
    ra.ranova
    ra.mauchly
%%
   [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'T_by_DT','bonferroni'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Tr-Dist','Tr-Time','ITr-Dist','ITr-Time'}; set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.1 0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'zMI',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_overall.pdf'),600);
    %%
    break;
end
%% compare PWs
while 1
%     oT = oIDo; dzMIT = dzMIoI;
    oT = oDo; respT1 = dzMIo.resp_D_g_T; respT2 = dzMIo.resp_T_g_D; 
    popDT = cell_list_op(respT1,oT.props.good,'and'); popTD = cell_list_op(respT2,oT.props.good,'and');
    zMIsD = reduce_Rs(oT.props.MFR,popDT); zMIsT = reduce_Rs(oT.props.MFR,popTD);
    
    mZMIsD = exec_fun_on_cell_mat(zMIsD,'nanmean'); mZMIsT = exec_fun_on_cell_mat(zMIsT,'nanmean');
    
    [within,dvn,xlabels] = make_within_table({'DT','Cond'},[2,3]);
    dataT = make_between_table({mZMIsD,mZMIsT},dvn);
    ra = RMA(dataT,within,0.05);

    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'DT','bonferroni'},[1 1 1]);
    hf = get_figure(5,[8 7 3.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'D-T','D-IT','T-T','T-IT'}; set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.35 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'zMI',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('avg_speed_anova.pdf'),600);
    ra.ranova
    ra.mauchly
    %%
    break;
end
%% find correlations across conditions
resp = resp_ORCS;
[CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,1,0);

outRemap = find_population_vector_corr_remap(Rs,mR,resp_ORCS);
%% find trial by trial comparison - correlations
trials = mat2cell([1:10]',ones(size([1:10]')));
resp = get_cell_list(resp_valsCS,[]);
out1 = find_population_vector_corr_remap_trials(Rs(:,1),resp,trials);
out2 = find_population_vector_corr_remap_trials(Rs(:,2),resp,trials);
out3 = find_population_vector_corr_remap_trials(Rs(:,3),resp,trials);
% save(fullfile(mData.pd_folder,'air_place_cells_figure.mat'),'out1','out2','out3');
n = 0;
%% Speed Figure
while 1
%     Rs = oDT.Rs;
    for ii = 1:length(ei)
        b1 = ei{ii}.b;
        for jj = 1:10
            alds(ii,jj) = b1.dist(b1.stim_r(jj+10)) - b1.dist(b1.air_puff_r(jj+20));
        end
    end
    ald = round(mean(alds(:)));
     ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 6],'spaceRowsCols',[0 0],'rightUpShifts',[-0.05 0.13],'widthHeightAdjustment',[-3 -350]);
    set(gcf,'color','w'); set(gcf,'Position',[5 5 5.60 1]);
    [Y,E] = discretize(1:49,3);
    all_speeds = []; cTxt = {'Ar-Trials','Ar-Inter-Trials','ArL-Trials','ArL-Inter-Trials','Ar*-Trials','Ar*-Inter-Trials'}; 
    for cn = 1:6
        mean_speed_over_trials = [];
        aThisSpeed = [];
        for an = 1:size(Rs,1)
            thisSpeed = nanmean(Rs{an,cn}.speed);
            mean_speed_over_trials(an,:) = thisSpeed;
            for bb = 1:3
                aThisSpeed(an,bb) = nanmean(thisSpeed(Y==bb));
            end
        end
        all_speeds = [all_speeds aThisSpeed];
        axes(ff.h_axes(1,cn));
        hold on;
        xs = Rs{1,cn}.xs; N = length(xs);
        mspeed = mean(mean_speed_over_trials(:,1:N)); semspeed = std(mean_speed_over_trials(:,1:N))/sqrt(5);
        plot(xs,mspeed);
        shadedErrorBar(xs,mspeed,semspeed);
        changePosition(gca,[0.1 0.15 -0.05 -0.15]);
        if cn == 1
            put_axes_labels(gca,{'',[0 0 0]},{{'Speed (cm/sec)'},[0 -1 0]});
        end
        if ismember(cn,[1 3 5])
            if cn == 3
                plot([ald ald],[0 30],'m:','linewidth',0.5);
            end
            plot([50 50],[0 30],'k--','linewidth',0.25);
            plot([100 100],[0 30],'k--','linewidth',0.25);
            bTxt = {'dB1','dB2','dB3'}; 
            xbTxt = [25 75 125]-7; ybTxt = 31;
            for ii = 1:length(bTxt)
                text(xbTxt(ii),ybTxt,bTxt{ii},'FontSize',5);
            end
            xlabel('Distance (cm)');
        else
            plot([5 5],[0 30],'k--','linewidth',0.25);
            plot([10 10],[0 30],'k--','linewidth',0.25);
            bTxt = {'tB1','tB2','tB3'}; xbTxt = [2.5 7.5 12.5]-1; ybTxt = 31;
            for ii = 1:length(bTxt)
                text(xbTxt(ii),ybTxt,bTxt{ii},'FontSize',5);
            end
            xlabel('Time (sec)');
        end
        text(xbTxt(1),ybTxt+10,cTxt{cn},'FontSize',5);
        ylim([0 30]);
        box off;
        format_axes(gca);
    end
    
    
    
    save_pdf(ff.hf,mData.pdf_folder,sprintf('speeds_345'),600);
%%
    all_speeds_Trials = all_speeds(:,[1:3 7:9 13:15]+0);
    [within,dvn,xlabels] = make_within_table({'Conds','Bins'},[3,3]);
    dataT = make_between_table({all_speeds_Trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Conds_by_Bins','hsd'},[1 0.5 1]);
    h(h==1) = 0;
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors(7:9); tcolors = repmat(tcolors,3,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 37]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[-0.01 0.02 0.05 -0.05])
    put_axes_labels(gca,{[],[0 0 0]},{'Avg. Speed (cm/sec)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('avg_speed_anova.pdf'),600);
    ra.ranova
    ra.mauchly
    %%
    all_speeds_Trials = all_speeds(:,[1:3 7:9 13:15]+0);
    [within,dvn,xlabels] = make_within_table({'Conds','Bins'},[3,3]);
    dataT = make_between_table({all_speeds_Trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Bins','hsd'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors(7:9);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 37]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'yticklabels',[]); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.35 -0.05])
%     put_axes_labels(gca,{[],[0 0 0]},{'Avg. Speed (cm/sec)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('avg_speed_anova_pooled.pdf'),600);
    %%
    all_speeds_Trials = all_speeds(:,[7:9]+0);
    [within,dvn,xlabels] = make_within_table({'Bins'},[3]);
    dataT = make_between_table({all_speeds_Trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Bins','hsd'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors(7:9);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 37]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'yticklabels',[]); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.35 -0.05])
%     put_axes_labels(gca,{[],[0 0 0]},{'Avg. Speed (cm/sec)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('avg_speed_anova_pooled.pdf'),600);
    %%
    all_speeds_Trials = all_speeds(:,[1:3 7:9 13:15]+3);
    [within,dvn,xlabels] = make_within_table({'Conds','Bins'},[3,3]);
    dataT = make_between_table({all_speeds_Trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Bins','hsd'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors([10 11 12]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 37]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'tB1','tB2','tB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.35 -0.05])
    put_axes_labels(gca,{[],[0 0 0]},{'Avg. Speed (cm/sec)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('avg_speed_anova_pooled_IT.pdf'),600);
    %%
    break;
end
%% Show rasters to view
an = 3; cn = 1;
Rs = oD_seq.Rs; RsT = oT_seq.Rs; Rs77T = o77T.Rs;
% ccsT = cell_list_op(dzMI_seq.resp_T_g_D,oT_seq.props.good_Gauss,'and'); 
% ccsD = cell_list_op(dzMI_seq.resp_D_g_T,oD_seq.props.good_Gauss,'and'); 
ccsT = cell_list_op(dzMI_seq.resp_T_g_D,oT_seq.props.good_Gauss,'and'); 
ccsD = cell_list_op(dzMI_seq.resp_D_g_T,oD_seq.props.good_zMI_FR,'and'); 
% ccsD = dzMI_seq.resp_D_g_T;
% ccsT = dzMI_seq.resp_T_g_D;
ccs = ccsT;
ccs = o77T.resp.vals;
% plotRasters_simplest(Rs{an,cn},find(ccs{an,cn}))
plotRasters_multi_simplest_line({Rs77T{an,cn}.sp_rasters,RsT{an,cn}.sp_rasters},find(ccs{an,cn}))
%% Show sample rasters
while 1
    Rs = oDT.Rs;
    props1 = get_props_Rs(Rs,5); 
    an = 2; cn = 1;
    good_FR = props1.good_FR(an,[1:2:6]);
    gfr = cell_list_op(good_FR,[],'and');
    % plotRasters_simplest(Rs{an,cn})
    % find(resp_valsC{an}(:,cn));
    ff = makeFigureRowsCols(2020,[10 4 5.1 2],'RowsCols',[2 6],'spaceRowsCols',[0.15 0.05],'rightUpShifts',[0.05 0.15],'widthHeightAdjustment',[-55 -230]);
%     ff = sample_rasters(Rs(an,:),[24 96 41],ff);
    ff = sample_rasters(Rs(an,:),[512 191],ff);%ff = sample_rasters(Rs(an,:),[512 451 191],ff);
    tff = ff; tff.h_axes = ff.h_axes(1,:); xls = get_obj(tff,'xlabel'); set_obj(xls,'string',''); 
    cTxt = {'Ar-t-D','Ar-i-T','ArL-t-D','ArL-i-T','Ar*-t-D','Ar*-i-T'}; 
    for ii = 1:length(tff.h_axes)
        axes(tff.h_axes(ii));
        ylims = ylim;
        text(12,14,cTxt{ii},'FontSize',8,'FontWeight','Bold');
    end
%     tff = ff; tff.h_axes = ff.h_axes(2,:); xls = get_obj(tff,'xlabel'); set_obj(xls,'string','');
    for ii = 1:2:6
        tff = ff; tff.h_axes = ff.h_axes(:,ii); xls = get_obj(tff,'xtick'); set_obj(tff,{'xtick',[1 25 49],'xticklabel',{'0','75','150'}});
    end
    for ii = 2:2:6
        tff = ff; tff.h_axes = ff.h_axes(:,ii); xls = get_obj(tff,'xtick'); set_obj(tff,{'xtick',[1 68 136],'xticklabel',{'0','7.5','15'}});
    end
    tff = ff; tff.h_axes = ff.h_axes(1,2:end); xls = get_obj(tff,'ylabel'); set_obj(xls,'string',''); 
    tff = ff; tff.h_axes = ff.h_axes(2,2:end); xls = get_obj(tff,'ylabel'); set_obj(xls,'string',''); 
%     tff = ff; tff.h_axes = ff.h_axes(3,2:end); xls = get_obj(tff,'ylabel'); set_obj(xls,'string',''); 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersDIT'),600);
    colormap parula
    break;
end

%% Show sample rasters for presentation
while 1
    Rs = oDT.Rs;
    props1 = get_props_Rs(Rs,5); 
    an = 4; cn = 1;
    good_FR = props1.good_FR_and_zMI(an,[1:2:6]);
%     gfr = cell_list_op(good_FR,[],'and');
    % plotRasters_simplest(Rs{an,cn})
    % find(resp_valsC{an}(:,cn));
    ff = makeFigureRowsCols(2020,[10 4 6.1 4],'RowsCols',[4 6],'spaceRowsCols',[0.15 0.05],'rightUpShifts',[0.05 0.15],'widthHeightAdjustment',[-55 -130]);
%     ff = sample_rasters(Rs(an,:),[24 96 41],ff);
    ff = sample_rasters(Rs(an,:),[23 35 46 62],ff);%ff = sample_rasters(Rs(an,:),[512 451 191],ff);
    tff = ff; tff.h_axes = ff.h_axes(1,:); xls = get_obj(tff,'xlabel'); set_obj(xls,'string',''); 
    cTxt = {'Ar-t-D','Ar-i-T','ArL-t-D','ArL-i-T','Ar*-t-D','Ar*-i-T'}; 
    for ii = 1:length(tff.h_axes)
        axes(tff.h_axes(ii));
        ylims = ylim;
        text(12,14,cTxt{ii},'FontSize',8,'FontWeight','Bold');
    end
%     tff = ff; tff.h_axes = ff.h_axes(2,:); xls = get_obj(tff,'xlabel'); set_obj(xls,'string','');
    for ii = 1:2:6
        tff = ff; tff.h_axes = ff.h_axes(:,ii); xls = get_obj(tff,'xtick'); set_obj(tff,{'xtick',[1 25 49],'xticklabel',{'0','75','150'}});
    end
    for ii = 2:2:6
        tff = ff; tff.h_axes = ff.h_axes(:,ii); xls = get_obj(tff,'xtick'); set_obj(tff,{'xtick',[1 68 136],'xticklabel',{'0','7.5','15'}});
    end
    tff = ff; tff.h_axes = ff.h_axes(1,2:end); xls = get_obj(tff,'ylabel'); set_obj(xls,'string',''); 
    tff = ff; tff.h_axes = ff.h_axes(2,2:end); xls = get_obj(tff,'ylabel'); set_obj(xls,'string',''); 
%     tff = ff; tff.h_axes = ff.h_axes(3,2:end); xls = get_obj(tff,'ylabel'); set_obj(xls,'string',''); 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersDIT'),600);
    colormap parula
    break;
end
%% Mutual Information Time versus Distance Distributions
while 1
    Rs = oD.Rs; RsT = oT.Rs;
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            R = Rs{rr,cc};
            zMIsC{rr,cc} = R.info_metrics.ShannonMI_Zsh';
            R = RsT{rr,cc};
            zMIsA{rr,cc} = R.info_metrics.ShannonMI_Zsh';
        end
    end
    ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 6],'spaceRowsCols',[0 0],'rightUpShifts',[-0.05 0.13],'widthHeightAdjustment',[-3 -350]);
    set(gcf,'color','w'); set(gcf,'Position',[5 5 6.99 1]);
    [Y,E] = discretize(1:49,3);
    all_speeds = []; bTxt = {'B1','B2','B3'}; cTxt = {'C3-Trial','C3-Inter-Trial','C4-Trial','C4-Inter-Trial','C3''-Trial','C3''-Inter-Trial'}; 
    for cn = 1:6
        CN = cn;
        axes(ff.h_axes(1,cn));
        hold on;
        tcolors = {'m','c'};
        distD(:,1) = zMIsC(:,CN);
        distD(:,2) = zMIsA(:,CN);
        [distDo,allVals] = getAveragesAndAllValues(distD);
        minBin = min(allVals);
        maxBin = max(allVals);
        incr = 0.001; 
        [ha,hb,hca] = plotDistributions(distD,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
%         plot([1.65 1.65],[0 100],'--k');
%         myxlims2 = [7 9 13];
%         xlims = xlim; xlim([xlims(1) myxlims2(CN)]);
        set(gca,'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
        changePosition(gca,[0.1 0.15 -0.05 -0.15]);
        if CN == 1
            put_axes_labels(gca,{'zMI',[0 0 0]},{{'Cells (%)'},[0 0 0]});
        else
            put_axes_labels(gca,{'zMI',[0 0 0]},{'',[1 0 0]});
        end
        if CN == 1
            legs = {'Dist (D)','Time (T)'};
            legs{end+1} = [4 2 60 5];
            putLegend(gca,legs,'colors',tcolors)
        end
    end
    

    save_pdf(ff.hf,mData.pdf_folder,sprintf('Distribution_zMI_%d',CN),600);
    break;
end

%% Mutual Information Distance-Time (minus) Distributions
while 1
    distD = [];
    Rs = oD.Rs; RsT = oT.Rs;
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            R = Rs{rr,cc};
            zMIsC{rr,cc} = R.info_metrics.ShannonMI_Zsh';
            R = RsT{rr,cc};
            zMIsA{rr,cc} = R.info_metrics.ShannonMI_Zsh';
            dzMI_raw{rr,cc} = zMIsC{rr,cc} - zMIsA{rr,cc};
        end
    end
    ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 6],'spaceRowsCols',[0 0],'rightUpShifts',[-0.05 0.13],'widthHeightAdjustment',[-3 -350]);
    set(gcf,'color','w'); set(gcf,'Position',[5 5 6.99 1]);
    [Y,E] = discretize(1:49,3);
    all_speeds = []; bTxt = {'B1','B2','B3'}; cTxt = {'C3-Trial','C3-Inter-Trial','C4-Trial','C4-Inter-Trial','C3''-Trial','C3''-Inter-Trial'}; 
    for cn = 1:6
        CN = cn;
        axes(ff.h_axes(1,cn));
        hold on;
        tcolors = {'k'};
        distD = dzMI_raw(:,CN);
        [distDo,allVals] = getAveragesAndAllValues(distD);
        minBin = min(allVals);
        maxBin = max(allVals);
        incr = 0.5; 
        [ha,hb,hca] = plotAverageDistributions(distD,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'pdf_or_cdf','pdf');
%         plot([1.65 1.65],[0 100],'--k');
%         myxlims2 = [7 9 13];
%         xlims = xlim; xlim([xlims(1) myxlims2(CN)]);
        set(gca,'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
        changePosition(gca,[0.1 0.15 -0.05 -0.15]);
        if CN == 1
            put_axes_labels(gca,{'zMI',[0 0 0]},{{'Cells (%)'},[0 0 0]});
        else
            put_axes_labels(gca,{'zMI',[0 0 0]},{'',[1 0 0]});
        end
%         if CN == 1
%             legs = {'Dist-Time'};
%             legs{end+1} = [4 2 60 5];
%             putLegend(gca,legs,'colors',tcolors)
%         end
    end
    

    save_pdf(ff.hf,mData.pdf_folder,sprintf('Distribution_diff_zMI'),600);
    break;
end
%% Mutual Information Time versus Distance bar graph
while 1
    mzMIsC = arrayfun(@(x) nanmean(x{1}),zMIsC);
    mzMIsA = arrayfun(@(x) nanmean(x{1}),zMIsA);
    [within,dvn,xlabels] = make_within_table({'DT','Cond'},[2,3]);
    dataT = make_between_table({mzMIsC,mzMIsA},dvn);
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0.25 1]);
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors = [mData.shades.m;mData.shades.c];
    % tcolors = repmat(tcolors,2,1)';
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    set(gca,'xlim',[0.35 xdata(end)+.65],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
    xticks = xdata; xticklabels = {'C3-D','C4-D','C3''-D','C3-T','C4-T','C3''-T'}; xticklabels = repmat(xticklabels,1,2);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.05 -0.011])
    put_axes_labels(gca,{[],[0 0 0]},{{'zMI'},[0 0 0]});

    save_pdf(hf,mData.pdf_folder,'zMI_bar_graph_all_cells.pdf',600);
break;
end
%% spatial and population vector correlation single animal
while 1
    an = 1; cn = 2;
    ff = makeFigureRowsCols(106,[1 0.5 6 1],'RowsCols',[2 2],...
        'spaceRowsCols',[0.09 -0.03],'rightUpShifts',[0.07 0.1],'widthHeightAdjustment',...
        [-100 -150]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 5 3.25 2]);
%     [resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC,resp_exc_inh] = get_responsive_fraction(Rs);
%     resp = get_cell_list(resp_valsC,[]);
    [CRc,aCRc,mRR,spCR] = find_population_correlations(Rs,mR,1,0);
    axes(ff.h_axes(1,1)); imagesc(mRR{an,cn}); colorbar; set(gca,'Ydir','Normal');
    axes(ff.h_axes(2,1)); imagesc(CRc{an,cn}); colorbar; set(gca,'Ydir','Normal');
    axes(ff.h_axes(1,2)); imagesc(spCR{an,cn}); colorbar; set(gca,'Ydir','Normal');
    set(ff.h_axes(2,2),'Visible','Off');
    set_colormap;
    set_obj(ff,{'FontWeight','Normal','FontSize',6,'LineWidth',0.25});
    ht = get_obj(ff,'title'); hyl = get_obj(ff,'ylabel'); changePosition(hyl(1,1),[4 0 0]);
%     set_obj(ht,'String',{'Pop. Activity','Pop. Activity','Pop. Activity';'Pop. Correlation','Pop.Correlation','Pop.Correlation'});
    set_obj(ht,{'FontSize',5,'FontWeight','Normal'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_correlations.pdf'),600);
    break;
end

%% Percentage of Responsive Cells
while 1
    [within,dvn,xlabels] = make_within_table({'FR','Cond'},[2,3]);
    dataT = make_between_table({resp_fractionCS*100,resp_fractionCSB*100},dvn);
    ra = RMA(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR_by_Cond','bonferroni'},[1 1 1]);
%     s = generate_shades(3); tcolors = s.g;
    tcolors = colors;%mData.colors(1:3);
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
    format_axes(gca);
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[0 25])
    xticks = [xdata(1:end)]; xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    any_mean = mean(100*resp_OR_fractionCS);    any_sem = std(100*resp_OR_fractionCS)/sqrt(5);
    pmchar=char(177); any_text = sprintf('%.0f%c%.0f%%',any_mean,pmchar,any_sem); text(0.75,20,any_text,'FontSize',6);
    any_mean = mean(100*resp_OR_fractionCSB);    any_sem = std(100*resp_OR_fractionCSB)/sqrt(5);
    pmchar=char(177); any_text = sprintf('%.0f%c%.0f%%',any_mean,pmchar,any_sem); text(4.75,20,any_text,'FontSize',6);
    changePosition(gca,[0.17 0.02 -0.1 -0.011])
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('percentage_PCs_responsive'),600);
    break;
end

%% Percentage of Responsive Cells Unique in Conditions
while 1
    percCells = [];
    percCells(:,1) = get_cell_list(resp_valsCS,[1 -2 -3],1)'; percCells(:,2) = get_cell_list(resp_valsCS,[-1 2 -3],1)'; percCells(:,3) = get_cell_list(resp_valsCS,[-1 -2 3],1)';
    within = make_within_table({'Cond'},3);
    dataT = make_between_table({percCells*100},{'C31','C4','C32'});
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
%     s = generate_shades(3); tcolors = s.g;
    tcolors = mData.colors;
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
    format_axes(gca);
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[0 25])
    xticks = [xdata(1:end)]; xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    any_mean = mean(100*resp_OR_fractionCS);    any_sem = std(100*resp_OR_fractionCS)/sqrt(5);
%     pmchar=char(177); any_text = sprintf('%.0f%c%.0f%%',any_mean,pmchar,any_sem); text(0.75,20,any_text,'FontSize',6);
    changePosition(gca,[0.17 0.02 -0.4 -0.011])
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('percentage_PCs_responsive_unique'),600);
    break;
end
%% Spatial correlation between adjacent trails
while 1
    [within,dvn,xlabels] = make_within_table({'Condition','TrialPairs'},[3,9]);
    var1 = arrayfun(@(x) mean(x{1}),out1.adj_SP_corr_diag); 
    var2 = arrayfun(@(x) mean(x{1}),out2.adj_SP_corr_diag);
    var3 = arrayfun(@(x) mean(x{1}),out3.adj_SP_corr_diag);
    between = make_between_table({var1,var2,var3},dvn);
    ra = RMA(between,within);
    %%
    hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_line_graph(mData,ra,0);
    tcolors = mData.colors;
    ii = 1; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 2; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linestyle','-.','linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 3; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    legs = {'C3','C4','C3'''};
    legs{end+1} = [1 0.2 0.25 0.25];
    putLegendH(gca,legs,tcolors)
    ylims = ylim;
    set(gca,'xlim',[0.5 xdata(end)+0.5],'ylim',[ylims(1) ylims(2)],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'T1-T2','T2-T3','T3-T4','T4-T5','T5-T6','T6-T7','T7-T8','T8-T9','T9-T10'};
    xticklabels = repmat(xticklabels,1,2);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(45);
    changePosition(gca,[0.05 0.12 0.05 -0.08]);
    put_axes_labels(gca,{'Trial-Pairs',[0 0 0]},{{'Spatial Correlation'},[0 -0.05 0]});
    save_pdf(hf,mData.pdf_folder,'spatial_correlation_trials',600);
    break;
end
%% Spatial correlation between adjacent trails (taking mean) 
while 1
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'TrialPairs','bonferroni'},[1 0.5 1]);
    hollowsep = 19;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.75 1],'color','w');
    hold on;
    tcolors = colors(1:9); tcolors = repmat(tcolors,[1 3]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
            'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
            'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    ylims = ylim;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[ylims(1) maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'T1-T2','T2-T3','T3-T4','T4-T5','T5-T6','T6-T7','T7-T8','T8-T9','T9-T10'};
    xticklabels = repmat(xticklabels,1,3);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(45);
%     changePosition(gca,[0.1 0.11 -0.2 -0.05]);
    changePosition(gca,[0.0 0.03 0.08 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatial','Correlation'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'spatial_correlation_trials_bar_graph',600);
    break;
end

%% Spatial correlation trails with mean of trials
while 1
    [within,dvn,xlabels] = make_within_table({'Condition','TrialPairs'},[3,10]);
    var1 = []; var2 = []; var3 = [];
    for ii = 1:length(out1.all_SP_corr_diag_with_mean)
        var1(ii,:) = (arrayfun(@(x) mean(x{1}),out1.all_SP_corr_diag_with_mean{ii}))'; 
        var2(ii,:) = (arrayfun(@(x) mean(x{1}),out2.all_SP_corr_diag_with_mean{ii}))'; 
        var3(ii,:) = (arrayfun(@(x) mean(x{1}),out3.all_SP_corr_diag_with_mean{ii}))'; 
    end
    between = make_between_table({var1,var2,var3},dvn);
    ra = repeatedMeasuresAnova(between,within);
    %%
    hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_line_graph(mData,ra,0);
    tcolors = mData.colors;
    ii = 1; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 2; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linestyle','-.','linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 3; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    legs = {'C3','C4','C3'''};
    legs{end+1} = [1 0.2 0.25 0.25];
    putLegendH(gca,legs,tcolors)
    ylims = ylim;
    set(gca,'xlim',[0.5 xdata(end)+0.5],'ylim',[ylims(1) ylims(2)],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'T1-T2','T2-T3','T3-T4','T4-T5','T5-T6','T6-T7','T7-T8','T8-T9','T9-T10'};
    xticklabels = repmat(xticklabels,1,2);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(45);
    changePosition(gca,[0.05 0.12 0.05 -0.08]);
    put_axes_labels(gca,{'Trial-Pairs',[0 0 0]},{{'Spatial Correlation'},[0 -0.05 0]});
    save_pdf(hf,mData.pdf_folder,'spatial_correlation_trials',600);
    break;
end
%% Spatial correlation trails with mean of trials (mean bar graph)
while 1
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0.5 1]);
    hollowsep = 19;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 3.25 1],'color','w');
    hold on;
    tcolors = colors(1:10); tcolors = repmat(tcolors,[1 3]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
            'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
            'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    ylims = ylim;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[ylims(1) maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'T1','T2','T3','T4-T5','T5-T6','T6-T7','T7-T8','T8-T9','T9-T10'};
    xticklabels = repmat(xticklabels,1,3);
    set(gca,'xtick',xticks,'xticklabels',xlabels);
    xtickangle(45);
%     changePosition(gca,[0.1 0.11 -0.2 -0.05]);
    changePosition(gca,[-0.02 0.03 0.1 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Correlation'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'spatial_correlation_trials_bar_graph',600);
    break;
end
%% Rate remap (Delta FR score) between adjacent trails
while 1
    [within,dvn,xlabels] = make_within_table({'Condition','TrialPairs'},[3,9]);
    var1 = arrayfun(@(x) nanmean(x{1}),out1.adj_RR_SP); 
    var2 = arrayfun(@(x) nanmean(x{1}),out2.adj_RR_SP);
    var3 = arrayfun(@(x) nanmean(x{1}),out3.adj_RR_SP);
    between = make_between_table({var1,var2,var3},dvn);
    ra = repeatedMeasuresAnova(between,within);
    %%
    hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_line_graph(mData,ra,0);
    tcolors = mData.colors;
    ii = 1; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 2; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linestyle','-.','linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 3; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    legs = {'C3','C4','C3'''};
    legs{end+1} = [1.5 0.3 1.02 0.25];
    putLegendH(gca,legs,tcolors)
    ylims = ylim;
    set(gca,'xlim',[0.5 xdata(end)+0.5],'ylim',[ylims(1)-0.01 ylims(2)],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'T1-T2','T2-T3','T3-T4','T4-T5','T5-T6','T6-T7','T7-T8','T8-T9','T9-T10'};
    xticklabels = repmat(xticklabels,1,2);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(45);
    changePosition(gca,[0.065 0.12 0.05 -0.08]);
    put_axes_labels(gca,{'Trial-Pairs',[0 0 0]},{{'\Delta FR Score'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'Delta_FR_score_trials',600);
    break;
end
%% Spatial correlation between conditions
while 1
    [within,dvn,xlabels] = make_within_table({'Cond'},2);
    var_CE = arrayfun(@(x) mean(x{1}),outRemap.adj_SP_corr_diag);
    dataT = make_between_table({var_CE},dvn);
    ra = repeatedMeasuresAnova(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0 1]);
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors ={colors{7};colors{8};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.02,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C3-C4','C4-C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     set(hbs(2),'facecolor','none','edgecolor',tcolors{2});
    xtickangle(45);
%     changePosition(gca,[0.2 0.03 -0.55 -0.05]);
    changePosition(gca,[0.1 0.03 -0.4 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatial Correlation'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'air_place_cell corr remap',600);
    break;
end

%% Delta FR score (RAte remap) between conditions
while 1
    [within,dvn,xlabels] = make_within_table({'Cond'},2);
    var_CE = arrayfun(@(x) nanmean(x{1}),outRemap.adj_RR_SP);
    dataT = make_between_table({var_CE},dvn);
    ra = repeatedMeasuresAnova(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0 1]);
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors ={colors{7};colors{8};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.02,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+0.1],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C3-C4','C4-C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     set(hbs(2),'facecolor','none','edgecolor',tcolors{2});
    xtickangle(45);
%     changePosition(gca,[0.2 0.03 -0.55 -0.05]);
    changePosition(gca,[0.1 0.03 -0.4 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'\Delta FR Score'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'air_place_cell rate remap (DFR)',600);
    break;
end

%% Popuation Vector correlation between conditions
while 1
    [within,dvn,xlabels] = make_within_table({'Cond'},2);
    var_CE = arrayfun(@(x) nanmean(x{1}),outRemap.adj_PV_corr_diag);
    dataT = make_between_table({var_CE},dvn);
    ra = repeatedMeasuresAnova(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0 1]);
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors ={colors{7};colors{8};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+0.1],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C3-C4','C4-C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     set(hbs(2),'facecolor','none','edgecolor',tcolors{2});
    xtickangle(45);
%     changePosition(gca,[0.2 0.03 -0.55 -0.05]);
    changePosition(gca,[0.2 0.03 -0.4 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Pop. Vec.','Correlation'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'air_place_cell PV corr',600);
    break;
end

%%
%% overlap RM bar graph
while 1
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    
    OI_E = get_overlap_index(resp_valsCS);
    for ii = 1:5
        C12(ii,1) = OI_E{ii}(1,2);
        C12(ii,2) = OI_E{ii}(2,3);
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},2);
    var_CE = C12;
    dataT = make_between_table({var_CE},dvn);
    ra = repeatedMeasuresAnova(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0 1]);
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors ={colors{7};colors{8};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+0.1],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C3-C4','C4-C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     set(hbs(2),'facecolor','none','edgecolor',tcolors{2});
    xtickangle(45);
%     changePosition(gca,[0.2 0.03 -0.55 -0.05]);
    changePosition(gca,[0.1 0.03 -0.4 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Overlap Index'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('air_overlap_stats_place_cells'),600);
    break;
end

%% Place Field Properties
props = {'Field Width (cm)','Field Center (cm)','Field FR (AU)'};
fileNames = {'PFW','PFC','PFFR'};
all_incr = [1 1 1];
% legs_v = 
for pri = 1%:3;
if 1
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            R = Rs{rr,cc};
            [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
            switch pri
                case 1
                    zMIsC{rr,cc} = PWs(respAnB{rr,cc})';
                case 2
                    zMIsC{rr,cc} = centers(respAnB{rr,cc})';
                case 3
                    zMIsC{rr,cc} = MFR(respAnB{rr,cc})';
            end
        end
    end

    tcolors = mData.colors;
    distD = zMIsC;
    [distDo,allVals] = getAveragesAndAllValues(distD);
    minBin = min(allVals);
    maxBin = max(allVals);
    incr = all_incr(pri); %maxBin =
    hf = get_figure(8,[5 7 1.25 1]);
    [ha,hb,~,bins] = plotAverageDistributions(distD,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
    format_axes(gca);
    if pri == 1
        legs = {'C3','C4','C3''',[50 3 70 5]}; putLegend(gca,legs,tcolors);
    end
    changePosition(gca,[0.1 0.13 -0.25 -0.13]);
    put_axes_labels(gca,{props{pri},[0 0 0]},{{'Neurons (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_%s',fileNames{pri}),600);
end

if 1
    mzMIsC = arrayfun(@(x) nanmean(x{1}),zMIsC);
    [within,dvn,xlabels] = make_within_table({'Cond'},[3]);
    dataT = make_between_table({mzMIsC},dvn);
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0.25 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.07 0.02 -0.4 -0.011])
    put_axes_labels(gca,{[],[0 0 0]},{props{pri},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells.pdf',fileNames{pri}),600);
end
end

%% Place Field Properties widths vs centers scatter plot and correlation
props = {'Field Width (cm)','Field Center (cm)','Field FR (AU)'};
fileNames = {'PFW','PFC','PFFR'};
all_incr = [1 1 1];

if 1
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            R = Rs{rr,cc};
            [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
            W{rr,cc} = PWs(respAnB{rr,cc})';
            C{rr,cc} = centers(respAnB{rr,cc})';
            CW_corr(rr,cc) = corr(C{rr,cc},W{rr,cc});
        end
    end
    tcolors = mData.colors;
    hf = get_figure(8,[5 7 1.25 1]);hold on;
    an = 3; 
    for cn = 1:3
        scatter(C{an,cn},W{an,cn},2,tcolors{cn});
        [x,p] = polyfit(C{an,cn},W{an,cn},1);
        y_hat = polyval(x,C{an,cn});
        plot(C{an,cn},y_hat,'color',tcolors{cn});
    end
    format_axes(gca);
    changePosition(gca,[0.1 0.13 -0.25 -0.13]);
    put_axes_labels(gca,{'Field Center (cm)',[0 0 0]},{'Field Width (cm)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('scatter_%s','PFW_PFC'),600);

    [within,dvn,xlabels] = make_within_table({'Cond'},[3]);
    dataT = make_between_table({CW_corr},dvn);
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0.25 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1)-0.1 maxY]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.23 0.02 -0.4 -0.011])
    put_axes_labels(gca,{[],[0 0 0]},{{'Corr. Field Width','Vs Field Center'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells.pdf','CW_corr'),600);
end

%% Place Field Properties on the belt
respAnB = get_props
respCells1 = get_cell_list(respAnB,[1 -2 -3]);     respCells2 = get_cell_list(respAnB,[-1 2 -3]);     respCells3 = get_cell_list(respAnB,[-1 -2 3]);
all_respCells = {respCells1,respCells2,respCells3};
props = {'Field Width (cm)',{'Spatially Tuned','Cells (%)'},'Field FR (AU)'};
fileNames = {'PFW','PFC','PFFR'}; yspacing = [10 5 5]; 
pri = 2;
if 1
    sRs = Rs(:,1:3);
    minBin = 0;     maxBin = 150;     incr = 50; % choosing more than 3 bins, give significant anova but not significant multcompare
    % three bins also make sense because LED in condition 2 comes ON at
    % 110cm
    bins = minBin:incr:maxBin;
    num_cells = []; PFW_belt = []; MFR_belt = [];
    for rr = 1:size(sRs,1)
        for cc = 1:size(sRs,2)
            R = sRs{rr,cc};
            [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
            resp = respAnB{rr,cc}; 
%             resp = all_respCells{cc}{rr};
            centers = centers(resp)'; PWs = PWs(resp)'; MFR = MFR(resp)';
            for ii = 1:(length(bins)-1)
                binS = bins(ii); binE = bins(ii+1);
                inds = find(centers >= binS & centers < binE);
                num_cells(rr,cc,ii) = 100*length(inds)/length(resp);
                PFW_belt(rr,cc,ii) = mean(PWs(inds));
                MFR_belt(rr,cc,ii) = mean(MFR(inds));
            end
        end
    end
    mzMIsC = [];
    for rr = 1:size(sRs,1)
        switch pri
            case 1
                mzMIsC(rr,:) = reshape(squeeze(PFW_belt(rr,:,:))',1,(length(bins)-1)*size(sRs,2));                
            case 2
                mzMIsC(rr,:) = reshape(squeeze(num_cells(rr,:,:))',1,(length(bins)-1)*size(sRs,2));
            case 3
                mzMIsC(rr,:) = reshape(squeeze(MFR_belt(rr,:,:))',1,(length(bins)-1)*size(sRs,2));
        end
    end
    [within,dvn,xlabels] = make_within_table({'Cond','Bins'},[size(sRs,2),(length(bins)-1)]);
    dataT = make_between_table({mzMIsC},dvn);
    ra = RMA(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_Bins','bonferroni'},[1 0.25 1]);
    hf = get_figure(5,[8 7 1.75 1]);
    s = generate_shades(length(bins)-1);
    tcolors = colors;%[s.m;s.c;s.y];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',yspacing(pri),'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    if pri == 2
        set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);
    else
        set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);
    end
    format_axes(gca);
    xticks = xdata; xticklabels = {'C3-B1','C3-B2','C3-B3','C4-B1','C4-B2','C4-B3','C3''-B1','C3''-B2','C3''-B3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    if pri == 2
        changePosition(gca,[0.11 0.02 -0.03 -0.011])
    else
        changePosition(gca,[0.01 0.02 -0.03 -0.011])
    end
    put_axes_labels(gca,{[],[0 0 0]},{props{pri},[0 -2 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells_belt.pdf',fileNames{pri}),600);
    ra.ranova
    ra.mauchly
end


%% percent of spatially tuned cells on the belt
if 1
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Bins','bonferroni'},[1 0.25 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    s = generate_shades(length(bins)-1);
    tcolors = colors;%[s.m;s.c;s.y];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',15,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);
    format_axes(gca);
    xticks = xdata; xticklabels = xlabels;
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0 0.02 -0.35 -0.05])
    put_axes_labels(gca,{[],[0 0 0]},{props{pri},[0 -2 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells_belt_pooled.pdf',fileNames{pri}),600);
    ra.ranova
    ra.mauchly
end

%% Place cell emergence disruption stability
for tC = 1:4
if 1
    txtT = {'Unique','New','Disrupted','Stable'};
    or_cells = get_cell_list(respAnB,[1;2;3]);
    percCells = []; 
    switch tC
        case 1 % unique place cells
            respCells1 = get_cell_list(respAnB,[1 -2 -3]);    respCells2 = get_cell_list(respAnB,[-1 2 -3]);    respCells3 = get_cell_list(respAnB,[-1 -2 3]);
        case 2 % new place cells
            respCells1 = get_cell_list(respAnB,[-1 2]);    respCells2 = get_cell_list(respAnB,[-2 3]);    respCells3 = get_cell_list(respAnB,[-1 3]);
        case 3 % disrupted place cells
            respCells1 = get_cell_list(respAnB,[1 -2]);    respCells2 = get_cell_list(respAnB,[2 -3]);    respCells3 = get_cell_list(respAnB,[1 -3]);
        case 4 % remained place cells
            respCells1 = get_cell_list(respAnB,[1 2]);    respCells2 = get_cell_list(respAnB,[2 3]);    respCells3 = get_cell_list(respAnB,[1 3]);
    end

    for ii = 1:length(respCells1)
        percCells(ii,1) = 100*sum(respCells1{ii})/length(or_cells{ii});%length(respCells1{ii});
        percCells(ii,2) = 100*sum(respCells2{ii})/length(or_cells{ii});%length(respCells2{ii});
        percCells(ii,3) = 100*sum(respCells3{ii})/length(or_cells{ii});%length(respCells3{ii});
    end

    [within,dvn,xlabels] = make_within_table({'Conds'},[3]);
    if tC > 1
       xlabels = {'C3-C4','C4-C3''','C3-C3'''};
    end
    dataT = make_between_table({percCells},dvn);
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0.25 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    tcolors = colors;%[s.m;s.c;s.y];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',15,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);
    format_axes(gca);
    xticks = xdata; xticklabels = xlabels;
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.1 0.02 -0.35 -0.05])
    put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    text(0.75,maxY,txtT{tC},'FontSize',6)
    save_pdf(hf,mData.pdf_folder,sprintf('%s_cells.pdf',txtT{tC}),600);
    ra.ranova
    ra.mauchly
end
end


%% Place cell emergence disruption stability - big combined test
if 1
    respCells1 = get_cell_list(respAnB,[1]); respCells2 = get_cell_list(respAnB,[2]);    respCells3 = get_cell_list(respAnB,[3]);
    respCells1o2 = get_cell_list(respAnB,[1;2]); respCells2o3 = get_cell_list(respAnB,[2;3]);    respCells1o3 = get_cell_list(respAnB,[1;3]);
    or_cells = get_cell_list(respAnB,[1;2;3])
    txtT = {'New','Disrupted'};
    percCells = [];
    % new place cells
    respCells12 = get_cell_list(respAnB,[-1 2]);    respCells23 = get_cell_list(respAnB,[-2 3]);    respCells13 = get_cell_list(respAnB,[-1 3]);
    for ii = 1:length(respCells1)
        percCells(ii,1) = sum(respCells12{ii})/length(or_cells{ii});
        percCells(ii,2) = sum(respCells23{ii})/length(or_cells{ii});
        percCells(ii,3) = sum(respCells13{ii})/length(or_cells{ii});
    end
    % disrupted place cells
    respCells12 = get_cell_list(respAnB,[1 -2]);    respCells23 = get_cell_list(respAnB,[2 -3]);    respCells13 = get_cell_list(respAnB,[1 -3]);
    for ii = 1:length(respCells1)
        percCells(ii,4) = sum(respCells12{ii})/length(or_cells{ii});%length(respCells1{ii});
        percCells(ii,5) = sum(respCells23{ii})/length(or_cells{ii});%length(respCells2{ii});
        percCells(ii,6) = sum(respCells13{ii})/length(or_cells{ii});%length(respCells3{ii});
    end
     % remained place cells
    respCells12 = get_cell_list(respAnB,[1 2]);    respCells23 = get_cell_list(respAnB,[2 3]);    respCells13 = get_cell_list(respAnB,[1 3]);
    for ii = 1:length(respCells1)
        percCells(ii,7) = sum(respCells12{ii})/length(or_cells{ii});%length(respCells1{ii});
        percCells(ii,8) = sum(respCells23{ii})/length(or_cells{ii});%length(respCells2{ii});
        percCells(ii,9) = sum(respCells13{ii})/length(or_cells{ii});%length(respCells3{ii});
    end
    

    [within,dvn,xlabels] = make_within_table({'Conds','Type'},[3 3]);
    dataT = make_between_table({percCells},dvn);
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
    hf = get_figure(5,[8 7 2.5 1]);
    tcolors = colors;%[s.m;s.c;s.y];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);
    format_axes(gca);
    xticks = xdata; xticklabels = xlabels;
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0 0.02 -0.15 -0.05])
    put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('dynamics_cells.pdf'),600);
    ra.ranova
    ra.mauchly
end

%% trial by trial comparison - place location gaussian fitting
[resp_fractionCS,respAnB,OIC,mean_OICS,resp_ORCS,resp_OR_fractionCS,resp_ANDCS,resp_AND_fractionCS] = get_responsive_fraction(Rs);
resp = get_cell_list(respAnB,[1;2;3]);
out_tt_PC = find_trial_by_trial_Gauss_Fit_Results(Rs,mR,resp);
[within,dvn,xlabels] = make_within_table({'Conds','DC'},[3,9]);
dataT1 = [];
for rr = 1:size(out_tt_PC.diff_centers_from_mean,1)
    thisD = [];
    for cc = 1:size(out_tt_PC.diff_centers_from_mean,2)
        thisDT = nanmean(out_tt_PC.diff_centers{rr,cc});
        thisD = [thisD thisDT];
    end
    dataT1(rr,:) = thisD;
end

dataT = make_between_table({dataT1},dvn);
ra = RMA(dataT,within,0.05);
[xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'DC','bonferroni'},[1 1 1]);
hf = get_figure(5,[8 7 9 1]);
% s = generate_shades(length(bins)-1);
tcolors = colors;%[s.m;s.c;s.y];
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
ylims = ylim;
set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]);
format_axes(gca);
xticks = xdata; xticklabels = xlabels;
set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
changePosition(gca,[0.05 0.02 -0.35 -0.05])
put_axes_labels(gca,{[],[0 0 0]},{'Jitter Centers (cm)',[0 0 0]});
save_pdf(hf,mData.pdf_folder,sprintf('jitter_centers.pdf'),600);
ra.ranova
ra.mauchly