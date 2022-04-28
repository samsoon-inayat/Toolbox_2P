function dist_dur_Hi_RF_1


%% check the difference in zMI for Dis and Dur 
while 1
    %%
    mean_dzMI = [];
    FD_Prop = dzMI_FD.diff_T_D; FT_Prop = dzMI_FT.diff_T_D; 
%     FD_Prop = dzMI_FD.rs.diff_T_D; FT_Prop = dzMI_FT.rs.diff_T_D; 
%     FD_Prop = dzMI_FD.HaFD.diff_T_D; FT_Prop = dzMI_FT.HaFD.diff_T_D; 
%     FD_Prop = dzMI_FD.HiFD.diff_T_D; FT_Prop = dzMI_FT.HiFD.diff_T_D; 
%     FD_Prop = props{rfi,1}.N_Resp_Trials; FT_Prop = props{rfi,3}.N_Resp_Trials;
%     FD_Prop = props{1}.mean_FR; FT_Prop = props{3}.mean_FR;
    for rfi = 4
        TD = FD_Prop;
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dur_cells_T)];
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dis_cells_T)];

        TD = FT_Prop;
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dur_cells_I)];
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dis_cells_I)];
    end
    [within,dvn,xlabels] = make_within_table({'TI','CT','Cond'},[2,2,3]);
    dataT = make_between_table({mean_dzMI},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova
    print_for_manuscript(ra)
    %%
    break;
end



%% check the difference in zMI for Dis and Dur 
while 1
    %%
    mean_dzMI = [];
    for rfi = 4
        TD = props{rfi,1}.zMI;
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dur_cells_T)];
        TD = props{rfi,2}.zMI;
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dis_cells_T)];

         TD = props{rfi,3}.zMI;
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dur_cells_I)];
        TD = props{rfi,4}.zMI;
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dis_cells_I)];
    end
    [within,dvn,xlabels] = make_within_table({'TI','CT','Cond'},[2,2,3]);
    dataT = make_between_table({mean_dzMI},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova
    print_for_manuscript(ra)
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_CT','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
    xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = repmat(mData.colors(1:2),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.75,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dur','Dis'};
    make_bars_hollow(hbs(3:4))
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    changePosition(gca,[0.04 0.01 -0.4 0]); put_axes_labels(gca,{[],[0 0 0]},{'zMI (A.U.)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_all.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
    xdata = make_xdata([2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = repmat(mData.colors(1:2),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dur','Dis'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    changePosition(gca,[0.04 0.01 -0.55 0]); put_axes_labels(gca,{[],[0 0 0]},{'zMI (A.U.)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_all.pdf'),600);
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
    xdata = make_xdata([2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = repmat(mData.colors(1:2),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'T','I'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    changePosition(gca,[0.04 0.01 -0.55 0]); put_axes_labels(gca,{[],[0 0 0]},{'zMI (A.U.)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_all.pdf'),600);
    %%
    break;
end

%%
