%% TI_CT_Cond
while 1
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_CT_Cond','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3 3 3 3],[1 1.5]);
    hf = get_figure(5,[8 7 3.5 1]);
    tcolors = repmat(mData.colors(1:9),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
%     changePosition(gca,[-0.04 0.01 0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
%     put_axes_labels(gca,{[],[0 0 0]},{'zMI Difference',[0 0 0]});
%     save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_all.pdf'),600);
break;
end
%% TI_CT
while 1
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_CT','hsd'},[1.5 1 1]);
    xdata = make_xdata([3 3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = repmat(mData.colors(1:3),1,2);
    MmVar = max(mVar);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',MmVar/1.25,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'tDur','tDis','Ind'};
    make_bars_hollow(hbs(4:end))
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.04 0.01 -0.1 0]); 
% put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
%     save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_all.pdf'),600);
break;
end
%% CT
while 1
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1.5 1 1]);
    xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1.25 1]);
    tcolors = repmat(mData.colors(1:9),1,2);
    MmVar = max(mVar);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',MmVar/0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'tDur','tDis','Ind'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.06 0.01 -0.35 0]); 
% put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
%     save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_all.pdf'),600);
break;
end
%% CT cond
while 1
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT_by_Cond','hsd'},[1.5 1 1]);
    xdata = make_xdata([3 3 3],[1 1.5]);
    hf = get_figure(5,[8 7 1.25 1]);
    tcolors = repmat(mData.colors(1:9),1,2);
    MmVar = max(mVar);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',MmVar/3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.06 0.01 -0.1 0]); 
%     put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
%     save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_all.pdf'),600);
break;
end
%% Cond CT
while 1
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_CT','hsd'},[1.5 1 1]);
    xdata = make_xdata([3 3 3],[1 2]);
    hf = get_figure(5,[8 7 1.75 1]);
    tcolors = repmat(mData.colors(1:3),1,3);
    MmVar = max(mVar);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',MmVar/3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dur','Dis','Vag'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.02 0.01 0.0 0]); 
%     put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
%     save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_all.pdf'),600);
break;
end
 %% TI
 while 1
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI','hsd'},[1.5 1 1]);
    xdata = make_xdata([2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(7:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.75,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'T','I','Mix'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.4 0]); 
    break;
 end
 %% TI cond
while 1
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_Cond','hsd'},[1.5 1 1]);
    xdata = make_xdata([3 3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = repmat(mData.colors(1:3),1,2);
    MmVar = max(mVar);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',MmVar/3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.06 0.01 -0.05 0]); 
%     put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
%     save_pdf(hf,mData.pdf_folder,sprintf('var_TI_by_Cond.pdf'),600);
break;
end

%% print values for stats in manuscript
print_for_manuscript(ra)

%% percentage format
figure(hf);
set(gca,'ytick',[10 20 30]);
ylabel('Cells (%)');

%% zMI format
figure(hf);
set(gca,'ytick',[-2 -1 0 1 2]);
ylim([-2.5 4]);
ylabel('dzMI (A.U.)');

%% save file
save_pdf(hf,mData.pdf_folder,sprintf('omnifile.pdf'),600);