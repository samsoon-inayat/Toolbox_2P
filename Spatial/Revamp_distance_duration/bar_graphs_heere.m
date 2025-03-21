%% Resp. two graphs, all and cond
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[10 5 2.75 1],'RowsCols',[1 2],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.08 0.36],'widthHeightAdjustment',[10 -510]);
switch varT
    case 1 % responsive cells 
        MY = 70; ysp = 4; mY = 0; titletxt = 'Responsiveness'; ylabeltxt = {'Cells (%)'}; % for all cells (vals) MY = 80
    case 2
        MY = 90; ysp = 4; mY = 0; titletxt = 'Response Fidelity'; ylabeltxt = {'% of Trials'};% for all cells (vals) MY = 70
    case 3
        MY = 0.1; ysp = 0.01; mY = -0.15; titletxt = 'Mutual Information'; ylabeltxt = {'Mutual Information','Z-Score'};
        MY = 3.5; ysp = 0.3; mY = 0; titletxt = 'Mutual Information'; ylabeltxt = {'(z-score)'};
    case 10
        MY = 200; ysp = 1; mY = -0.15; titletxt = 'Tuning Width'; ylabeltxt = {'cm'};
end
stp = 0.25*magfac; widths = ([1.85 0.65 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = repmat(mData.colors(1:2),1,6);
axes(ff.h_axes(1,1));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_Cond_by_DT','hsd'},[1.5 1 1]);
    xdata = make_xdata([2 2 2 2 2 2],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
    h(h==1) = 0;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,[],[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'Ti','Di'};set(gca,'xtick',xticks,'xticklabels',xticklabels);% xtickangle(30);
make_bars_hollow(hbs(7:end));
[~,hyl] = put_axes_labels(gca,{'',[]},{ylabeltxt,[]}); set(hyl,'FontWeight','bold');
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'C3','C4','C5','C3','C4','C5'},{[0 0.03]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,6,{'AOn','AOff'},{[-0.1 -0.012]});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.01 0 0]);set(ht,'FontWeight','Normal');
format_axes(gca);

axes(ff.h_axes(1,2));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_DT','hsd'},[1.5 1 1]);
    xdata = make_xdata([4],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
%     h(h==1) = 0;
tcolors = repmat(mData.colors(1:2),1,6);
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'Ti','Di'};set(gca,'xtick',xticks,'xticklabels',xticklabels); %xtickangle(30);
make_bars_hollow(hbs(3:end));
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Trials','InterTrials'},{[-0.1 -0.012]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Pooled'},{[0 0]});
format_axes(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% RF two graphs, all and cond
ff = makeFigureRowsCols(107,[10 5 2.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.36],'widthHeightAdjustment',[10 -510]);
switch varT
    case 1 % responsive cells 
        MY = 70; ysp = 4; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'% of cells'}; % for all cells (vals) MY = 80
    case 2
        MY = 80; ysp = 4; mY = 0; titletxt = 'Response Fidelity'; ylabeltxt = {'% of Trials'};% for all cells (vals) MY = 70
    case 3
        MY = 0.1; ysp = 0.01; mY = -0.15; titletxt = 'Mutual Information'; ylabeltxt = {'Mutual Information','Z-Score'};
        MY = 3.5; ysp = 0.3; mY = 0; titletxt = 'Mutual Information'; ylabeltxt = {'(z-score)'};
    case 9
        MY = 20; ysp = 1; mY = -0.15; titletxt = 'Tuning Width'; ylabeltxt = {'cm'};
end
stp = 0.225*magfac; widths = ([0.4 0.55 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = repmat(mData.colors(1:2),1,6);
axes(ff.h_axes(1,1));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'DT','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([2],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
%     h(h==1) = 0;
tcolors = repmat(mData.colors(1:2),1,6);
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'Ti','Di'};set(gca,'xtick',xticks,'xticklabels',xticklabels);% xtickangle(30);
make_bars_hollow(hbs(3:end));
[~,hyl] = put_axes_labels(gca,{'',[]},{ylabeltxt,[]}); set(hyl,'FontWeight','bold');
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Trials','InterTrials'},{[-0.1 -0.012]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0 0]});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[-0.05 -0.01 0.5 0]);set(ht,'FontWeight','Bold');
format_axes(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% zMI two graphs, all and cond
ff = makeFigureRowsCols(107,[10 5 2.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.36],'widthHeightAdjustment',[10 -510]);
switch varT
    case 1 % responsive cells 
        MY = 70; ysp = 4; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'% of cells'}; % for all cells (vals) MY = 80
    case 2
        MY = 90; ysp = 4; mY = 0; titletxt = 'Response Fidelity'; ylabeltxt = {'% of Trials'};% for all cells (vals) MY = 70
    case 3
        MY = 0.1; ysp = 0.01; mY = -0.15; titletxt = 'Mutual Information'; ylabeltxt = {'Mutual Information','Z-Score'};
        MY = 2.5; ysp = 0.3; mY = 0; titletxt = 'Mutual Information'; ylabeltxt = {'z-score'};
    case 9
        MY = 20; ysp = 1; mY = -0.15; titletxt = 'Tuning Width'; ylabeltxt = {'cm'};
end
stp = 0.225*magfac; widths = ([0.55 0.55 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = repmat(mData.colors(1:2),1,6);
axes(ff.h_axes(1,1));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_DT','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([4],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
%     h(h==1) = 0;
tcolors = repmat(mData.colors(1:2),1,6);
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'Ti','Di'};set(gca,'xtick',xticks,'xticklabels',xticklabels); %xtickangle(30);
[~,hyl] = put_axes_labels(gca,{'',[]},{ylabeltxt,[]}); set(hyl,'FontWeight','bold');
make_bars_hollow(hbs(3:end));
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Trials','InterTrials'},{[-0.1 -0.012]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Pooled'},{[0 0]});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[-0.05 -0.01 0.2 0]);set(ht,'FontWeight','Bold');
format_axes(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% RS two graphs, all and cond
ff = makeFigureRowsCols(107,[10 5 2.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.36],'widthHeightAdjustment',[10 -510]);
switch varT
    case 1 % responsive cells 
        MY = 70; ysp = 4; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'% of cells'}; % for all cells (vals) MY = 80
    case 2
        MY = 90; ysp = 4; mY = 0; titletxt = 'Response Fidelity'; ylabeltxt = {'% of Trials'};% for all cells (vals) MY = 70
    case 3
        MY = 0.1; ysp = 0.01; mY = -0.15; titletxt = 'Mutual Information'; ylabeltxt = {'Mutual Information','Z-Score'};
        MY = 3.5; ysp = 0.3; mY = 0; titletxt = 'Mutual Information'; ylabeltxt = {'(z-score)'};
    case 4
        MY = 1.3; ysp = 0.14; mY = 0; titletxt = 'Goodness of Fit'; ylabeltxt = {'r-squared'};
end
stp = 0.27*magfac; widths = ([0.55 0.55 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = repmat(mData.colors(1:2),1,6);
axes(ff.h_axes(1,1));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_DT','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([4],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
%     h(h==1) = 0;
tcolors = repmat(mData.colors(1:2),1,6);
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'Ti','Di'};set(gca,'xtick',xticks,'xticklabels',xticklabels); %xtickangle(30);
[~,hyl] = put_axes_labels(gca,{'',[]},{ylabeltxt,[]}); set(hyl,'FontWeight','bold');
make_bars_hollow(hbs(3:end));
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Trials','InterTrials'},{[-0.1 -0.012]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Pooled'},{[0 0]});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[-0.05 -0.01 0.2 0]);set(ht,'FontWeight','Bold');
format_axes(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% zMI two graphs, all and cond --- dzMI
ff = makeFigureRowsCols(107,[10 5 2.25 1.25],'RowsCols',[1 2],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.36],'widthHeightAdjustment',[10 -510]);
switch varT
    case 1 % responsive cells 
        MY = 70; ysp = 4; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'% of cells'}; % for all cells (vals) MY = 80
    case 2
        MY = 90; ysp = 4; mY = 0; titletxt = 'Response Fidelity'; ylabeltxt = {'% of Trials'};% for all cells (vals) MY = 70
    case 3
        MY = 0.1; ysp = 0.01; mY = -0.15; titletxt = 'Mutual Information'; ylabeltxt = {'Mutual Information','Z-Score'};
        MY = 2.5; ysp = 0.3; mY = 0; titletxt = 'Mutual Information'; ylabeltxt = {'(z-score)'};
    case 9
        MY = 20; ysp = 1; mY = -0.15; titletxt = 'Tuning Width'; ylabeltxt = {'cm'};
end
stp = 0.225*magfac; widths = ([0.55 0.55 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});

axes(ff.h_axes(1,1));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'DT','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([2],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
%     h(h==1) = 0;
tcolors = repmat(mData.colors(1:2),1,6);
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'T','D'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30);
[~,hyl] = put_axes_labels(gca,{'',[]},{ylabeltxt,[]}); set(hyl,'FontWeight','bold');
make_bars_hollow(hbs(3:end));
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Trials','InterTrials'},{[-0.1 -0.012]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0 0]});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.051 0.2 0]);set(ht,'FontWeight','Bold');
format_axes(gca);

axes(ff.h_axes(1,2));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI','hsd'},[1.5 1 1]);
    xdata = make_xdata([2],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
%     h(h==1) = 0;
tcolors = repmat(mData.colors(3:4),1,6);
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'TR','IT'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30);
% [~,hyl] = put_axes_labels(gca,{'',[]},{ylabeltxt,[]}); set(hyl,'FontWeight','bold');
make_bars_hollow(hbs(3:end));
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Trials','InterTrials'},{[-0.1 -0.012]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0 0]});
format_axes(gca);

save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);


%% TI_CT for responsivity of Dur, Dis, and Ind cells
ff = makeFigureRowsCols(107,[10 5 1.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.36],'widthHeightAdjustment',[10 -510]);
MY = 100; mY = 0;
stp = 0.3*magfac; widths = ([1.1 0.55 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = repmat(mData.colors(4:6),1,6);
axes(ff.h_axes(1,1));

% [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_CT','hsd'},[1.5 1 1]);
[xdata,mVar,semVar,combs,p,h] = get_vals_RMA(mData,ra,{'TI:CT','hsd'},[1 2]);
% xdata = make_xdata([3 3],[1 1.5]);
%     hf = get_figure(5,[8 7 1.5 1]);
% tcolors = repmat(mData.colors(1:3),1,2);
MmVar = max(mVar);
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',MmVar/1.95,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
maxY = maxY + 0;
ylims = ylim;
format_axes(gca);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);
xticks = xdata; xticklabels = {'TE','DE','IE'};
make_bars_hollow(hbs(4:end))
set(gca,'xtick',xticks,'xticklabels',xticklabels); %xtickangle(45)
%     changePosition(gca,[0.04 0.01 -0.1 0]); 
put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
ht = set_axes_top_text_no_line(gcf,gca,'Responsiveness',[-0.05 -0.01 0.2 0]);set(ht,'FontWeight','Bold');
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'AOn','AOff'});
format_axes(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% CT:Cond for responsivity of Dur, Dis, and Ind cells
ff = makeFigureRowsCols(107,[10 5 1.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.36],'widthHeightAdjustment',[10 -510]);
MY = 100; mY = 0;
stp = 0.245*magfac; widths = ([1.25 0.55 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes(ff.h_axes(1,1));

% [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT_by_Cond','bonferroni'},[1.5 1 1]);
[xdata,mVar,semVar,combs,p,h] = get_vals_RMA(mData,ra,{'CT:Cond','hsd'},[1 2]);

tcolors = repmat(mData.colors(7:9),3,1);
MmVar = max(mVar);
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',MmVar/1.95,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
maxY = maxY + 0;
ylims = ylim;
format_axes(gca);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);
xticks = xdata; xticklabels = {'C3','C4','C5'};
make_bars_hollow(hbs(4:end))
set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[]); %xtickangle(45)
%     changePosition(gca,[0.04 0.01 -0.1 0]); 
% put_axes_labels(gca,{[],[0 0 0]},{'% of Cells',[0 0 0]});
% ht = set_axes_top_text_no_line(gcf,gca,'Responsiveness',[-0.05 -0.01 0.2 0]);set(ht,'FontWeight','Bold');
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'TE','DE','IE'});
format_axes(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% TI_CT for szMI of Dur, Dis, and Ind cells
ff = makeFigureRowsCols(107,[10 5 1.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.36],'widthHeightAdjustment',[10 -510]);
MY = 7; mY = -3; ys = 1;
stp = 0.245*magfac; widths = ([1.25 0.55 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = repmat(mData.colors(1:2),1,6);
axes(ff.h_axes(1,1));

% [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_CT','bonferroni'},[1.5 1 1]);
% xdata = make_xdata([3 3],[1 1.5]);
[xdata,mVar,semVar,combs,p,h] = get_vals_RMA(mData,ra,{'TI:CT','hsd'},[1 2]);

%     hf = get_figure(5,[8 7 1.5 1]);
tcolors = repmat(mData.colors(4:6),1,2);
MmVar = max(mVar);
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ys,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);

format_axes(gca);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);
xticks = xdata; xticklabels = {'TE','DE','IE'};
make_bars_hollow(hbs(4:end))
set(gca,'xtick',xticks,'xticklabels',xticklabels); %xtickangle(45)
%     changePosition(gca,[0.04 0.01 -0.1 0]); 
put_axes_labels(gca,{[],[0 0 0]},{'{\Delta} z-score',[0 0 0]});
ht = set_axes_top_text_no_line(gcf,gca,'Mutual Information',[-0.05 -0.01 0.2 0]);set(ht,'FontWeight','Bold');
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Air','No-Air'});
format_axes(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% CT zMI
ff = makeFigureRowsCols(107,[10 5 1.55 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.36],'widthHeightAdjustment',[10 -510]);
  MY = 2; mY = -2; ys = 1;
  stp = 0.27*magfac; widths = ([0.75 0.55 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
  adjust_axes(ff,[mY MY],stp,widths,gap,{''});
  tcolors = repmat(mData.colors(1:2),1,6);
  axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([3],[1 1.5]);
%     hf = get_figure(5,[8 7 1.25 1]);
    tcolors = repmat(mData.colors(1:9),2,1);
    MmVar = max(mVar);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ys,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);
    xticks = xdata; xticklabels = {'TE','DE','IE'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);% xtickangle(45)
put_axes_labels(gca,{[],[0 0 0]},{'z-score',[0 0 0]});
ht = set_axes_top_text_no_line(gcf,gca,'Mutual Information',[-0.05 -0.01 0.2 0]);set(ht,'FontWeight','Bold');
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Air','No-Air'});
format_axes(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% CT
ff = makeFigureRowsCols(107,[10 5 1.55 1],'RowsCols',[1 2],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.36],'widthHeightAdjustment',[10 -510]);
  MY = 0.25; mY = -0.35; ys = 2;
  stp = 0.27*magfac; widths = ([0.75 0.55 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
  adjust_axes(ff,[mY MY],stp,widths,gap,{''});
  tcolors = repmat(mData.colors(1:2),1,6);
  axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','bonferroni'},[1.5 1 1]);
%     xdata = make_xdata([3],[1 1.5]);
%     hf = get_figure(5,[8 7 1.25 1]);
    tcolors = repmat(mData.colors(4:6),2,1);
    MmVar = max(mVar);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);
    xticks = xdata; xticklabels = {'TE','DE','IE'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);% xtickangle(45)
put_axes_labels(gca,{[],[0 0 0]},{'{\Delta} r-squared',[0 0 0]});
ht = set_axes_top_text_no_line(gcf,gca,'Goodness-of-Fit',[-0.05 -0.01 0.2 0]);set(ht,'FontWeight','Bold');
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Air','No-Air'});
format_axes(gca);
% save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

  axes(ff.h_axes(1,2));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([2],[1 1.5]);
%     hf = get_figure(5,[8 7 1.25 1]);
    tcolors = repmat(mData.dcolors(7:9),2,1);
    MmVar = max(mVar);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Air','No Air'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30)

% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Air','No-Air'});
format_axes(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% CT cond
while 1
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([3],[1 1.5]);
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
   ff = makeFigureRowsCols(107,[10 5 1.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.36],'widthHeightAdjustment',[10 -510]);
  MY = 25; mY = -1; ys = 2;
  stp = 0.245*magfac; widths = ([0.55 0.55 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
  adjust_axes(ff,[mY MY],stp,widths,gap,{''});
  tcolors = repmat(mData.colors(1:2),1,6);
  axes(ff.h_axes(1,1));
%     [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI','bonferroni'},[1.5 1 1]);
%     xdata = make_xdata([2],[1 1.5]);
    [xdata,mVar,semVar,combs,p,h] = get_vals_RMA(mData,ra,{'TI','hsd'},[1 2]);

%     hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(7:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ys,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);
    xticks = xdata; xticklabels = {'AOn','AOff'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    put_axes_labels(gca,{[],[0 0 0]},{'{\Delta} % of trials',[0 0 0]});
ht = set_axes_top_text_no_line(gcf,gca,'Response Fidelity',[-0.05 -0.01 0.2 0]);set(ht,'FontWeight','Bold');
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'});
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
 %% TI cond
while 1
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_Cond','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = repmat(mData.colors(1:2),1,3);
    MmVar = max(mVar);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',MmVar/3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'D','T'};
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
% set(gca,'ytick',[-1 0 1 ]);
% ylim([-2 4]);
ylabel('dzMI (A.U.)');

%% save file
save_pdf(hf,mData.pdf_folder,sprintf('omnifile.pdf'),600);


%% CT zMI
ff = makeFigureRowsCols(107,[10 5 1.55 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.36],'widthHeightAdjustment',[10 -510]);
  MY = 150; mY = 0; ys = 0.1;
  stp = 0.27*magfac; widths = ([0.75 0.55 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
  adjust_axes(ff,[mY MY],stp,widths,gap,{''});
  tcolors = repmat(mData.colors(1:2),1,6);
  axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([2],[1 1.5]);
%     hf = get_figure(5,[8 7 1.25 1]);
    tcolors = repmat(mData.colors(1:9),2,1);
    MmVar = max(mVar);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ys,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Time','Dist'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);% xtickangle(45)
put_axes_labels(gca,{[],[0 0 0]},{'z-score',[0 0 0]});
ht = set_axes_top_text_no_line(gcf,gca,'Mutual Information',[-0.05 -0.01 0.2 0]);set(ht,'FontWeight','Bold');
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Air','No-Air'});
format_axes(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% CT zMI
ff = makeFigureRowsCols(107,[10 5 1.55 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.36],'widthHeightAdjustment',[10 -510]);
  MY = 5; mY = 0; ys = 1;
  stp = 0.27*magfac; widths = ([0.75 0.55 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
  adjust_axes(ff,[mY MY],stp,widths,gap,{''});
  tcolors = repmat(mData.colors(1:2),1,6);
  axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([2],[1 1.5]);
%     hf = get_figure(5,[8 7 1.25 1]);
    tcolors = repmat(mData.colors(1:9),2,1);
    MmVar = max(mVar);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ys,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);
    xticks = xdata; xticklabels = {'TE','DE','IE'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);% xtickangle(45)
put_axes_labels(gca,{[],[0 0 0]},{'z-score',[0 0 0]});
ht = set_axes_top_text_no_line(gcf,gca,'Mutual Information',[-0.05 -0.01 0.2 0]);set(ht,'FontWeight','Bold');
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Air','No-Air'});
format_axes(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);