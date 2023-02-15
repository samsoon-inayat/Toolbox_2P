   %% TI_CT for responsivity of Dur, Dis, and Ind cells
switch varT
  case 1
    MY = 1; mY = 0; ysp = 0.05; titletext = 'Responsiveness'; y_label = '% of cells';
  case 2
    MY = 100; mY = 0; ysp = 10; titletext = 'Response Fidelity'; y_label = '% of trials';
  case 3
    MY = 3; mY = 0; ysp = 0.4; titletext = 'Mutual Information'; y_label = 'z-score';
  case 4
    MY = 1.1; mY = 0; ysp = 0.12; titletext = 'Goodness-of-Fit'; y_label = 'r-squared';
end
ff = makeFigureRowsCols(107,[10 5 2.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.36],'widthHeightAdjustment',[10 -510]);
stp = 0.3*magfac; widths = ([1.95 0.55 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = repmat(mData.colors(4:5),1,6);

axes(ff.h_axes(1,1));

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TD_by_PH','bonferroni'},[1.5 1 1]);
xdata = make_xdata([2 2 2],[1 1.5]);
%     hf = get_figure(5,[8 7 1.5 1]);
% tcolors = repmat(mData.colors(1:3),1,2);
MmVar = max(mVar);
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.015);
maxY = maxY + 0;
ylims = ylim;
format_axes(gca);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);
xticks = xdata; xticklabels = {'Time','Dist'};
make_bars_hollow(hbs(3:end))
set(gca,'xtick',xticks,'xticklabels',xticklabels); %xtickangle(45)
%     changePosition(gca,[0.04 0.01 -0.1 0]); 
put_axes_labels(gca,{[],[0 0 0]},{y_label,[0 0 0]});
ht = set_axes_top_text_no_line(gcf,gca,titletext,[-0.05 -0.01 0.2 0]);set(ht,'FontWeight','Bold');
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Air','No-Air'});
format_axes(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
