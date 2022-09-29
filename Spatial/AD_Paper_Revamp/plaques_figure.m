function plaques_figure

Plaques_APP = [54.07282243
35.39646519
45.35599232
37.59057383
92.69148972
84.53510672
];

Plaques_Ctrl = [0
0
0
0
0
0
];

mData = evalin('base','mData');
n = 0;

mVar = [mean(Plaques_Ctrl) mean(Plaques_APP)];
semVar = [std(Plaques_Ctrl)/sqrt(length(Plaques_Ctrl)) std(Plaques_APP)/sqrt(length(Plaques_APP))];

[h,p,ci,stats] = ttest2(Plaques_Ctrl,Plaques_APP);
combs = [1 2];

ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
set(gcf,'color','w');    set(gcf,'Position',[10 3 1.9 1]);
MY = 80; ysp = 8; mY = 0; % responsive cells
stp = 0.3; widths = [0.55 1.3 1.3 1.3 1.3 0.5 0.5 0.5]-0.15; gap = 0.16;
adjust_axes(ff,[mY MY],stp,widths,gap,{'R-squared'});
colors = {'k','r'};
tcolors = {colors{1};colors{2}};

% [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_Cond','hsd'},[1.5 1 1]);
xdata = make_xdata([2],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;

[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
    xticklabels = {'Control','APP'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30);
    make_bars_hollow(hbs(5:end))
    ylabel('count/mm^2');
    set_axes_top_text_no_line(ff.hf,gca,'Plaque Density',[-0.1 0 0.1 0]);
    save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);