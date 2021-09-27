function hf = plot_bar_graph(mData,ra,facs,gaps)

if ~exist('gaps','var')
    gaps = [1 1 1];
end

[xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,facs,gaps);
% hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 5 1.25],'color','w');
hold on;
[hbs,maxY] = plot_bars_sig_lines(mVar,semVar,combs,[h p],'colors',colors,'sigColor','k',...
        'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);