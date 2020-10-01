function figure1_Distributions_all_conditions

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
protocol = '15';
selAnimals = eval(sprintf('mData.selAnimals%s_d1',protocol));
ra1 = figure1_Distributions_all_conditions_15_1(selAnimals);
selAnimals = eval(sprintf('mData.selAnimals%s_d2',protocol));
ra2 = figure1_Distributions_all_conditions_15_1(selAnimals);

data = [ra1.rm.BetweenDesign;ra2.rm.BetweenDesign];
data = [table([ones(5,1);2*ones(5,1)]) data];
data.Properties.VariableNames{1} = 'Day';
within = ra1.rm.WithinDesign; varNames = data.Properties.VariableNames;
ra = repeatedMeasuresAnova(data{:,:},varNames,within,2);

rm = ra.rm;
mcTI = find_sig_mctbl(multcompare(rm,'Day','By','Condition','ComparisonType','bonferroni'),6);
mcDays = find_sig_mctbl(multcompare(rm,'Condition','By','Day','ComparisonType','bonferroni'),6);

dataR = [data{1:5,2:end} data{6:10,2:end}];
[mVar semVar] = findMeanAndStandardError(dataR);
[combs,h,p] = populate_multcomp_h_p(dataR,within,[],[]);

xdata = [1:1.5:(10*size(dataR,2))]; xdata = xdata(1:size(dataR,2)); 
ind = 1;
for ii = 1:2
    for jj = 1:3
        tcolors{ind} = colors{ii};
        ind = ind + 1;
    end
end
hf = figure(15);clf;set(gcf,'Units','Inches');set(gcf,'Position',[3 3 6 1.5],'color','w');
hold on;
[hbs,maxYr] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
set(gca,'xlim',[0.25 max(xdata)+.75],'ylim',[0 maxYr],'FontSize',7,'FontWeight','Bold','TickDir','out');
xticks = xdata; 
%     xticklabels = repmat(xticklabels,length(all_rts),length(all_conds));
set(gca,'xtick',xticks,'xticklabels',xticklabels);
xtickangle(30);
changePosition(gca,[0 0.02 0.03 -0.011])
put_axes_labels(gca,{[],[0 0 0]},{{'Mutual Information','(z-score)'},[0 0 0]});
save_pdf(hf,mData.pdf_folder,sprintf('Mean zMI all Conditions'),600);
