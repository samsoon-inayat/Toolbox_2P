function plotAverageSpeedConditions_R(b,markers1,markers2,fn)
%%
n = 0;
%%
ei = evalin('base','ei');

mData = evalin('base','mData');

conds = [3 4 5];

for an = 1:length(ei)
    for cci = 1:length(conds)
        cc = conds(cci);
        thisspeed = ei{an}.plane{1}.contexts(cc).rasters.airT.speed;
        meanSpeedTrials(:,cci) = nanmean(thisspeed,2);
        thisspeed = ei{an}.plane{1}.contexts(cc).rasters.airIT.speed;
        meanSpeedTrialsI(:,cci) = nanmean(thisspeed,2);
    end
    meanSpeedTrialsAnimalsT(an,:) = mean(meanSpeedTrials);
    meanSpeedTrialsAnimalsIT(an,:) = mean(meanSpeedTrialsI);
end

moas = meanSpeedTrialsAnimalsT;
moasi = meanSpeedTrialsAnimalsIT;
for ii = 1:size(moas,2)
    varNames{ii} = sprintf('Trials_Cond%d',ii);
end
for ii = 1:size(moasi,2)
    varNamesI{ii} = sprintf('InterTrials_Cond%d',ii);
end

data = NaN(size(moas)); data = repmat(data,1,2);
data(:,1:2:size(data,2)) = moas; data(:,2:2:size(data,2)) = moasi;

data1 = data;%(:,[1 5 2 6 3 7 4 8]);
[within,dvn,xlabels] = make_within_table({'Conds','TI'},[3 2]);
dataT = make_between_table({data1},dvn);
ra = RMA(dataT,within,{'hsd'});
ra.ranova
% raO = repeatedMeasuresAnova(dataT,within);

[xdata,mVar,semVar,combs,p,h,colorsi,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Conds_by_TI','hsd'},[1 1 1]);
colors = [mData.colors(1);mData.colors(1);mData.colors(2);mData.colors(2);mData.colors(3);mData.colors(3)];
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
hold on;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',colors,'sigColor','k',...
    'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.7,'sigLinesStartYFactor',0.1);

maxY = maxY + 8;
set(gca,'xlim',[0.25 8.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
xticks = xdata; xticklabels = {'T','I'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
make_bars_hollow(hbs(2:2:end))
xtickangle(45);
changePosition(gca,[0.11 0.05 -0.15 -0.1]);
put_axes_labels(gca,{[],[0 0 0]},{{'Avg. Speed','(cm/sec)'},[0 0 0]});
format_axes(gca);
legs = {'C3','C4','C5',[0.5 0.75 43 1]}; 
putLegendH(gca,legs,mData.colors(1:3),'sigR',{[],'anova',[],6});
save_pdf(hf,mData.pdf_folder,'AverageSpeedConditions_C',600);
