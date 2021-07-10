
function plotAverageSpeedConditions(b,markers1,markers2,fn)
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
dataC = data;
dataC = [ones(size(dataC,1),1) dataC];
dataT = array2table(dataC(:,2:end));
dataT.Properties.VariableNames = {varNames{1} varNamesI{1} varNames{2} varNamesI{2} varNames{3} varNamesI{3}};
dataT
within = table([1 1 2 2 3 3]',[1 2 1 2 1 2]');
within.Properties.VariableNames = {'Condition','TI'};
within.TI = categorical(within.TI);
within.Condition = categorical(within.Condition);

ra = repeatedMeasuresAnova(dataT,within,0.05);
%%
mVar = ra.est_marginal_means.Mean;
semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
xdata = [1 2 4 5 7 8]; maxY = 25;
colors = mData.colors;
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.35 1],'color','w');
hold on;
tcolors = {colors{3};colors{3};colors{4};colors{4};colors{5};colors{5};colors{6};colors{6}};
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'maxY',maxY,'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.7,'sigLinesStartYFactor',0.1);
for ii = 2:2:length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
end
maxY = maxY + 6;
set(gca,'xlim',[0.25 8.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
xticks = [1.5 4.5 7.5]; xticklabels = {'C3','C4','C3'''};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
xtickangle(30);
changePosition(gca,[0.1 0.05 -0.07 -0.1]);
put_axes_labels(gca,{[],[0 0 0]},{{'Avg. Speed (cm/sec)'},[0 0 0]});
rectangle(gca,'Position',[0.75 maxY-4 1 2],'edgecolor','k','facecolor','k');
text(1.85,maxY-3,'Trials','FontSize',5);
rectangle(gca,'Position',[4 maxY-4 1 2],'edgecolor','k');
text(5.2,maxY-3,'Inter-Trials','FontSize',5);
save_pdf(hf,mData.pdf_folder,'AverageSpeedConditions_C',600);

