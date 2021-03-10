
function plotAverageSpeedConditions(b,markers1,markers2,fn)
%%
n = 0;
%%
% ei = evalin('base','ei10_C');
ei1 = evalin('base','ei10_A');
ei2 = evalin('base','ei10_C1');
ei3 = evalin('base','ei10_C');
ei = [ei3 ei1(1:5) ei2([3])];
% ei = [ei3 ei1(2:5) ei2([1 2 3])];

mData = evalin('base','mData');

for an = 1:length(ei)
    for cc = 1:4
        thisspeed = ei{an}.plane{1}.contexts(cc).rasters.airT.speed;
        meanSpeedTrials(:,cc) = nanmean(thisspeed,2);
        thisspeed = ei{an}.plane{1}.contexts(cc).rasters.airIT.speed;
        meanSpeedTrialsI(:,cc) = nanmean(thisspeed,2);
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
data = [moas moasi];
dataT = table(data(:,1),data(:,5),data(:,2),data(:,6),data(:,3),data(:,7),data(:,4),data(:,8));
dataT.Properties.VariableNames = {varNames{1} varNamesI{1} varNames{2} varNamesI{2} varNames{3} varNamesI{3} varNames{4} varNamesI{4}};
dataT
within = table([varNames';varNamesI']);
columnText = cell(size(within,1),1);columnText(1:2:end)= varNames';columnText(2:2:end)= varNamesI';
within = table([varNames';varNamesI'],columnText);
within = table([1 1 2 2 3 3 4 4]',[1 2 1 2 1 2 1 2]');
within.Properties.VariableNames = {'Condition','TI'};
within.TI = categorical(within.TI);
within.Condition = categorical(within.Condition);

ra = repeatedMeasuresAnova(dataT,within);

mVar = ra.est_marginal_means.Mean;
semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
inds86 = (combs(:,2)== 8 | combs(:,2)== 6);
inds87 = (combs(:,2)== 8 & combs(:,1)== 7);
inds65 = (combs(:,2)== 6 & combs(:,1)== 5);
inds = inds86 & ~inds87 & ~inds65;
h(inds) = 0;
xdata = [1 2 4 5 7 8 10 11]; maxY = 25;
colors = mData.colors;
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
hold on;
tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'maxY',maxY,'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.7,'sigLinesStartYFactor',0.1);
for ii = 2:2:length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
end
maxY = maxY + 6;
set(gca,'xlim',[0.25 11.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
xticks = [1.5 4.5 7.5 10.5]; xticklabels = {'C1','C2','C3','C4'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
changePosition(gca,[0.1 0.02 -0.03 -0.011]);
put_axes_labels(gca,{[],[0 0 0]},{'Avg. Speed (cm/sec)',[0 0 0]});
rectangle(gca,'Position',[0.75 maxY-4 1 2],'edgecolor','k','facecolor','k');
text(1.85,maxY-3,'Trials','FontSize',5);
rectangle(gca,'Position',[6 maxY-4 1 2],'edgecolor','k');
text(7.2,maxY-3,'Inter-Trials','FontSize',5);
save_pdf(hf,mData.pdf_folder,'AverageSpeedConditions_C',600);

