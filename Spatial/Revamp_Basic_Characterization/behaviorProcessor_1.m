function behaviorProcessor_1
temp = evalin('base','training_data');
numberOfTrials = findNumberOfTrials(temp);
[rr,cc] = find(numberOfTrials < 20);

ei1 = temp.bs;
for rr = 1:size(ei1,1)
    for cc = 1:size(ei1,2)
        this_b = ei1{rr,cc};
        if isempty(this_b)
            continue;
        end
        this_b.animal_id = temp.animalIDs(rr);
        this_b.belt_length = temp.belt_lengths(rr);
        ei1{rr,cc} = this_b;
    end
end
mData = evalin('base','mData');
colors = mData.colors;

selRows = [1:5]; selCols = [1:3];
% for ii = 1:length(selRows)
%     for jj = 1:3
%         ei(ii,jj) = ei1(selRows(ii),selColsAll(ii,jj));
%     end
% end
ei = ei1;
temp.animalIDs(selRows)
aids = temp.animalIDs(selRows);
% td = temp.training_days(selRows,selColsAll);
ii = 1;
moas = NaN(size(ei));
moasi = moas; mott = moas; motd = moas; momaxtd = moas;
pa7 = moas;
for iii = 1:size(ei,1)
    if iii == 2
        n = 0;
    end
    disp(temp.animalIDs(iii));
    out = behaviorProcessor_2(ei(iii,:));
    allbs{ii} = out;
    as{ii} = out.asT; mas{ii} = out.masT; semas{ii} = out.semasT;
    ml{ii} = out.mlT; mml{ii} = out.mmlT; semml{ii} = out.semmlT;
    mtt{ii} = out.mttT; semtt{ii} = out.semttT; mtd{ii} = out.mtdT; semtd{ii} = out.semtdT;
    maxtd{ii} = out.maxtdT;
    as1{ii} = out.asIT; mas1{ii} = out.masIT; semas1{ii} = out.semasIT;
    masD{ii} = out.masD; semasD{ii} = out.semasD; trialT{ii} = out.trial_times; trialD{ii} = out.trial_dist;
    ii = ii + 1;
    moas(iii,:) = mas{ii-1}; mott(iii,:) = mtt{ii-1}; motd(iii,:) = mtd{ii-1};
    momaxtd(iii,:) = maxtd{ii-1};
    moasi(iii,:) = mas1{ii-1};
    moml(iii,:) = mml{ii-1};
    pa7(iii,:) = out.percentAbove7;
end
num_trials = nan(size(moas));
for ii = 1:length(as)
    this_as = as{ii};
    for jj = 1:length(this_as)
        if length(this_as{jj})<20
            continue;
        end
        num_trials(ii,jj) = length(this_as{jj});
    end
end
selColsAll = [1 2 3
              1 2 3
              1 2 3
              1 2 3
              1 2 3
              ];
for ii = 1:size(num_trials,1)
    if isnan(selColsAll(ii,1))
        continue;
    end
    disp(num_trials(ii,selColsAll(ii,:)));
end


n = 0;
%% peercent of trials with greater than 7 cm/sec
runthis = 1;
if runthis
    thesecolors = distinguishable_colors(size(pa7,1),'w');
    ff = makeFigureWindow__one_axes_only(5,[10 4 2 1],[0.19 0.2 0.79 0.75]);
    axes(ff.ha);hold on;
    for ii = 1:size(pa7,1)
        if ismember(ii,[2])
            continue;
        end
        cols_to_pick = ~isnan(num_trials(ii,:));
        ys = pa7(ii,cols_to_pick);
        xs = temp.training_days(ii,cols_to_pick)+1;
        if xs(1) > 1
            xs = xs - xs(1)+1;
        end
%         temppos = find(xs == 1);
%         lastpos = temppos(find(temppos>1,1,'first'))-1;
%         xs = xs(1:lastpos);
%         ys = ys(1:lastpos);
        plot(xs,ys,'.-','color',thesecolors(ii,:));
    end
    set(gca,'xlim',[0 12.5],'ylim',[-5 101],'FontSize',6,'FontWeight','Bold','TickDir','out');
    set(gca,'xtick',1:13);
    xlabel('Training Day');
    ylabel({'Percent of trials','above 7cm/sec'});
    changePosition(gca,[0.02 0.06 -0.05 -0.06]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('percenta7.pdf'),600);
    return;
end

%%
runthis =1;
if runthis
thisCols_all = mData.colors;
    selRowi = 5;
for selRowi = 1:11
    selRowi
    these_cols = selColsAll(selRowi,:)
    ass = as1{selRowi};
    ff = makeFigureWindow__one_axes_only(5,[10 4 1.25 1],[0.19 0.2 0.79 0.75]);
    axes(ff.ha);hold on;
%     ass = as1{selRowi};
    for ii = 1:length(these_cols)
        if isnan(these_cols(ii))
            continue;
        end
        this = ass{these_cols(ii)};
        plot(1:length(this),this,'linewidth',0.5,'color',thisCols_all{ii});
        lenTs(ii) = length(this);
    end
    plot(1:max(lenTs),ones(size(1:max(lenTs)))*7,'m','linewidth',0.5);
    set(gca,'xlim',[0 max(lenTs)],'ylim',[0 30],'FontSize',6,'FontWeight','Bold','TickDir','out');
    changePosition(gca,[0.03 0.09 -0.03 -0.1])
    put_axes_labels(gca,{'Intertrial Number',[0 0 0]},{'Speed (cm/sec)',[0 0 0]});
    legs = [];
    for ii = 1:length(these_cols)
        legs{ii} = sprintf('Day %1d',ii);
    end
    legs{ii+1} = [5 3 30 4];
    putLegendH(ff.ha,legs,'colors',mData.colors,'sigR',{[],'anova',[],5});
    legs = {sprintf('Animal %d',temp.animalIDs(selRows(selRowi))),[22 0 25 4]};
%     putLegend(ff.ha,legs,'colors',{'k'});
%     title(sprintf('Animal %d',temp.animalIDs(selRows(selRowi))));
    changePosition(gca,[0 -0.05 0 0]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('Figure_1_Speed_vs_InterTrials_Training_%d.pdf',temp.animalIDs(selRows(selRowi))),600);
end
return;
end

%%
runthis = 1;
if runthis
this_moas = get_sel_values(moas,[1:5],selColsAll);
this_moasi = get_sel_values(moasi,[1:5],selColsAll);
% this_moas = moas([1 5:11],1:3); this_moas(3,:) = moas(6,[2 3 4]);
% this_moasi = moasi([1 5:11],1:3); this_moasi(3,:) = moasi(6,[2 3 4]);
for ii = 1:size(this_moas,2)
    varNames{ii} = sprintf('Trials_Day%d',ii);
end
for ii = 1:size(this_moasi,2)
    varNamesI{ii} = sprintf('InterTrials_Day%d',ii);
end
data = [this_moas this_moasi];
dataT = table(data(:,1),data(:,4),data(:,2),data(:,5),data(:,3),data(:,6));
dataT.Properties.VariableNames = {varNames{1} varNamesI{1} varNames{2} varNamesI{2} varNames{3} varNamesI{3}};
within = table([1 1 2 2 3 3]',[1 2 1 2 1 2]');
within.Properties.VariableNames = {'Day','TI'};
within.TI = categorical(within.TI);
within.Day = categorical(within.Day);

ra = repeatedMeasuresAnova(dataT,within);

% writetable(between,'Training_Data.xls');
rm = fitrm(dataT,'Trials_Day1,InterTrials_Day1,Trials_Day2,InterTrials_Day2,Trials_Day3,InterTrials_Day3 ~ 1','WithinDesign',within,'WithinModel','Day*TI');
rtable = ranova(rm,'WithinModel',rm.WithinModel);
mauchlytbl = mauchly(rm);
% multcompare(rm,'Day','ComparisonType','bonferroni')
mcTI = find_sig_mctbl(multcompare(rm,'TI','By','Day','ComparisonType','bonferroni'),6);
mcDays = find_sig_mctbl(multcompare(rm,'Day','By','TI','ComparisonType','bonferroni'),6);

mVar = ra.est_marginal_means.Mean;
semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
xdata = [1 2 4 5 7 8]; maxY = 50;

hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
hold on;
tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3}};
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'maxY',maxY,'ySpacing',6,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.7,'sigLinesStartYFactor',0.1);
for ii = 2:2:length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
end
set(gca,'xlim',[0.25 8.75],'ylim',[0 maxY+7],'FontSize',6,'FontWeight','Bold','TickDir','out');
xticks = [1.5 4.5 7.5]; xticklabels = {'Day1','Day2','Day3'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
changePosition(gca,[0.1 0.02 -0.03 -0.011])
put_axes_labels(gca,{[],[0 0 0]},{'Speed (cm/sec)',[0 0 0]});
rectangle(gca,'Position',[0.75 maxY+3 1 3],'edgecolor','k','facecolor','k');
text(1.85,maxY+5,'Trials','FontSize',5);
rectangle(gca,'Position',[4 maxY+3 1 3],'edgecolor','k');
text(5.2,maxY+5,'Inter-Trials','FontSize',5);

save_pdf(hf,mData.pdf_folder,'Figure_1_behavior_anova.pdf',600);
return;
end

%% check speed differences between groups
runthis = 1;
if runthis
    selColsAll1 = selColsAll;
this_moas = get_sel_values(moas,[1:5],selColsAll1(:,1:2));
this_moasi = get_sel_values(moasi,[1:5],selColsAll1(:,1:2));
% this_moas = moas([1 5:11],1:3); this_moas(3,:) = moas(6,[2 3 4]);
% this_moasi = moasi([1 5:11],1:3); this_moasi(3,:) = moasi(6,[2 3 4]);
for ii = 1:size(this_moas,2)
    varNames{ii} = sprintf('Trials_Day%d',ii);
end
for ii = 1:size(this_moasi,2)
    varNamesI{ii} = sprintf('InterTrials_Day%d',ii);
end
data0 = [this_moas this_moasi];
data = data0(:,[1 3 2 4]);
dataT = array2table([[ones(3,1);(2*ones(8,1))] data]);
dataT.Properties.VariableNames = {'Group' varNames{1} varNamesI{1} varNames{2} varNamesI{2}};
dataT.Group = categorical(dataT.Group);
within = table([1 1 2 2]',[1 2 1 2]');
within.Properties.VariableNames = {'Day','TI'};
within.TI = categorical(within.TI);
within.Day = categorical(within.Day);

ra = repeatedMeasuresAnova(dataT,within);

% writetable(between,'Training_Data.xls');
rm = fitrm(dataT,'Trials_Day1,InterTrials_Day1,Trials_Day2,InterTrials_Day2,Trials_Day3,InterTrials_Day3 ~ 1','WithinDesign',within,'WithinModel','Day*TI');
rtable = ranova(rm,'WithinModel',rm.WithinModel);
mauchlytbl = mauchly(rm);
% multcompare(rm,'Day','ComparisonType','bonferroni')
mcTI = find_sig_mctbl(multcompare(rm,'TI','By','Day','ComparisonType','bonferroni'),6);
mcDays = find_sig_mctbl(multcompare(rm,'Day','By','TI','ComparisonType','bonferroni'),6);

mVar = ra.est_marginal_means.Mean;
semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
xdata = [1 2 4 5 7 8]; maxY = 50;

hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
hold on;
tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3}};
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'maxY',maxY,'ySpacing',6,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.7,'sigLinesStartYFactor',0.1);
for ii = 2:2:length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
end
set(gca,'xlim',[0.25 8.75],'ylim',[0 maxY+7],'FontSize',6,'FontWeight','Bold','TickDir','out');
xticks = [1.5 4.5 7.5]; xticklabels = {'Day1','Day2','Day3'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
changePosition(gca,[0.1 0.02 -0.03 -0.011])
put_axes_labels(gca,{[],[0 0 0]},{'Speed (cm/sec)',[0 0 0]});
rectangle(gca,'Position',[0.75 maxY+3 1 3],'edgecolor','k','facecolor','k');
text(1.85,maxY+5,'Trials','FontSize',5);
rectangle(gca,'Position',[4 maxY+3 1 3],'edgecolor','k');
text(5.2,maxY+5,'Inter-Trials','FontSize',5);

save_pdf(hf,mData.pdf_folder,'Figure_1_behavior_anova_speed_groups.pdf',600);
return;
end

%%
runthis = 1;
if runthis
% this_moas = pa7([1 5:11],1:3); this_moas(3,:) = pa7(6,[2 3 4]);
this_moas = get_sel_values(pa7,[1:5],selColsAll);
for ii = 1:size(this_moas,2)
    varNames{ii} = sprintf('Trials_Day%d',ii);
end
data = [this_moas];
dataT = table(data(:,1),data(:,2),data(:,3));
dataT.Properties.VariableNames = varNames;
within = table([1 2 3]');
within.Properties.VariableNames = {'Day'};
within.Day = categorical(within.Day);
ra = repeatedMeasuresAnova(dataT,within);

mVar = ra.est_marginal_means.Mean;
semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
xdata = [1 2 3]; maxY = 50;

hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
hold on;
tcolors = {colors{1};colors{2};colors{3}};
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'maxY',maxY,'ySpacing',15,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);

set(gca,'xlim',[0.25 3.75],'ylim',[0 maxY+7],'FontSize',6,'FontWeight','Bold','TickDir','out');
xticks = xdata; xticklabels = {'Day1','Day2','Day3'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
% xtickangle(30)
changePosition(gca,[0.2 0.02 -0.15 -0.011])
put_axes_labels(gca,{[],[0 0 0]},{{'Percent trials','above 7cm/sec'},[0 0 0]});

save_pdf(hf,mData.pdf_folder,'Figure_1_behavior_anova_perc_trials.pdf',600);
return;
end

%%
runthis = 1;
if runthis
% this_moas = motd([1 5:11],1:3); this_moas(3,:) = motd(6,[2 3 4]);
this_moas = get_sel_values(motd,[1:5],selColsAll);
for ii = 1:size(this_moas,2)
    varNames{ii} = sprintf('Trials_Day%d',ii);
end
data = [this_moas];
dataT = table(data(:,1),data(:,2),data(:,3));
dataT.Properties.VariableNames = varNames;
within = table([1 2 3]');
within.Properties.VariableNames = {'Day'};
within.Day = categorical(within.Day);
ra = repeatedMeasuresAnova(dataT,within);

mVar = ra.est_marginal_means.Mean;
semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
xdata = [1 2 3]; maxY = 50;

hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
hold on;
tcolors = {colors{1};colors{2};colors{3}};
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'maxY',maxY,'ySpacing',20,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);

set(gca,'xlim',[0.25 3.75],'ylim',[0 maxY+3],'FontSize',6,'FontWeight','Bold','TickDir','out');
xticks = xdata; xticklabels = {'Day1','Day2','Day3'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
% xtickangle(30)
changePosition(gca,[0.2 0.02 -0.15 -0.011])
put_axes_labels(gca,{[],[0 0 0]},{{'Trial Distance','(% of belt length)'},[0 0 0]});

save_pdf(hf,mData.pdf_folder,'Figure_1_behavior_anova_trial_distance.pdf',600);
return;
end


%%
runthis = 1;
if runthis
% this_moas = motd([1 5:11],1:3); this_moas(3,:) = motd(6,[2 3 4]);
this_moas = get_sel_values(momaxtd,[1:5],selColsAll);
for ii = 1:size(this_moas,2)
    varNames{ii} = sprintf('Trials_Day%d',ii);
end
data = [this_moas];
dataT = table(data(:,1),data(:,2),data(:,3));
dataT.Properties.VariableNames = varNames;
within = table([1 2 3]');
within.Properties.VariableNames = {'Day'};
within.Day = categorical(within.Day);
ra = repeatedMeasuresAnova(dataT,within);

mVar = ra.est_marginal_means.Mean;
semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
xdata = [1 2 3]; maxY = 50;

hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
hold on;
tcolors = {colors{1};colors{2};colors{3}};
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'maxY',maxY,'ySpacing',20,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);

set(gca,'xlim',[0.25 3.75],'ylim',[0 maxY+3],'FontSize',6,'FontWeight','Bold','TickDir','out');
xticks = xdata; xticklabels = {'Day1','Day2','Day3'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
% xtickangle(30)
changePosition(gca,[0.22 0.02 -0.15 -0.011])
put_axes_labels(gca,{[],[0 0 0]},{{'Max Trial Distance','(% of belt length)'},[0 0 0]});

save_pdf(hf,mData.pdf_folder,'Figure_1_behavior_anova_max_trial_distance.pdf',600);
return;
end

%%
runthis = 1;
if runthis
% this_moas = mott([1 5:11],1:3); this_moas(3,:) = mott(6,[2 3 4]);
this_moas = get_sel_values(mott,[1:5],selColsAll);
for ii = 1:size(this_moas,2)
    varNames{ii} = sprintf('Trials_Day%d',ii);
end
data = [this_moas];
dataT = table(data(:,1),data(:,2),data(:,3));
dataT.Properties.VariableNames = varNames;
within = table([1 2 3]');
within.Properties.VariableNames = {'Day'};
within.Day = categorical(within.Day);
ra = repeatedMeasuresAnova(dataT,within);

mVar = ra.est_marginal_means.Mean;
semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
xdata = [1 2 3]; maxY = 50;

hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
hold on;
tcolors = {colors{1};colors{2};colors{3}};
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'maxY',maxY,'ySpacing',20,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);

set(gca,'xlim',[0.25 3.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
xticks = xdata; xticklabels = {'Day1','Day2','Day3'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
% xtickangle(30)
changePosition(gca,[0.1 0.02 -0.15 -0.011])
put_axes_labels(gca,{[],[0 0 0]},{{'Trial Time (sec)'},[0 0 0]});

save_pdf(hf,mData.pdf_folder,'Figure_1_behavior_anova_trial_time.pdf',600);
return;
end

%%
runthis = 1;
if runthis
% this_moas = moml([1 5:11],1:3); this_moas(3,:) = moml(6,[2 3 4]);
this_moas = get_sel_values(moml,[1:5],selColsAll);
for ii = 1:size(this_moas,2)
    varNames{ii} = sprintf('Trials_Day%d',ii);
end
data = [this_moas];
dataT = table(data(:,1),data(:,2),data(:,3));
dataT.Properties.VariableNames = varNames;
within = table([1 2 3]');
within.Properties.VariableNames = {'Day'};
within.Day = categorical(within.Day);
ra = repeatedMeasuresAnova(dataT,within);

mVar = ra.est_marginal_means.Mean;
semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
xdata = [1 2 3]; maxY = 50;

hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
hold on;
tcolors = {colors{1};colors{2};colors{3}};
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'maxY',maxY,'ySpacing',20,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);

set(gca,'xlim',[0.25 3.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
xticks = xdata; xticklabels = {'Day1','Day2','Day3'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
% xtickangle(30)
changePosition(gca,[0.12 0.02 -0.15 -0.011])
put_axes_labels(gca,{[],[0 0 0]},{{'Movement','Latency (sec)'},[0 0 0]});

save_pdf(hf,mData.pdf_folder,'Figure_1_behavior_anova_mov_lat.pdf',600);
return;
end

%%
runthis = 1;
if runthis
% this_moas = moml([1 5:11],1:3); this_moas(3,:) = moml(6,[2 3 4]);
this_moas = get_sel_values(temp.weight,[1:5],repmat([1 2 3 4],5,1));
this_moas = 100 * this_moas./this_moas(:,1);
for ii = 1:size(this_moas,2)
    varNames{ii} = sprintf('Trials_Day%d',ii);
end
data = [this_moas];
dataT = array2table(data);
dataT.Properties.VariableNames = varNames;
within = table([1 2 3 4]');
within.Properties.VariableNames = {'Day'};
within.Day = categorical(within.Day);
ra = repeatedMeasuresAnova(dataT,within);

mVar = ra.est_marginal_means.Mean;
semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
xdata = [1 2 3 5]; maxY = 50;

hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
hold on;
tcolors = {colors{1};colors{2};colors{3};colors{4}};
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);

set(gca,'xlim',[0.25 5.75],'ylim',[95 maxY-10],'FontSize',6,'FontWeight','Bold','TickDir','out');
xticks = xdata; xticklabels = {'Day1','Day2','Day3','Test-Day'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
xtickangle(20)
changePosition(gca,[0.2 0.02 -0.15 -0.011])
put_axes_labels(gca,{[],[0 0 0]},{{'Percentage of','weight on day 1'},[0 0 0]});

save_pdf(hf,mData.pdf_folder,'Figure_1_behavior_anova_weight.pdf',600);
return;
end


function out = behaviorProcessor_2(ei)
as = cell(1,length(ei)); lenT = NaN(1,length(ei)); mas = lenT; semas = lenT; percentAbove7 = lenT;
mtt = mas; mtd = mas; semtt = semas; semtd = semas; maxtd = mas;
ml = as; mml = mas; semml = semas;
as1 = as; lenT1 = lenT; mas1 = mas; semas1 = semas; 
asd = as; lenTd = lenT; masd = mas; semasd = semas; 

for ii = 1:length(ei)
    b = ei{ii};
    if isempty(b)
        continue;
    end
%     figure(101);clf;
%     plot(b.ts,b.air_puff_raw,'color','b','linewidth',1.5);hold on;
%     plot(b.ts,0.25*b.photo_sensor_raw+0.55,'color','r','linewidth',1.5);
%     plot(b.ts,0.5*b.fSpeed/max(b.fSpeed),'color','m','linewidth',1.5);
%     xlabel('Time (sec)');
%     set(gca,'FontSize',12,'FontWeight','Bold');
    [as{ii},tt{ii},td{ii},ml{ii}] = findAverageSpeedTrials(b);
    lenT(ii) = length(as{ii});
    mas(ii) = mean(as{ii}); 
    mtt(ii) = mean(tt{ii}); 
    mtd(ii) = mean(td{ii});
    maxtd(ii) = max(td{ii});
    mml(ii) = mean(ml{ii});
    semas(ii) = std(as{ii})/sqrt(lenT(ii)); 
    semtt(ii) = std(tt{ii})/sqrt(lenT(ii)); 
    semtd(ii) = std(td{ii})/sqrt(lenT(ii));
    semml(ii) = std(ml{ii})/sqrt(lenT(ii));
    percentAbove7(ii) = 100*sum(as{ii}>7)/sum(~isnan(as{ii}));
    
    as1{ii} = findAverageSpeedInterTrials(b);
    lenT1(ii) = length(as1{ii});
    mas1(ii) = mean(as1{ii});
    semas1(ii) = std(as1{ii})/sqrt(lenT1(ii));
    
    asd{ii} = findDiffTrialsInterTrials(b);
    lenTd(ii) = length(asd{ii});
    masd(ii) = mean(asd{ii});
    semasd(ii) = std(asd{ii})/sqrt(lenTd(ii));
    n = 0;
end
out.asT = as; out.lenT = lenT; out.masT = mas; out.semasT = semas;
out.asIT = as1; out.lenIT = lenT1; out.masIT = mas1; out.semasIT = semas1;
out.masD = masd; out.semasD = semasd;
out.percentAbove7 = percentAbove7; out.trial_times = tt; out.trial_dist = td;
out.mttT = mtt; out.maxtdT = maxtd; out.mtdT = mtd; out.semttT = semtt; out.semtdT = semtd;
out.mlT = ml; out.mmlT = mml; out.semmlT = semml;
n = 0;


function as = findDiffTrialsInterTrials (b)
if isfield(b,'stim_f')
    b.air_puff_r = b.stim_r;
    b.air_puff_f = b.stim_f;
end
as = [];
for ii = 1:(length(b.air_puff_r)-1)
    speeds = b.fSpeed(b.air_puff_r(ii):b.air_puff_f(ii));
    if sum(speeds<0) > 0
        speeds(speeds < 0) = NaN;
    end
    ass = nanmean(speeds);
    speeds = b.fSpeed(b.air_puff_f(ii):b.air_puff_r(ii+1));
    if sum(speeds<0) > 0
        speeds(speeds < 0) = NaN;
    end
    assi = nanmean(speeds);
    as(ii) = ass - assi;
end


function [as,trial_time,trial_dist,mov_lat] = findAverageSpeedTrials (b)
if isfield(b,'stim_f')
    b.air_puff_r = b.stim_r;
    b.air_puff_f = b.stim_f;
end
% if length(b.air_puff_r) < 30
%     as = NaN;
%     return;
% end
as = [];
BL = b.belt_length;
for ii = 1:length(b.air_puff_r)
    speeds = b.fSpeed(b.air_puff_r(ii):b.air_puff_f(ii));
    trial_time(ii) = b.ts(b.air_puff_f(ii)) - b.ts(b.air_puff_r(ii));
    trial_dist(ii) = round(100*(b.dist(b.air_puff_f(ii)) - b.dist(b.air_puff_r(ii)))/BL);
    if sum(speeds<0) > 0
        speeds(speeds < 0) = NaN;
%         speeds = fillmissing(speeds,'linear',2,'EndValues','nearest');
        n = 0;
    end
    ind = find(speeds > 0,1,'first');
    if isempty(ind)
        mov_lat(ii) = NaN;
    else
        mov_lat(ii) = b.ts(b.air_puff_r(ii)+ind)-b.ts(b.air_puff_r(ii));
    end
    as(ii) = nanmean(speeds);
end
n = 0;

function as = findAverageSpeedInterTrials (b)
if isfield(b,'stim_f')
    b.air_puff_r = b.stim_r;
    b.air_puff_f = b.stim_f;
end
as = [];
for ii = 1:(length(b.air_puff_r)-1)
    speeds = b.fSpeed(b.air_puff_f(ii):b.air_puff_r(ii+1));
    if sum(speeds<0) > 0
        speeds(speeds < 0) = NaN;
%         speeds = fillmissing(speeds,'linear',2,'EndValues','nearest');
        n = 0;
    end
    as(ii) = nanmean(speeds);
end


function numberOfTrials = findNumberOfTrials(training_data)
bs = training_data.bs;
numberOfTrials = NaN(size(bs));
for ii = 1:size(bs,1)
    for jj = 1:size(bs,2)
        thisb = bs{ii,jj};
        if isempty(thisb)
            continue;
        end
        numberOfTrials(ii,jj) = length(thisb.air_puff_r);
    end
end

function vals = get_sel_values(var,selrows,selcols)
for ii = 1:length(selrows)
    iii = selrows(ii);
    vals(ii,:) = var(iii,selcols(iii,:));
end