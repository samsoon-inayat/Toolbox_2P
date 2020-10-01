function behaviorProcessor_1
temp = evalin('base','training_data_C');
numberOfTrials = findNumberOfTrials(temp);
[rr,cc] = find(numberOfTrials < 20);
selColsAll = [1 2 3
              1 3 4;
              1 3 4;
              1 2 3;
              1 3 4;
              1 2 3;
              1 2 3;
              1 2 3;
              1 2 3;
              1 2 3;
              ];
ei1 = temp.bs;
mData = evalin('base','mData');
colors = mData.colors;
% selRows = [4 5 6 8 9]; selCols = [1:4];
selRows = [1:6 8:10]; selCols = [1:3];
for ii = 1:length(selRows)
    for jj = 1:3
        ei(ii,jj) = ei1(selRows(ii),selColsAll(ii,jj));
    end
end
temp.animalIDs(selRows)
aids = temp.animalIDs(selRows);
td = temp.training_days(selRows,selColsAll);
ii = 1;
moas = NaN(size(ei));
moasi = moas;
for iii = 1:size(ei,1)
    if iii == 4
        n = 0;
    end
    out = behaviorProcessor_2(ei(iii,:));
    as{ii} = out.asT; mas{ii} = out.masT; semas{ii} = out.semasT;
    as1{ii} = out.asIT; mas1{ii} = out.masIT; semas1{ii} = out.semasIT;
    masD{ii} = out.masD; semasD{ii} = out.semasD;
    ii = ii + 1;
    moas(iii,:) = mas{ii-1};
    moasi(iii,:) = mas1{ii-1};
end
n = 0;
%%
runthis =1;
if runthis
thisCols_all = mData.colors;
    selRowi = 5;
for selRowi = 1:9
    selRowi
    ass = as{selRowi};
    ff = makeFigureWindow__one_axes_only(5,[10 4 1.25 1.25],[0.19 0.2 0.79 0.75]);
    axes(ff.ha);hold on;
%     ass = as1{selRowi};
    for ii = 1:length(ass)
        this = ass{ii};
        plot(1:length(this),this,'linewidth',0.5,'color',thisCols_all{ii});
        lenTs(ii) = length(this);
    end
    plot(1:max(lenTs),ones(size(1:max(lenTs)))*7,'m','linewidth',0.5);
    set(gca,'xlim',[0 max(lenTs)],'ylim',[0 30],'FontSize',6,'FontWeight','Bold','TickDir','out');
    changePosition(gca,[0.03 0.09 -0.03 -0.1])
    put_axes_labels(gca,{'Trial Number',[0 0 0]},{'Speed (cm/sec)',[0 0 0]});
    legs = [];
    for ii = 1:length(ass)
        legs{ii} = sprintf('Day %1d',ii);
    end
    legs{ii+1} = [5 3 30 4];
    putLegendH(ff.ha,legs,'colors',mData.colors,'sigR',{[],'anova',[],5});
    legs = {sprintf('Animal %d',temp.animalIDs(selRows(selRowi))),[22 0 25 4]};
%     putLegend(ff.ha,legs,'colors',{'k'});
    title(sprintf('Animal %d',temp.animalIDs(selRows(selRowi))));
    changePosition(gca,[0 -0.05 0 0]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('Figure_1_Speed_vs_Trials_Training_%d.pdf',temp.animalIDs(selRows(selRowi))),600);
end
return;
end


%%
runthis = 0;
if runthis
thisCols_all = mData.colors;
ff = makeFigureWindow__one_axes_only(5,[10 4 1.25 1],[0.19 0.2 0.79 0.75]);
axes(ff.ha);hold on;
% for ii = 1:length(mas)
%     plot(1:length(as{ii}),mas{ii},'linewidth',0.5,'color',thisCols_all{ii});
%     errorbar(1:length(as{ii}),mas{ii},semas{ii},'linewidth',0.25,'color',thisCols_all{ii},'CapSize',1);
% end
% for ii = 1:length(mas)
%     plot(1:length(as1{ii}),mas1{ii},'-.','linewidth',0.5,'color',thisCols_all{ii});
%     errorbar(1:length(as1{ii}),mas1{ii},semas1{ii},'--','linewidth',0.25,'color',thisCols_all{ii},'CapSize',1);
% end
for ii = 1:length(masD)
    plot(1:length(as{ii}),masD{ii},'linewidth',0.5,'color',thisCols_all{ii});
    errorbar(1:length(as{ii}),masD{ii},semasD{ii},'linewidth',0.25,'color',thisCols_all{ii},'CapSize',1);
end

set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out');
hxl = xlabel('Training Day'); changePosition(hxl,[0 -3 0]);
hyl = ylabel('Speed (cm/sec)');get(hyl,'Position');get(hyl,'Units')
changePosition(hyl,[0 7 0]);
get(hyl,'Position')
get(hyl,'Units')
changePosition(hyl,[-0.2 -1 0]);
get(hyl,'Position')
get(hyl,'Units')
xlim([0.75 3.25]);
ylim([0 20]);

legs = [];
for ii = 1:length(mas)
    legs{ii} = sprintf('Animal %1d',ii);
end
legs{ii+1} = [0.85 0 46 5.5];
% putLegend(ff.ha,legs,'colors',mData.colors,'sigR',{[],'anova',[],5});

% text(2.25,45,'-  Trials','FontSize',5)
% text(2.25,40,'-- InterTrials','FontSize',5)
changePosition(gca,[0.03 0.09 -0.03 -0.1])
save_pdf(ff.hf,mData.pdf_folder,'Figure_1_behavior.pdf',600);
return;
end

%%
runthis = 1;
if runthis
% tempTxt = {'Trials'}; TrialsInterTrials = repmat(tempTxt,size(moas,1),1);
% tempTxt = {'InterTrials'}; TrialsInterTrials = [TrialsInterTrials;repmat(tempTxt,size(moas,1),1)];
% for ii = 1:size(moas,2)
%     varNames{ii} = sprintf('Day%d',ii);
% end
% data = [moas;moasi];
% 
% between = table(TrialsInterTrials,data(:,1),data(:,2),data(:,3));
% between.Properties.VariableNames = {'TI','Day1','Day2','Day3'};
% within = table(varNames');
% within.Properties.VariableNames = {'Day'};
% 
% % writetable(between,'Training_Data.xls');
% 
% rm = fitrm(between,'Day1-Day3 ~ TI','WithinDesign',within,'WithinModel','Day');
% rtable = ranova(rm,'WithinModel',rm.WithinModel);
% mauchlytbl = mauchly(rm);
% mcTI = find_sig_mctbl(multcompare(rm,'TI','By','Day','ComparisonType','bonferroni'),6);
% mcDays = find_sig_mctbl(multcompare(rm,'Day','By','TI','ComparisonType','bonferroni'),6);

for ii = 1:size(moas,2)
    varNames{ii} = sprintf('Trials_Day%d',ii);
end
for ii = 1:size(moasi,2)
    varNamesI{ii} = sprintf('InterTrials_Day%d',ii);
end
data = [moas moasi];
dataT = table(data(:,1),data(:,4),data(:,2),data(:,5),data(:,3),data(:,6));
dataT.Properties.VariableNames = {varNames{1} varNamesI{1} varNames{2} varNamesI{2} varNames{3} varNamesI{3}};
within = table([varNames';varNamesI']);
columnText = cell(size(within,1),1);columnText(1:2:end)= varNames';columnText(2:2:end)= varNamesI';
within = table([varNames';varNamesI'],columnText);
within = table([1 1 2 2 3 3]',[1 2 1 2 1 2]');
within.Properties.VariableNames = {'Day','TI'};
within.TI = categorical(within.TI);
within.Day = categorical(within.Day);

% writetable(between,'Training_Data.xls');
rm = fitrm(dataT,'Trials_Day1,InterTrials_Day1,Trials_Day2,InterTrials_Day2,Trials_Day3,InterTrials_Day3 ~ 1','WithinDesign',within,'WithinModel','Day*TI');
rtable = ranova(rm,'WithinModel',rm.WithinModel);
mauchlytbl = mauchly(rm);
% multcompare(rm,'Day','ComparisonType','bonferroni')
mcTI = find_sig_mctbl(multcompare(rm,'TI','By','Day','ComparisonType','bonferroni'),6);
mcDays = find_sig_mctbl(multcompare(rm,'Day','By','TI','ComparisonType','bonferroni'),6);

[mVarT,semVarT] = findMeanAndStandardError(moas);
[mVarIT,semVarIT] = findMeanAndStandardError(moasi);
mVar = NaN(1,length(mVarT)+length(mVarIT));
semVar = mVar;
mVar(1:2:length(mVar)) = mVarT;semVar(1:2:length(mVar)) = semVarT;
mVar(2:2:length(mVar)) = mVarIT;semVar(2:2:length(mVar)) = semVarIT;
combs = nchoosek(1:6,2); p = ones(size(combs,1),1); h = logical(zeros(size(combs,1),1));
row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{1,6}; h(ii) = 1; 
row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 
row = [1 5]; ii = ismember(combs,row,'rows'); p(ii) = mcDays{2,6}; h(ii) = 1; 
row = [1 3]; ii = ismember(combs,row,'rows'); p(ii) = mcDays{1,6}; h(ii) = 1; 
% row = [2 6]; ii = ismember(combs,row,'rows'); p(ii) = mcDays{3,6}; h(ii) = 1; 

xdata = [1 2 4 5 7 8]; maxY = 50;

hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
hold on;
tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3}};
hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'maxY',maxY,'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0);
for ii = 2:2:length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
end
set(gca,'xlim',[0.25 8.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
xticks = [1.5 4.5 7.5]; xticklabels = {'Day1','Day2','Day3'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
changePosition(gca,[0.1 0.02 -0.03 -0.011])
put_axes_labels(gca,{[],[0 0 0]},{'Speed (cm/sec)',[0 0 0]});
rectangle(gca,'Position',[0.75 47 1 3],'edgecolor','k','facecolor','k');
text(1.85,49,'Trials','FontSize',5);
rectangle(gca,'Position',[4 47 1 3],'edgecolor','k');
text(5.2,49,'Inter-Trials','FontSize',5);

save_pdf(hf,mData.pdf_folder,'Figure_1_behavior_anova.pdf',600);
return;
end



function out = behaviorProcessor_2(ei)
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
    as{ii} = findAverageSpeedTrials(b);
    lenT(ii) = length(as{ii});
    mas(ii) = mean(as{ii});
    semas(ii) = std(as{ii})/sqrt(lenT(ii));
    
    as1{ii} = findAverageSpeedInterTrials(b);
    lenT1(ii) = length(as1{ii});
    mas1(ii) = mean(as1{ii});
    semas1(ii) = std(as1{ii})/sqrt(lenT1(ii));
    
    asd{ii} = findDiffTrialsInterTrials(b);
    lenT(ii) = length(asd{ii});
    masd(ii) = mean(asd{ii});
    semasd(ii) = std(asd{ii})/sqrt(lenT(ii));
    
    n = 0;
end
out.asT = as; out.lenT = lenT; out.masT = mas; out.semasT = semas;
out.asIT = as1; out.lenIT = lenT1; out.masIT = mas1; out.semasIT = semas1;
out.masD = masd; out.semasD = semasd;
n = 0;


function as = findDiffTrialsInterTrials (b)
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


function as = findAverageSpeedTrials (b)
as = [];
for ii = 1:length(b.air_puff_r)
    speeds = b.fSpeed(b.air_puff_r(ii):b.air_puff_f(ii));
    if sum(speeds<0) > 0
        speeds(speeds < 0) = NaN;
%         speeds = fillmissing(speeds,'linear',2,'EndValues','nearest');
        n = 0;
    end
    as(ii) = nanmean(speeds);
end
% n = 0;

function as = findAverageSpeedInterTrials (b)
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