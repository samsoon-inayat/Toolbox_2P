function behaviorProcessor_1_AD_Ctrl
mData = evalin('base','mData');
colors = mData.colors;
data_C = get_training_data_C;
data_A = get_training_data_A;
n = 0;
% %%
% xl_fileName = 'Training_data_C_A.xlsx';
% [numbers, strings, raw] = xlsread(xl_fileName);
% [within,dvn,xlabels] = make_within_table({'Day','Trials'},[3,2]);
% dataT = make_between_table({var_C;var_A},dvn);
% ra = RMA(dataT,within,{0.05,{'hsd'}});
% ra.ranova
% print_for_manuscript(ra)
% n = 0;

%%
moas = data_C.moas;
moasi = data_C.moasi;
moas_A = data_A.moas;
moasi_A = data_A.moasi;
var_C = [moas(:,1) moasi(:,1) moas(:,2) moasi(:,2) moas(:,3) moasi(:,3)]; 
var_A = [moas_A(:,1) moasi_A(:,1) moas_A(:,2) moasi_A(:,2) moas_A(:,3) moasi_A(:,3)];


var_C = [moas(:,1) moasi(:,1) moas(:,2) moasi(:,2) moas(:,3) moasi(:,3)]; var_A = [moas_A(:,1) moasi_A(:,1) moas_A(:,2) moasi_A(:,2) moas_A(:,3) moasi_A(:,3)];
[within,dvn,xlabels] = make_within_table({'Day','Trials'},[3,2]);
dataT = make_between_table({var_C;var_A},dvn);
ra = RMA(dataT,within,{0.05,{'bonferroni','hsd','lsd'}});
ra.ranova
print_for_manuscript(ra)

%
redF = [2]; redV = {2};
[dataTR,withinR] = reduce_within_between(dataT,within,redF,redV);
raR = RMA(dataTR,withinR,{0.05,{'hsd'}});
raR.ranova
print_for_manuscript(raR)
n = 0;
%%
    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[10 3 3.5 1],'RowsCols',[1 2],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.35],'widthHeightAdjustment',[10 -410]);
    MY = 47; ysp = 0.5; mY = 0; titletxt = ''; ylabeltxt = {'Average','Speed (cm/s)'};
    stp = 0.36*magfac; widths = ([1.9 1 1.3 1.3 1.3 0.5 0.5 0.5])*magfac; gap = 0.16*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{''});

   [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_Day_Trials','hsd'},[1.5 1 1]);
    xdata = make_xdata([2 2 2 2 2 2],[1 1.5]);   
    axes(ff.h_axes(1,1));    hold on;
    tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3}};
    tcolors = repmat(tcolors,2,1)';
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,[],[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);

    set(gca,'xlim',[0.25 14.75],'ylim',[0 MY],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = [1.5 3.5 5.5 9.5 11.5 13.5]; xticklabels = {'Day1','Day2','Day3'}; xticklabels = repmat(xticklabels,2,1)';
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    for ii = 2:2:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,6,{'C-TG','A-TG'});
%     rectangle(gca,'Position',[0.75 30.5 1 3.5],'edgecolor','k','facecolor','k');
%     text(1.85,30.5,'Trials','FontSize',5);
%     rectangle(gca,'Position',[6 30.5 1 3.5],'edgecolor','k');
%     text(7.2,30.5,'Inter-Trials','FontSize',5);
% For pooled
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Day_by_Trials','lsd'},[1.5 1 1]);
    xdata = make_xdata([2 2 2],[1 1.5]);   
    axes(ff.h_axes(1,2));
    hold on;
    colors = mData.colors;
    tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3}};
    tcolors = repmat(tcolors,2,1)';
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);

    set(gca,'xlim',[0.25 max(xdata)+0.75],'ylim',[0 MY],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = [1.5 4.5 7.5]; xticklabels = {'Day1','Day2','Day3'}; 
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    for ii = 2:2:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,6,{'Pooled'});
    save_pdf(ff.hf,mData.pdf_folder,'Figure_1_behavior_anova_speed.pdf',600);
%%
runthis =0;
if runthis
thisCols_all = mData.colors;
    selRowi = 5;
for selRowi = 1:6
    selRowi
%     ass = as{selRowi};
    ff = makeFigureWindow__one_axes_only(5,[10 4 1.25 1.25],[0.19 0.2 0.79 0.75]);
    axes(ff.ha);hold on;
    ass = as1{selRowi};
    for ii = 1:length(ass)
        this = ass{ii};
        plot(1:length(this),this,'linewidth',0.5,'color',thisCols_all{ii});
        lenTs(ii) = length(this);
    end
    plot(1:max(lenTs),ones(size(1:max(lenTs)))*7,'m','linewidth',0.5);
    set(gca,'xlim',[0 max(lenTs)],'ylim',[0 30],'FontSize',6,'FontWeight','Bold','TickDir','out');
    changePosition(gca,[0.03 0.09 -0.03 -0.1])
    put_axes_labels(gca,{'Inter-Trial Number',[0 0 0]},{'Speed (cm/sec)',[0 0 0]});
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
    save_pdf(ff.hf,mData.pdf_folder,sprintf('Figure_1_Speed_vs_InterTrials_Training_%d.pdf',temp.animalIDs(selRows(selRowi))),600);
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
n = 0;
moas = data_C.moas;
moasi = data_C.moasi;
moas_A = data_A.moas;
moasi_A = data_A.moasi;
for ii = 1:size(moas,2)
    varNames{ii} = sprintf('Trials_Day%d',ii);
end
for ii = 1:size(moasi,2)
    varNamesI{ii} = sprintf('InterTrials_Day%d',ii);
end
data_C = [moas moasi];
data_A = [moas_A moasi_A];
group = categorical([ones(size(data_C,1),1);(2*ones(size(data_A,1),1))]);
data = [data_C;data_A];
dataT = table(group,data(:,1),data(:,4),data(:,2),data(:,5),data(:,3),data(:,6));
dataT.Properties.VariableNames = {'Group' varNames{1} varNamesI{1} varNames{2} varNamesI{2} varNames{3} varNamesI{3}};
writetable(dataT,'Training_data_C_A.xlsx')
within = table([varNames';varNamesI']);
columnText = cell(size(within,1),1);columnText(1:2:end)= varNames';columnText(2:2:end)= varNamesI';
within = table([varNames';varNamesI'],columnText);
within = table([1 1 2 2 3 3]',[1 2 1 2 1 2]');
within.Properties.VariableNames = {'Day','TI'};
within.TI = categorical(within.TI);
within.Day = categorical(within.Day);

% writetable(between,'Training_Data.xls');
rm = fitrm(dataT,'Trials_Day1,InterTrials_Day1,Trials_Day2,InterTrials_Day2,Trials_Day3,InterTrials_Day3 ~ Group','WithinDesign',within,'WithinModel','Day*TI');
rtable = ranova(rm,'WithinModel',rm.WithinModel);
mauchlytbl = mauchly(rm);
% multcompare(rm,'Day','ComparisonType','bonferroni')
mcGroup = find_sig_mctbl(multcompare(rm,'Group','By','Day','ComparisonType','bonferroni'),6);
mcTI = find_sig_mctbl(multcompare(rm,'TI','By','Day','ComparisonType','bonferroni'),6);
mcDays = find_sig_mctbl(multcompare(rm,'Day','By','TI','ComparisonType','bonferroni'),6);

[mVarT,semVarT] = findMeanAndStandardError(moas);
[mVarIT,semVarIT] = findMeanAndStandardError(moasi);

[mVarT_A,semVarT_A] = findMeanAndStandardError(moas_A);
[mVarIT_A,semVarIT_A] = findMeanAndStandardError(moasi_A);

mVar = NaN(1,2*(length(mVarT)+length(mVarIT)));
semVar = mVar;
mVar(1:2:6) = mVarT;semVar(1:2:6) = semVarT;
mVar(2:2:6) = mVarIT;semVar(2:2:6) = semVarIT;
mVar(7:2:12) = mVarT_A;semVar(7:2:12) = semVarT_A;
mVar(8:2:12) = mVarIT_A;semVar(8:2:12) = semVarIT_A;
combs = nchoosek(1:14,2); p = ones(size(combs,1),1); h = logical(zeros(size(combs,1),1));
% row = [7 8]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{1,6}; h(ii) = 1; 
% row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
% row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

xdata = [1 2 4 5 7 8 [13 14 16 17 19 20]-1]; maxY = 22;
colors = mData.colors;
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 2.25 1],'color','w');
hold on;
tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{1};colors{1};colors{2};colors{2};colors{3};colors{3}};
hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'maxY',maxY,'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',-0.1);
for ii = 2:2:length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
end
% plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+1],'FontSize',6,'FontWeight','Bold','TickDir','out');
xticks = xdata(1:2:end)+0.5; xticklabels = {'D1','D2','D3','D1','D2','D3'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
changePosition(gca,[0.02 0.03 0.02 -0.11]);
put_axes_labels(gca,{[],[0 0 0]},{'Avg. Speed (cm/sec)',[0 0 0]});
rectangle(gca,'Position',[0.75 21 1 2],'edgecolor','k','facecolor','k');
text(1.85,22,'Trials','FontSize',5);
rectangle(gca,'Position',[6 21 1 2],'edgecolor','k');
text(7.2,22,'Inter-Trials','FontSize',5);
text(3.5,18,'Control','FontSize',7);
text(18.5,18,'APP','FontSize',7);
% applyhatch_plusC(gcf

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


function out1 = get_training_data_C
temp = evalin('base','training_data_C');
numberOfTrials = findNumberOfTrials(temp);
[rr,cc] = find(numberOfTrials < 20);
% selColsAll = [1 2 3
%               1 3 4;
%               1 3 4;
%               1 2 3;
%               1 3 4;
%               1 2 3;
%               1 2 3;
%               1 2 3;
%               1 2 3;
%               1 2 3;
%               ];
ei1 = temp.bs;
mData = evalin('base','mData');
colors = mData.colors;
% selRows = [4 5 6 8 9]; selCols = [1:4];
selRows = [1 4 5 6 7 8:10]; selCols = [1:3];

for ii = 1:length(selRows)
    for jj = 1:3
        if ii == 3
            ei(ii,jj) = ei1(selRows(ii),jj+1);
        else
            ei(ii,jj) = ei1(selRows(ii),jj);
        end
%         ei(ii,jj) = ei1(selRows(ii),selColsAll(ii,jj));
    end
end
temp.animalIDs(selRows)
aids = temp.animalIDs(selRows);
td = temp.training_days(selRows,1:3);
ii = 1;
moas = NaN(size(ei));
moasi = moas;
for iii = 1:size(ei,1)
    if iii == 3
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
out1.moas = moas;
out1.moasi = moasi;
n = 0;

function out1 = get_training_data_A
temp = evalin('base','training_data_A');
numberOfTrials = findNumberOfTrials(temp);
[rr,cc] = find(numberOfTrials < 20);
% selColsAll = [1 2 3
%               1 2 3;
%               1 2 3;
%               1 2 3;
%               1 2 3;
%               1 2 3;
%               ];
ei1 = temp.bs;
mData = evalin('base','mData');
colors = mData.colors;
% selRows = [4 5 6 8 9]; selCols = [1:4];
selRows = [1 2 3 4:6]; selCols = [1:3];
for ii = 1:length(selRows)
    for jj = 1:3
        ei(ii,jj) = ei1(selRows(ii),jj);
%         ei(ii,jj) = ei1(selRows(ii),selColsAll(ii,jj));
    end
end
temp.animalIDs(selRows)
aids = temp.animalIDs(selRows);
td = temp.training_days(selRows,1:3);
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
out1.moas = moas;
out1.moasi = moasi;
n = 0;