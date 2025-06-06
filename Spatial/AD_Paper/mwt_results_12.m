function mwt_results

mData = evalin('base','mData');
colors = mData.colors;

filename = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\PDFs\APP Figure 1 MWT Data.xlsx';

marker_styles = {'s','d','^','s','d','^','^','*'};
line_styles = {'-',':','-.','-',':','-.'};
mSize = 5;
tcolors = {'k','r'};
n=0;
%% 12months speed
if 0
[num,strings,raw] = xlsread(filename,1,'A1:I24');
% data = array2table(num);
% data.Properties.VariableNames = strings;
% data.Group = categorical(data.Group);
% within = [1:8]';
% within = array2table(within);
% within.Properties.VariableNames = {'Day'};
% within.Day = categorical(within.Day);
% ra = repeatedMeasuresAnova(data,within,0.05);

[within,dvn,xlabels] = make_within_table({'Day'},[8]);
dataT = make_between_table({num(num(:,1)==1,2:end);num(num(:,1)==2,2:end)},dvn);
ra = RMA(dataT,within,{0.05,{'bonferroni'}});
ra.ranova
print_for_manuscript(ra)

n = 0;


    mVar = ra.est_marginal_means.Mean;
    semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < -0.05;
    xdata = 1:8;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.9 1],'color','w');
    hold on;
    plot(xdata,mVar(1:8),'.','color',tcolors{1},'markersize',mSize,'marker','.','linestyle','-');
    errorbar(xdata,mVar(1:8),semVar(1:8),'.','color',tcolors{1},'markersize',mSize);
    plot(xdata,mVar(9:16),'.','color',tcolors{2},'markersize',mSize,'marker','.','linestyle','-');
    errorbar(xdata,mVar(9:16),semVar(9:16),'.','color',tcolors{2},'markersize',mSize);

    set(gca,'xlim',[0 9],'ylim',[0 0.25],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = [1:8]; %%xticklabels = {'Day1','Day2','Day3'}; xticklabels = repmat(xticklabels,2,1)';
    set(gca,'xtick',xticks);
%     xtickangle(30);
    changePosition(gca,[0.08 0.13 -0.01 -0.1])
    put_axes_labels(gca,{'Days',[0 0 0]},{{'Average','Speed (m/s)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'Figure_1_mwt_anova_speed.pdf',600);
return;
end

%% 12months latency
if 1
[num,strings,raw] = xlsread(filename,2,'A1:I24');
txt_filename = fullfile('E:\PostProcessing\AD_paper','MWT.txt');
adata = make_text_file_for_nlme_R(num,txt_filename);
group = categorical(adata(:,4));
dv = dummyvar(group);
X = [adata(:,2) dv]; y = adata(:,3);
modelfun = @(b,x)X(:,2).*(b(1).*2.^b(2)*X(:,1)) + X(:,3).*(b(3).*2.^b(4)*X(:,1));
beta0 = [60 0.2 60 0.2];

mdl = fitnlm(X,y,modelfun,beta0);


X = adata(:,2); y = adata(:,3); group = adata(:,4);
% yfit = modelfun([60 0.1],,VFUN)
model = @(PHI,d)PHI(1).*exp(PHI(2)*d);
beta0 = [60 0.1];% fegd(:,:,1) = beta0; fegd(:,:,2) = beta0;
options = statset('nlmefit');
clc
group = categorical(group);
options = statset(options,'MaxIter',350,'TolX',1e-8,'Display','final');
[beta,PSI,stats,B] = nlmefit(X,y,group,[],model,beta0,'Options',options,'RefineBeta0','On');
%%
% data = array2table(num);
% data.Properties.VariableNames = strings;
% data.Group = categorical(data.Group);
% within = [1:8]';
% within = array2table(within);
% within.Properties.VariableNames = {'Day'};
% within.Day = categorical(within.Day);
% ra = repeatedMeasuresAnova(data,within,0.05);
[within,dvn,xlabels] = make_within_table({'Day'},[8]);
dataT = make_between_table({num(num(:,1)==1,2:end);num(num(:,1)==2,2:end)},dvn);
ra = RMA(dataT,within,{0.05,{'bonferroni'}});
ra.ranova
print_for_manuscript(ra)
n = 0;


    mVar = ra.EM.Group_by_Day.Mean;
    semVar = ra.EM.Group_by_Day.Formula_StdErr;
%     combs = ra.MC.Group_by_Day.combs; p = ra.MC.Group_by_Day.p; h = ra.MC.Group_by_Day.p < -0.05;
    xdata = 1:8;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.9 1],'color','w');
    hold on;
    plot(xdata,mVar(1:8),'.','color',tcolors{1},'markersize',mSize,'marker','.','linestyle','-');
    errorbar(xdata,mVar(1:8),semVar(1:8),'.','color',tcolors{1},'markersize',mSize);
    plot(xdata,mVar(9:16),'.','color',tcolors{2},'markersize',mSize,'marker','.','linestyle','-');
    errorbar(xdata,mVar(9:16),semVar(9:16),'.','color',tcolors{2},'markersize',mSize);

    set(gca,'xlim',[0 9],'ylim',[0 70],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = [1:8]; %%xticklabels = {'Day1','Day2','Day3'}; xticklabels = repmat(xticklabels,2,1)';
    set(gca,'xtick',xticks);
%     xtickangle(30);
    changePosition(gca,[0.08 0.13 -0.01 -0.1])
    put_axes_labels(gca,{'Days',[0 0 0]},{{'Average','Latency (s)'},[0 0 0]});
    putLegend(gca,{'C-TG','A-TG',[5 1 4 50]});
    save_pdf(hf,mData.pdf_folder,'Figure_1_mwt_anova_latency.pdf',600);
return;
end

%% 12months probe time
if 1
[num,strings,raw] = xlsread(filename,3,'A1:B24');
[h,p,stat] = ttest2(num(1:11,2),num(12:23,2))

[within,dvn,xlabels] = make_within_table({'Day'},[1]);
dataT = make_between_table({num(1:11,2);num(12:23,2)},dvn);
% ra = RMA(dataT,within,{0.05,{''}});
% ra.ranova
% print_for_manuscript(ra)
ra = anova1(dataT{:,2},dataT{:,1})
% [hc,pc,statc] = ttest(num(1:11,2),25)
% [ha,pa,stata] = ttest(num(12:23,2),25)
n = 0;
[mC,semC] = findMeanAndStandardError(num(1:11,2));
[mA,semA] = findMeanAndStandardError(num(12:23,2));
    mVar = [mC mA];
    semVar = [semC semA];
    combs = [1 2]; 
    xdata = 1:2;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.5,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    plot([0 3],[25 25],':','color','k')
    set(gca,'xlim',[0.25 2.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = [1 2]; xticklabels = {'C-TG','A-TG'}; 
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    for ii = 1:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
%     xtickangle(30);
%     rectangle(gca,'Position',[0.75 30.5 1 3.5],'edgecolor','k','facecolor','k');
%     text(1.85,30.5,'Trials','FontSize',5);
%     rectangle(gca,'Position',[6 30.5 1 3.5],'edgecolor','k');
%     text(7.2,30.5,'Inter-Trials','FontSize',5);
    changePosition(gca,[0.03 0.13 -0.2 -0.1])
    put_axes_labels(gca,{[],[0 0 0]},{{'Probe Time (%)'},[0 0 0]});
    
    save_pdf(hf,mData.pdf_folder,'Figure_1_mwt_probe_time.pdf',600);
return;
end

%% 12months platform proximity
if 1
[num,strings,raw] = xlsread(filename,4,'A1:B24');
[h,p,stat] = ttest2(num(1:11,2),num(12:23,2))
n = 0;

dataT = make_between_table({num(1:11,2);num(12:23,2)},dvn);
% ra = RMA(dataT,within,{0.05,{''}});
% ra.ranova
% print_for_manuscript(ra)
ra = anova1(dataT{:,2},dataT{:,1})

[mC,semC] = findMeanAndStandardError(num(1:11,2));
[mA,semA] = findMeanAndStandardError(num(12:23,2));
    mVar = [mC mA];
    semVar = [semC semA];
    combs = [1 2]; 
    xdata = 1:2;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.001);
    plot([0 3],[25 25],':','color','k')
    set(gca,'xlim',[0.25 2.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = [1 2]; xticklabels = {'C-TG','A-TG'}; 
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    for ii = 1:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
%     xtickangle(30);
%     rectangle(gca,'Position',[0.75 30.5 1 3.5],'edgecolor','k','facecolor','k');
%     text(1.85,30.5,'Trials','FontSize',5);
%     rectangle(gca,'Position',[6 30.5 1 3.5],'edgecolor','k');
%     text(7.2,30.5,'Inter-Trials','FontSize',5);
    changePosition(gca,[0.14 0.13 -0.2 -0.1])
    put_axes_labels(gca,{[],[0 0 0]},{{'Platform','Proximity (m)'},[0 0 0]});
    
    save_pdf(hf,mData.pdf_folder,'Figure_1_mwt_platform_proximity.pdf',600);
return;
end