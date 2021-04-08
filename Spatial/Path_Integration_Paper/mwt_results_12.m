function mwt_results

mData = evalin('base','mData');
colors = mData.colors;

filename = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\PDFs\APP Figure 1 MWT Data.xlsx';

marker_styles = {'s','d','^','s','d','^','^','*'};
line_styles = {'-',':','-.','-',':','-.'};
mSize = 5;
%% 12months speed
if 1
[num,strings,raw] = xlsread(filename,1,'A1:I24');
data = array2table(num);
data.Properties.VariableNames = strings;
data.Group = categorical(data.Group);
within = [1:8]';
within = array2table(within);
within.Properties.VariableNames = {'Day'};
within.Day = categorical(within.Day);
ra = repeatedMeasuresAnova(data,within,0.05);
n = 0;


    mVar = ra.est_marginal_means.Mean;
    semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < -0.05;
    xdata = 1:8;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.9 1],'color','w');
    hold on;
    tcolors = {'b','r'};
    plot(xdata,mVar(1:8),'.','color',tcolors{1},'markersize',mSize,'marker','.','linestyle','-');
    errorbar(xdata,mVar(1:8),semVar(1:8),'.','color',tcolors{1},'markersize',mSize);
    plot(xdata,mVar(9:16),'.','color',tcolors{2},'markersize',mSize,'marker','.','linestyle',':');
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
if 0
[num,strings,raw] = xlsread(filename,2,'A1:I24');
data = array2table(num);
data.Properties.VariableNames = strings;
data.Group = categorical(data.Group);
within = [1:8]';
within = array2table(within);
within.Properties.VariableNames = {'Day'};
within.Day = categorical(within.Day);
ra = repeatedMeasuresAnova(data,within,0.05);
n = 0;


    mVar = ra.est_marginal_means.Mean;
    semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < -0.05;
    xdata = 1:8;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.9 1],'color','w');
    hold on;
    tcolors = {'b','r'};
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
    save_pdf(hf,mData.pdf_folder,'Figure_1_mwt_anova_latency.pdf',600);
return;
end

%% 12months probe time
if 1
[num,strings,raw] = xlsread(filename,3,'A1:B24');
[h,p,stat] = ttest2(num(1:8,2),num(9:16,2))
n = 0;


    mVar = ra.est_marginal_means.Mean;
    semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < -0.05;
    xdata = 1:8;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.9 1],'color','w');
    hold on;
    tcolors = {'b','r'};
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

    save_pdf(hf,mData.pdf_folder,'Figure_1_mwt_anova_latency.pdf',600);
return;
end

