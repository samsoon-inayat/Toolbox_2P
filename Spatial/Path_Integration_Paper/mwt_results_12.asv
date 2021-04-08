function mwt_results

mData = evalin('base','mData');
colors = mData.colors;

filename = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\PDFs\APP Figure 1 MWT Data.xlsx';

marker_styles = {'s','d','^','s','d','^','^','*'};
line_styles = {'-',':','-.','-',':','-.'};
mSize = 8;
%% 12months speed



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
%%
if 1
    mVar = ra.est_marginal_means.Mean;
    semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < -0.05;
    xdata = 1:8;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.9 1],'color','w');
    hold on;
    ii = 1;
    plot(xdata,mVar(1:8),'.','color',colors{ii},'markersize',mSize,'marker',marker_styles{ii});
    plot(xdata,mVar(9:16),'linewidth',1.5,'color',colors{ii+1},'linestyle',line_styles{ii});
    errorbar(xs1,allTrials(:,ii),allSEM(:,ii),'.','color',thisCols{ii},'markersize',mSize);

    set(gca,'xlim',[0.25 14.75],'ylim',[0 46.1434],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = [1.5 3.5 5.5 9.5 11.5 13.5]; xticklabels = {'Day1','Day2','Day3'}; xticklabels = repmat(xticklabels,2,1)';
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    for ii = 2:2:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
%     rectangle(gca,'Position',[0.75 30.5 1 3.5],'edgecolor','k','facecolor','k');
%     text(1.85,30.5,'Trials','FontSize',5);
%     rectangle(gca,'Position',[6 30.5 1 3.5],'edgecolor','k');
%     text(7.2,30.5,'Inter-Trials','FontSize',5);
    changePosition(gca,[0.07 0.02 -0.01 -0.011])
    put_axes_labels(gca,{[],[0 0 0]},{{'Average','Speed (cm/s)'},[0 0 0]});

    save_pdf(hf,mData.pdf_folder,'Figure_1_behavior_anova_speed.pdf',600);
return;
end
