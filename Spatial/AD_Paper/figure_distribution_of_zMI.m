function figure_distribution_of_zMI

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_C = evalin('base','ei10_C1'); 
ei_A = evalin('base','ei10_A1'); 

selContexts = [1 2 3 4];
rasterNames = {'airD','airD','airD','airD'};

RsC = get_rasters_data(ei_C,selContexts,rasterNames);
mRsC = calc_mean_rasters(RsC,1:10);
RsC = find_responsive_rasters(RsC,1:10);
% view_population_vector(Rs,mRs,300);
[resp_fractionC,resp_valsC,OIC,mean_OIC] = get_responsive_fraction(RsC);

RsA = get_rasters_data(ei_A,selContexts,rasterNames);
mRsA = calc_mean_rasters(RsA,1:10);
RsA = find_responsive_rasters(RsA,1:10);
% view_population_vector(Rs,mRs,400);
[resp_fractionA,resp_valsA,OIA,mean_OIA] = get_responsive_fraction(RsA);

n = 0;
%%
for rr = 1:size(RsC,1)
    for cc = 1:size(RsC,2)
        R = RsC{rr,cc};
        zMIsC{rr,cc} = R.info_metrics.ShannonMI_Zsh';
        R = RsA{rr,cc};
        zMIsA{rr,cc} = R.info_metrics.ShannonMI_Zsh';
    end
end
n = 0;
%%
if 1
CN = 3;
tcolors = {'k','r'};
distD(:,1) = zMIsC(:,CN);
distD(:,2) = zMIsA(:,CN);
[distDo,allVals,allValsG] = plotDistributions(distD);
minBin = min(allVals);
maxBin = max(allVals);

tcolors = {'k','r'};
incr = 0.001; %maxBin =
hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
hold on;
%    [ha,hb,hca] = plotDistributions(distD,'colors',tcolors,'maxY',maxBin,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
[ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
plot([1.65 1.65],[0 100],'--k');
xlim([-5 30]);
set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
changePosition(gca,[0.15 0.13 -0.2 -0.13]);
if CN > 1
    put_axes_labels(gca,{'Mutual Information (z-score)',[0 0 0]},{{'Percentage','of Neurons'},[0 0 0]});
else
    put_axes_labels(gca,{'Mutual Information (z-score)',[0 0 0]},{{'Percentage','of Neurons'},[1 0 0]});
end
save_pdf(hf,mData.pdf_folder,sprintf('Distribution_zMI_%d',CN),600);
end
%%
n = 0;
for rr = 1:size(RsC,1)
    for cc = 1:size(RsC,2)
        R = RsC{rr,cc};
        mzMIsC(rr,cc) = nanmean(zMIsC{rr,cc});
        R = RsA{rr,cc};
        mzMIsA(rr,cc) = nanmean(zMIsA{rr,cc});
    end
end

%%
dataT = array2table([[1;1;1;1;1;2;2;2;2;2] [mzMIsC;mzMIsA]]);
dataT.Properties.VariableNames = {'Group','C1','C2','C3','C4'};
within = array2table([1 2 3 4]');
within.Properties.VariableNames = {'Cond'};
within.Cond = categorical(within.Cond);
dataT.Group = categorical(dataT.Group);
ra = repeatedMeasuresAnova(dataT,within,0.05);
% writetable(dataT,fullfile(mData.pdf_folder,'zMI_all_cells.xls'));
%%
mVar = ra.est_marginal_means.Mean;
semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
xdata = [1 2 3 4 6 7 8 9]; 

hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.9 1],'color','w');
hold on;
tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
% tcolors = repmat(tcolors,2,1)';
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);

set(gca,'xlim',[0.25 9.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
xticks = xdata; xticklabels = {'C1','C2','C3','C4'}; xticklabels = repmat(xticklabels,1,2);
set(gca,'xtick',xticks,'xticklabels',xticklabels);

for ii = 5:length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
end
%     rectangle(gca,'Position',[0.75 30.5 1 3.5],'edgecolor','k','facecolor','k');
%     text(1.85,30.5,'Trials','FontSize',5);
%     rectangle(gca,'Position',[6 30.5 1 3.5],'edgecolor','k');
%     text(7.2,30.5,'Inter-Trials','FontSize',5);
changePosition(gca,[0.07 0.02 -0.01 -0.011])
put_axes_labels(gca,{[],[0 0 0]},{{'Mutual Information','(z-score)'},[0 0 0]});

save_pdf(hf,mData.pdf_folder,'zMI_bar_graph_all_cells.pdf',600);


%%
mVar = ra.est_marginal_means_wf.Mean;
semVar = ra.est_marginal_means_wf.Formula_StdErr;
combs = ra.mcs_wf.combs; p = ra.mcs_wf.p; h = ra.mcs_wf.p < 0.05;
xdata = [1 2 3 4]; 

hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 7 1.25 1],'color','w');
hold on;
tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
% tcolors = repmat(tcolors,2,1)';
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);

set(gca,'xlim',[0.25 4.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
xticks = xdata; xticklabels = {'C1','C2','C3','C4'}; xticklabels = repmat(xticklabels,1,2);
set(gca,'xtick',xticks,'xticklabels',xticklabels);

changePosition(gca,[0.07 0.02 -0.2 -0.011])
put_axes_labels(gca,{[],[0 0 0]},{'',[0 0 0]});

save_pdf(hf,mData.pdf_folder,'zMI_bar_graph_all_cells_pooled.pdf',600);