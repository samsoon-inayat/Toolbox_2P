function overall_population_vectors

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_C = evalin('base','ei10_C'); 
ei_A = evalin('base','ei10_A'); 


selContexts = [1 2 3 4];
rasterNames = {'airD','airD','airD','airD'};

RsC = get_rasters_data(ei_C,selContexts,rasterNames);
mRsC = calc_mean_rasters(RsC,1:10);
RsC = find_responsive_rasters(RsC,1:10);
% view_population_vector(Rs,mRs,300);
[resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC] = get_responsive_fraction(RsC);

RsA = get_rasters_data(ei_A,selContexts,rasterNames);
mRsA = calc_mean_rasters(RsA,1:10);
RsA = find_responsive_rasters(RsA,1:10);
% view_population_vector(Rs,mRs,400);
[resp_fractionA,resp_valsA,OIA,mean_OIA,resp_ORA,resp_OR_fractionA,resp_ANDA,resp_AND_fractionA] = get_responsive_fraction(RsA);

dataT = array2table([[1 1 1 1 1 2 2 2 2 2]' 100*[resp_fractionC;resp_fractionA]]);
dataT.Properties.VariableNames = {'Group','C1','C2','C3','C4'};
within = array2table([1 2 3 4]');
within.Properties.VariableNames = {'Cond'};
within.Cond = categorical(within.Cond);
dataT.Group = categorical(dataT.Group);
ra = repeatedMeasuresAnova(dataT,within,0.05);
% writetable(dataT,fullfile(mData.pdf_folder,'percentage_of_place_cells.xls'));
n = 0;
%%
mVar = ra.est_marginal_means.Mean;
semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
xdata = [1 2 3 4 6 7 8 9]; 

hf = figure(7);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.9 1],'color','w');
hold on;
tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
% tcolors = repmat(tcolors,2,1)';
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);

set(gca,'xlim',[0.25 9.75],'ylim',[0 30],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
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
    put_axes_labels(gca,{[],[0 0 0]},{{'Percentage of','Spatially Tuned Cells'},[0 0 0]});

    save_pdf(hf,mData.pdf_folder,'Perc Of Place Cells_contexts_10',600);
%%
%%
mVar = ra.est_marginal_means_wf.Mean;
semVar = ra.est_marginal_means_wf.Formula_StdErr;
combs = ra.mcs_wf.combs; p = ra.mcs_wf.p; h = ra.mcs_wf.p < 0.05;
xdata = [1 2 3 4]; 

hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 7 1.25 1],'color','w');
hold on;
tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
% tcolors = repmat(tcolors,2,1)';
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);

set(gca,'xlim',[0.25 4.75],'ylim',[0 30],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
xticks = xdata; xticklabels = {'C1','C2','C3','C4'}; xticklabels = repmat(xticklabels,1,2);
set(gca,'xtick',xticks,'xticklabels',xticklabels);

changePosition(gca,[0.07 0.02 -0.2 -0.011])
put_axes_labels(gca,{[],[0 0 0]},{'',[0 0 0]});

save_pdf(hf,mData.pdf_folder,'Perc Of Place Cells_contexts_10_pooled.pdf',600);

%% unique cells
    perc_cells_or_C = resp_OR_fractionC*100;
    perc_cells_or_A = resp_OR_fractionA*100;
    [h,p,ci,stats] = ttest2(perc_cells_or_C,perc_cells_or_A)
    
    mVar = [mean(perc_cells_or_C) mean(perc_cells_or_A)]; semVar = [std(perc_cells_or_C)/sqrt(5) std(perc_cells_or_A)/sqrt(5)];
    combs = [1 2]; %p = ra.mcs.p; h = ra.mcs.p < 0.05;
    xdata = [1:2];
%     xdata = [1 2 3 4];
    tcolors = {'k','r'};
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
%     tcolors ={colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.03,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    for ii = 1:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = xdata(1:end)+0; xticklabels = {'Control','APP'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30);
    changePosition(gca,[0.15 0.05 -0.35 -0.1])
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned Cells','in any Condition (%)'},[0 0 0]});
    
    save_pdf(hf,mData.pdf_folder,sprintf('Place_cells_percentage_overall'),600);

%% unique cells
    perc_cells_or_C = resp_AND_fractionC*100;
    perc_cells_or_A = resp_AND_fractionA*100;
    [h,p,ci,stats] = ttest2(perc_cells_or_C,perc_cells_or_A)
    
    mVar = [mean(perc_cells_or_C) mean(perc_cells_or_A)]; semVar = [std(perc_cells_or_C)/sqrt(5) std(perc_cells_or_A)/sqrt(5)];
    combs = [1 2]; %p = ra.mcs.p; h = ra.mcs.p < 0.05;
    xdata = [1:2];
%     xdata = [1 2 3 4];
    tcolors = {'k','r'};
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
%     tcolors ={colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.03,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    for ii = 1:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = xdata(1:end)+0; xticklabels = {'Control','APP'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30);
    changePosition(gca,[0.15 0.05 -0.35 -0.1])
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned Cells','common Condition (%)'},[0 0 0]});
    
    save_pdf(hf,mData.pdf_folder,sprintf('Place_cells_percentage_common'),600);