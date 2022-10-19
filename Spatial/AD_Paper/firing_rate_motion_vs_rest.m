function firing_rate_motion_vs_rest

% mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
% ei_C = evalin('base','ei10_C'); 
% ei_A = evalin('base','ei10_A'); 

% selContexts = [1 2 3 4];
% rasterNames = {'airD','airD','airD','airD'};

Rs_C = oC.Rs;% get_rasters_data(ei_C,selContexts,rasterNames);
Rs_A = oA.Rs;% get_rasters_data(ei_A,selContexts,rasterNames);
% typeP = {'all','vals'
thr = -1;
%%
fileName = fullfile(mData.pd_folder,sprintf('%s_%s_AD_new',mfilename,typeP));
if 1
    cellListC = cell_list_op(sel_pop_C,[],'or',1);
    cellListA = cell_list_op(sel_pop_A,[],'or',1);
    out_C = get_spike_rate(ei_C,thr,cellListC);
    out_A = get_spike_rate(ei_A,thr,cellListA);
    save(fileName,'out_C','out_A','thr');
else
    temp = load(fileName);
    out_C = temp.out_C;
    out_A = temp.out_A;
    thr = temp.thr;
end
tcolors = {'k','r','k','r'};
n=0;
%%
if 1
    perc_cells_or_C = out_C.m_sp_animal_level_resp_cells;
    perc_cells_or_A = out_A.m_sp_animal_level_resp_cells;
    [h,p,ci,stats] = ttest2(perc_cells_or_C,perc_cells_or_A)
    
    mVar = [mean(perc_cells_or_C) mean(perc_cells_or_A)]; semVar = [std(perc_cells_or_C)/sqrt(3) std(perc_cells_or_A)/sqrt(5)];
    combs = [1 2]; %p = ra.mcs.p; h = ra.mcs.p < 0.05;
    xdata = [1:2];
%     xdata = [1 2 3 4];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
%     tcolors ={colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.03,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 1:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 0.08],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = xdata(1:end)+0; xticklabels = {'C-TG','A-TG'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30);
    changePosition(gca,[0.15 0.05 -0.35 -0.1])
    put_axes_labels(gca,{[],[0 0 0]},{{'Average Firing','Rate (AU)'},[0 0 0]});
    
    save_pdf(hf,mData.pdf_folder,sprintf('Firing_Rate_overall'),600);
return;
end
%%
if 1
    perc_cells_or_C = out_C.m_sp_animal_level_resp_cells;
    perc_cells_or_A = out_A.m_sp_animal_level_resp_cells;
    [h,p,ci,stats] = ttest2(perc_cells_or_C,perc_cells_or_A)
    
    mVar = [mean(perc_cells_or_C) mean(perc_cells_or_A)]; semVar = [std(perc_cells_or_C)/sqrt(3) std(perc_cells_or_A)/sqrt(5)];
    combs = [1 2]; %p = ra.mcs.p; h = ra.mcs.p < 0.05;
    xdata = [1:2];
%     xdata = [1 2 3 4];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
%     tcolors ={colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.03,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 1:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 0.08],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = xdata(1:end)+0; xticklabels = {'C-TG','A-TG'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30);
    changePosition(gca,[0.15 0.05 -0.25 -0.1])
    put_axes_labels(gca,{[],[0 0 0]},{{'Average Firing','Rate (AU)'},[0 0 0]});
    
    save_pdf(hf,mData.pdf_folder,sprintf('Firing_Rate_overall'),600);
return;
end
% %%
% if 1
%     perc_cells_or_C = out_C.percent_motion;
%     perc_cells_or_A = out_A.percent_motion;
%     [h,p,ci,stats] = ttest2(perc_cells_or_C,perc_cells_or_A)
%     
%     mVar = [mean(perc_cells_or_C) mean(perc_cells_or_A)]; semVar = [std(perc_cells_or_C)/sqrt(3) std(perc_cells_or_A)/sqrt(5)];
%     combs = [1 2]; %p = ra.mcs.p; h = ra.mcs.p < 0.05;
%     xdata = [1:2];
% %     xdata = [1 2 3 4];
%     colors = mData.colors;
%     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
%     hold on;
%     tcolors ={colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
%     [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
%         'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
%         'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
%     for ii = 5:length(hbs)
%         set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
%     end
%     % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
%     set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
%     xticks = xdata(1:end)+0; xticklabels = {'RSEG','PSEG'};
%     set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30);
%     changePosition(gca,[0.2 0 -0.4 -0.09]);
%     put_axes_labels(gca,{[],[0 0 0]},{{'Percent','Motion (%)'},[0 0 0]});
%     
%     save_pdf(hf,mData.pdf_folder,sprintf('Firing_Rate_Percent_Motion'),600);
% return;
% end
%%
if 1
    tcolors = {'k','k','k','k'};

    data_C = [out_C.m_sp_animal_level_rest' out_C.m_sp_animal_level_motion'];
    data_A = [out_A.m_sp_animal_level_rest' out_A.m_sp_animal_level_motion'];
    data = [data_C;data_A];
    dataT = [table([ones(size(data_C,1),1);(2*ones(size(data_A,1),1))]) array2table(data)];
    dataT.Properties.VariableNames = {'Group','Rest','Motion'};
    dataT.Group = categorical(dataT.Group);

    within = table([1;2]);
    within.Properties.VariableNames = {'Type'};
    within.Type = categorical(within.Type);
    ra = repeatedMeasuresAnova(dataT,within);
%     writetable(dataT,fullfile(mData.pdf_folder,'FR_motion_rest_Table.xls'));
    percentChangeC = 100*(data_C(:,2)-data_C(:,1))./data_C(:,1);
    percentChangeA = 100*(data_A(:,2)-data_A(:,1))./data_A(:,1);
    [h,p,ci,stats] = ttest2(percentChangeC,percentChangeA)
    [mpcC,sempcC] = findMeanAndStandardError(percentChangeC);
    [mpcA,sempcA] = findMeanAndStandardError(percentChangeA);

    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = 2*ra.mcs.p; h = ra.mcs.p < 0.05;
    combs = circshift(combs,1); p = circshift(p,1); h = circshift(h,1);
     xdata = [1 2 4 5]; 
    colors = mData.colors;
    zMI_Th = NaN;
    if isnan(zMI_Th)
       hf = figure(7);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    else
        hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 9 1.5 1],'color','w');
    end
    hold on;
%     tcolors ={colors{1};colors{2};colors{1};colors{2};colors{3};colors{4}};
  h(end) = 0;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.005,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.0001);
    for ii = 1:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    set(hbs(1),'facecolor',[0.5 0.5 0.5],'edgecolor','k');
    set(hbs(3),'facecolor',[0.5 0.5 0.5],'edgecolor','k');
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 0.08],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = [1.5 4.5]; xticklabels = {'C-TG','A-TG'};
%     xtickangle(20)
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.17 0.05 -0.08 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Average Firing', 'Rate (Hz)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Firing_Rate_Motion_vs_Rest'),600);
return;
end


%%
if 1
    tcolors = {'k','r'};
   distD(:,1) = out_C.allVals_an';
   distD(:,2) = out_A.allVals_an';
   [distDo,allVals,allValsG] = plotDistributions(distD);
   minBin = min(allVals);
   maxBin = max(allVals);
   [h,p,ks2stat] = kstest2(allValsG{1},allValsG{2});
   [h,p,cd,ks2stat] = ttest2(allValsG{1},allValsG{2});
   %%
   incr = 0.001; %maxBin =
   hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.85 1],'color','w');
   hold on;
%    [ha,hb,hca] = plotDistributions(distD,'colors',tcolors,'maxY',maxBin,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
   [ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
   set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
   changePosition(gca,[0.09 0.13 -0.05 -0.13]);
    put_axes_labels(gca,{'Average Firing Rate (AU)',[0 0 0]},{{'Percentage','of Neurons'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_firing_rate'),600);
   %%
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    for ii = 1:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    set(hbs(1),'facecolor',[0.5 0.5 0.5],'edgecolor','k');
    set(hbs(3),'facecolor',[0.5 0.5 0.5],'edgecolor','k');
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = [1.5 4.5]; xticklabels = {'Control','APP'};
%     xtickangle(20)
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.14 0.05 -0.05 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Average Firing', 'Rate (Hz)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Firing_Rate_Motion_vs_Rest'),600);
return;
end

