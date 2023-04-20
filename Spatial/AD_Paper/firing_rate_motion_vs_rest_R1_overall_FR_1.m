function firing_rate_motion_vs_rest

% mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
% ei_C = evalin('base','ei10_C'); 
% ei_A = evalin('base','ei10_A'); 

% selContexts = [1 2 3 4];
% rasterNames = {'airD','airD','airD','airD'};

Rs_C = oC.Rs;% get_rasters_data(ei_C,selContexts,rasterNames);
Rs_A = oA.Rs;% get_rasters_data(ei_A,selContexts,rasterNames);
% typeP = {'all','vals'
typeP  = 'all';
thr = -1;
%%
fileName = fullfile(mData.pd_folder,sprintf('%s_%s_AD_new',mfilename,typeP));
if 1
    cellListC = cell_list_op(sel_pop_C,[],'or',1);
    cellListA = cell_list_op(sel_pop_A,[],'or',1);
    out_C_Ca = get_spike_rate_Ca(ei_C,thr,cellListC);
    out_A_Ca = get_spike_rate_Ca(ei_A,thr,cellListA);
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

%% average firing rate

    perc_cells_or_C = out_C.m_sp_animal_level_resp_cells;
    perc_cells_or_A = out_A.m_sp_animal_level_resp_cells;
    [h,p,ci,stats] = ttest2(perc_cells_or_C,perc_cells_or_A);
    esiz = computeCohen_d(perc_cells_or_C,perc_cells_or_A)
    
    mVar = [mean(perc_cells_or_C) mean(perc_cells_or_A)]; semVar = [std(perc_cells_or_C)/sqrt(5) std(perc_cells_or_A)/sqrt(5)];
    combs = [1 2]; %p = ra.mcs.p; h = ra.mcs.p < 0.05;
    xdata = [1:2];
%     xdata = [1 2 3 4];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
%     tcolors ={colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.02,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.051,...
        'raw_data',[perc_cells_or_C' perc_cells_or_A'],'dots','o','dot_size',1);
    for ii = 1:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 0.08],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = xdata(1:end)+0; xticklabels = {'C-TG','A-TG'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.25 0.05 -0.45 -0.1])
    put_axes_labels(gca,{[],[0 0 0]},{{'Avg. Firing','Rate (A.U.)'},[0 0 0]});
    format_axes_b(gca);
    
    save_pdf(hf,mData.pdf_folder,sprintf('Firing_Rate_overall'),600);
    
%% rest vs motion FR
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
%     ra = RMA(dataT,within,{0.05,{'bonferroni'}});
%     print_for_manuscript(ra)
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
    changePosition(gca,[0.17 0.05 -0.2 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Average Firing', 'Rate (A.U.)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Firing_Rate_Motion_vs_Rest'),600);
return;
end


%% rest vs motion FR
    tcolors = {'k','k','k','k'};
    rest_tr_C = out_C.m_sp_animal_level_rest; motion_tr_C = out_C.m_sp_animal_level_motion;
    rest_tr_A = out_A.m_sp_animal_level_rest; motion_tr_A = out_A.m_sp_animal_level_motion;
    
    varC = [rest_tr_C' motion_tr_C'];
    varA = [rest_tr_A' motion_tr_A'];
    [within,dvn,xlabels] = make_within_table({'RM'},[2]);
    dataT = make_between_table({varC;varA},dvn);
    ra = RMA(dataT,within,{0.05,{'bonferroni'}});
    ra.ranova
    print_for_manuscript(ra)
    
    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[10 3 1.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.31 0.40],'widthHeightAdjustment',[10 -450]);
    MY = 0.1; ysp = 0.01; mY = 0; titletxt = ''; ylabeltxt = {'Avgerage Firing','Rate (A.U.)'};
    stp = 0.485*magfac; widths = ([1 1.3 1.3 1.3 1.3 0.5 0.5 0.5]-0)*magfac; gap = 0.16*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{''});
    colors = mData.dcolors(7:end);
    tcolors = repmat({colors{1};colors{2};colors{1};colors{2}},4);

    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_RM','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([2 2],[1 1.5]); 
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
    xticklabels = {'Rest','Motion'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(25);
    make_bars_hollow(hbs(3:end));
    put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
    ht = set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'C-TG','A-TG'},{[-0.03 -0]}); 
    for ii = 1:length(ht) 
      set(ht(ii),'FontWeight','Bold');
    end
    format_axes_b(gca);
    save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);



%% rest vs motion transients
    tcolors = {'k','k','k','k'};
    rest_tr_C = exec_fun_on_cell_mat(out_C.m_tr_animal_rest,'mean'); motion_tr_C = exec_fun_on_cell_mat(out_C.m_tr_animal_motion,'mean');
    rest_tr_A = exec_fun_on_cell_mat(out_A.m_tr_animal_rest,'mean'); motion_tr_A = exec_fun_on_cell_mat(out_A.m_tr_animal_motion,'mean');
    
    varC = [rest_tr_C' motion_tr_C'];
    varA = [rest_tr_A' motion_tr_A'];
    [within,dvn,xlabels] = make_within_table({'RM'},[2]);
    dataT = make_between_table({varC;varA},dvn);
    ra = RMA(dataT,within,{0.05,{'bonferroni'}});
    ra.ranova
    print_for_manuscript(ra)
    
    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[13 5 1.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.31 0.40],'widthHeightAdjustment',[10 -450]);
    MY = 700; ysp = 100.5; mY = 0; titletxt = ''; ylabeltxt = {'Avg. # of','Trans./min'};
    stp = 0.45*magfac; widths = ([1 1.3 1.3 1.3 1.3 0.5 0.5 0.5]-0)*magfac; gap = 0.16*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{''});
    colors = mData.dcolors(7:end);
    tcolors = repmat({colors{1};colors{2};colors{1};colors{2}},4);

    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_RM','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([2 2],[1 1.5]); 
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
    xticklabels = {'Rest','Motion'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(25);
    make_bars_hollow(hbs(3:end));
    put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
    ht = set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'C-TG','A-TG'},{[-0.03 -0]}); 
    for ii = 1:length(ht) 
      set(ht(ii),'FontWeight','Bold');
    end
    format_axes_b(gca);
    save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);


%%
if 1
    tcolors = {'k','r'};
   distD(:,1) = out_C.allVals_an';
   distD(:,2) = out_A.allVals_an';
   [distDo,allVals,allValsG] = plotDistributions(distD);
   minBin = min(allVals);
   maxBin = max(allVals);
   [h,p,ks2stat] = kstest2(allValsG{1},allValsG{2}); n1 = length(allValsG{1}); n2 = length(allValsG{2});
   n      =  n1 * n2 /(n1 + n2);
%    [h,p,cd,ks2stat] = ttest2(allValsG{1},allValsG{2});
   %%
   incr = 0.001; %maxBin =
   hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
   hold on;
%    [ha,hb,hca] = plotDistributions(distD,'colors',tcolors,'maxY',maxBin,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
   [ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
   set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
   changePosition(gca,[0.129 0.15 -0.09 -0.13]);
   ylim([0 100]); xlim([minBin maxBin]);
    put_axes_labels(gca,{'Avg. Firing Rate (A.U.)',[0 0 0]},{{'Cells (%)'},[0 0 0]});
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

%%
ppm_C = find_transients_per_minute(ei_C);
ppm_A = find_transients_per_minute(ei_A);
%%
    tcolors = {'k','r'};
   distD(:,1) = ppm_C';
   distD(:,2) = ppm_A';
   [distDo,allVals,allValsG] = plotDistributions(distD);
   minBin = min(allVals);
   maxBin = max(allVals);
   [h,p,ks2stat] = kstest2(allValsG{1},allValsG{2});
%    [h,p,cd,ks2stat] = ttest2(allValsG{1},allValsG{2});
   %%
   incr = 0.1; %maxBin =
   hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
   hold on;
%    [ha,hb,hca] = plotDistributions(distD,'colors',tcolors,'maxY',maxBin,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
   [ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
   set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
   changePosition(gca,[0.129 0.15 -0.09 -0.13]);
   ylim([0 100]); xlim([minBin maxBin]);
    put_axes_labels(gca,{{'Avg. # of Trans. per min'},[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_firing_rate'),600);
    %%
    perc_cells_or_C = exec_fun_on_cell_mat(ppm_C,'mean');
    perc_cells_or_A = exec_fun_on_cell_mat(ppm_A,'mean');
    [h,p,ci,stats] = ttest2(perc_cells_or_C,perc_cells_or_A)
    
    mVar = [mean(perc_cells_or_C) mean(perc_cells_or_A)]; semVar = [std(perc_cells_or_C)/sqrt(5) std(perc_cells_or_A)/sqrt(5)];
    combs = [1 2]; %p = ra.mcs.p; h = ra.mcs.p < 0.05;
    xdata = [1:2];
%     xdata = [1 2 3 4];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
%     tcolors ={colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.02,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.051,...
        'raw_data',[perc_cells_or_C' perc_cells_or_A'],'dots','o','dot_size',1);
    for ii = 1:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 20],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = xdata(1:end)+0; xticklabels = {'C-TG','A-TG'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.25 0.05 -0.45 -0.1])
    put_axes_labels(gca,{[],[0 0 0]},{{'Avg. # of','Trans./min'},[0 0 0]});
    format_axes_b(gca);
    
    save_pdf(hf,mData.pdf_folder,sprintf('Firing_Rate_overall'),600);