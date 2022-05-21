function light_figure
%%
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei = evalin('base','ei'); 
    selContexts = [1 2 3 3 3 3 4 4 4 4 5 5 5 5 6 7 4 0 0];
    rasterNames = {'light22T','air55T','airD','airT','airID','airIT','airD','airT','airID','airIT','airD','airT','airID','airIT','light22T','air55T','light22T','motionOnsets','motionOffsets'};
    o = get_data(ei,selContexts,rasterNames);
    selContexts = [1 4 6];
    rasterNames = {'light22T','light22T','light22T'};
    oLight = get_data(ei,selContexts,rasterNames);
    break;
end
n = 0;
%% Percentage of Responsive Cells overall for all conditions
while 1
    o.props = get_props_Rs(o.Rs);
%     resp_vals = cell_list_op(o.resp.vals,o.resp_FR_ei,'and');
    resp_test = (cell_list_op(o.props.good_zMI,o.resp.vals,'and'));
    resp = cell_list_op(o.props.good_zMI,o.props.good_FR,'and');
%     resp = o.resp.vals;
    resp_OR = cell_list_op(resp,[],'or');
    per_resp = find_percent(resp);
    [m_per_resp_OR,sem_per_resp_OR] = findMeanAndStandardError(find_percent(resp_OR));
    [within,dvn,xlabels] = make_within_table({'Cond'},[19]);
    dataT = make_between_table({per_resp},dvn);
    ra = RMA(dataT,within,0.05);

    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','bonferroni'},[1 1 1]);
    hf = get_figure(5,[8 7 3.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'C1','C2','C3-DT','C3-TT','C3-DIT','C3-TIT','C4-DT','C4-TT','C4-DIT','C4-TIT','C3''-DT','C3''-TT','C3''-DIT','C3''-TIT','C1''','C2''','C4L','M-On','M-Off'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.15 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('all_responsive_cells.pdf'),600);
    ra.ranova
    ra.mauchly
    [m_per_resp_OR(1),sem_per_resp_OR(1)]
    break;
end
%%
% an = 5; cn = 2;
% plotRasters_simplest(Rs{an,cn})
%%
an = 3; cn = 1;
% plotRasters_simplest(Rs{an,cn})
% find(resp_valsC{an}(:,cn));
ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
    'spaceRowsCols',[0.15 0.03],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
    [-50 -375]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[10 4 3.25 1]);
ff = sample_rasters(Rs{an,cn},[191 11 96 41],ff);
save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersD'),600);
%% population vector and correlation single animal
an = 3;
ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
    'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.07 0.1],'widthHeightAdjustment',...
    [0.01 -60]);
set(gcf,'color','w');
set(gcf,'Position',[5 5 3.25 2]);
ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
save_pdf(ff.hf,mData.pdf_folder,sprintf('air_population_vector_corrD.pdf'),600);

%% average correlation of all animals
ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 3],...
    'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.075 0.2],'widthHeightAdjustment',...
    [0.01 -220]);
set(gcf,'color','w');
set(gcf,'Position',[5 5 3.25 1]);
ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
save_pdf(ff.hf,mData.pdf_folder,sprintf('air_average_population_vector_corrD.pdf'),600);

%%
dataT = array2table(resp_fractionC*100);
dataT.Properties.VariableNames = {'A1','A2','A3'};
within = array2table([1 2 3]');
within.Properties.VariableNames = {'Cond'};
within.Cond = categorical(within.Cond);
ra = repeatedMeasuresAnova(dataT,within,0.05);

%%
mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
    % row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
    % row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

    xdata = [1 2 3]; 
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
%     tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
    tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.05);
    for ii = 5:length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
%     maxY = maxY + 5 + 2;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30)
    changePosition(gca,[0.2 0.03 -0.3 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('fraction_air_responsiveD'),600);
%% overlap of cells
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[7 7 1.25 1],'color','w');
    hold on;
    imagesc(mean_OIC);
    axis equal
    colorbar;
    xlim([0.5 3.5]);
    ylim([0.5 3.5]);
    changePosition(gca,[0.1 0.03 0 0]);
    xticklabels = {'C3','C4','C3'''};
    set(gca,'XTick',[1 2 3],'XTickLabels',xticklabels,'YTick',[1 2 3],'YTickLabels',(xticklabels));
    set(gca,'Ydir','reverse','linewidth',0.5,'FontSize',6,'FontWeight','Bold');
    save_pdf(hf,mData.pdf_folder,sprintf('air_overlap_imageD'),600);
%% overlap RM bar graph
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    for ii = 1:5
        C12(ii) = OIC{ii}(1,2);
        C13(ii) = OIC{ii}(1,3);
        C23(ii) = OIC{ii}(2,3);
    end
    dataT = array2table([C12' C13' C23']);
    dataT.Properties.VariableNames = {'A1','A2','A3'};
    within = array2table([1 2 3]');
    within.Properties.VariableNames = {'Cond'};
    within.Cond = categorical(within.Cond);
    ra = repeatedMeasuresAnova(dataT,within,0.05);

%%
mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
    % row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
    % row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

    xdata = [1 2 3]; 
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
%     tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
    tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.05);
    for ii = 5:length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
%     maxY = maxY + 5 + 2;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'C3-C4','C3-C3''','C4-C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30)
    changePosition(gca,[0.2 0.03 -0.3 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Overlap Index'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('air_overlap_statsD'),600);
    
    
%%
    [mVar,semVar] = findMeanAndStandardError(resp_OR_fractionC*100);
    combs = []; p = 1; h = p < 0.05;
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
    % row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
    % row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

    xdata = [1]; 
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
%     tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
    tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);

    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {''''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30)
    changePosition(gca,[0.3 0.03 -0.5 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned Cells','in any Condition (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('air_Unique_place_Cells'),600);