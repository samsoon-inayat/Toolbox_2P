function figure_place_remapping_AD(fn,allRs,ccs)

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_C = evalin('base','ei10_C1'); 
ei_A = evalin('base','ei10_A1'); 

selContexts = [1 2 3 4];
rasterNames = {'airD','airD','airD','airD'};

RsC = get_rasters_data(ei_C,selContexts,rasterNames);
mRsC = calc_mean_rasters(RsC,1:10);
RsC = find_responsive_rasters(RsC,1:10);
[CR_C,aCR_C] = find_population_vector_corr(RsC,mRsC,1,0);
[resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC] = get_responsive_fraction(RsC);

RsA = get_rasters_data(ei_A,selContexts,rasterNames);
mRsA = calc_mean_rasters(RsA,1:10);
RsA = find_responsive_rasters(RsA,1:10);
[CR_A,aCR_A] = find_population_vector_corr(RsA,mRsA,1,0);
[resp_fractionA,resp_valsA,OIA,mean_OIA,resp_ORA,resp_OR_fractionA,resp_ANDA,resp_AND_fractionA] = get_responsive_fraction(RsA);

respC_pop = get_cell_list(resp_valsC,[]);
respA_pop = get_cell_list(resp_valsA,[]);

respC = get_cell_list(resp_valsC,[1;2;3;4]);
respA = get_cell_list(resp_valsA,[1;2;3;4]);

% [out.allP_an,out.allC_an,out.avg_C_conds,out.mean_rasters_T,out.all_corr_an,out.all_corr_cell_an] = get_pop_vector_corr(out,conditionsAndRasterTypes,min(out.sz(:)),cellSel_C);
out_C = find_population_vector_corr_remap(RsC,mRsC,respC);
out_A = find_population_vector_corr_remap(RsA,mRsA,respA);

out_C_pop = find_population_vector_corr_remap(RsC,mRsC,respC_pop);
out_A_pop = find_population_vector_corr_remap(RsA,mRsA,respA_pop);

selC = out_C_pop;
selA = out_A_pop;
n = 0;
%% Correlations across conditions
typeCorr = {'Spatial Correlation','Pop. Vec. Correlation','\Delta FR Score'};
FF = {'SP','PV','RR'};
ysp = [0.05 0.05 0.1];
for ci = 1;
if 1
    [within,dvn,xlabels] = make_within_table({'Cond'},3);
    switch ci
        case 1
            var_C = arrayfun(@(x) mean(x{1}),selC.adj_SP_corr_diag);var_A = arrayfun(@(x) mean(x{1}),selA.adj_SP_corr_diag);
        case 2
            var_C = arrayfun(@(x) mean(x{1}),selC.adj_PV_corr_diag);var_A = arrayfun(@(x) mean(x{1}),selA.adj_PV_corr_diag);
        case 3
            var_C = arrayfun(@(x) nanmean(x{1}),selC.adj_RR_SP);var_A = arrayfun(@(x) nanmean(x{1}),selA.adj_RR_SP);
    end
    dataT = make_between_table({var_C;var_A},dvn);
    ra = repeatedMeasuresAnova(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
    colors = mData.colors;
    hf = get_figure(5,[3 7 1.5 1]);
    tcolors = mData.colors(1:3); tcolors = repmat(tcolors,1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp(ci),'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C12','C23','C34','C12','C23','C34'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    for ii = 4:6
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    xtickangle(45);
    changePosition(gca,[0.1 0.03 -0.1 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{typeCorr{ci},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_correlation',FF{ci}),600);
end
end
%%
varC = selC.adj_SP_corr_diag;
varA = selA.adj_SP_corr_diag;
if 1
    CN = 1;
    tcolors = {'k','r'};
    distD(:,1) = varC(:,CN);
    distD(:,2) = varA(:,CN);
    [distDo,allVals,allValsG] = plotDistributions(distD);
    minBin = min(allVals);
    maxBin = max(allVals);

    tcolors = {'k','r'};
    incr = 0.001; %maxBin =
    hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    %    [ha,hb,hca] = plotDistributions(distD,'colors',tcolors,'maxY',maxBin,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    [ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
    % plot([1.65 1.65],[0 100],'--k');
    % xlim([-5 30]);
    set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    changePosition(gca,[0.15 0.13 -0.2 -0.13]);
    put_axes_labels(gca,{'Spatial Correlation',[0 0 0]},{{'Percentage','of Neurons'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_SP_Corr_%d',CN),600);
return;
end


%%
varC = selC.adj_PV_corr_diag;
varA = selA.adj_PV_corr_diag;
if 1
    CN = 1;
    tcolors = {'k','r'};
    distD(:,1) = varC(:,CN);
    distD(:,2) = varA(:,CN);
    [distDo,allVals,allValsG] = plotDistributions(distD);
    minBin = min(allVals);
    maxBin = max(allVals);

    tcolors = {'k','r'};
    incr = 0.001; %maxBin =
    hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    %    [ha,hb,hca] = plotDistributions(distD,'colors',tcolors,'maxY',maxBin,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    [ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
    % plot([1.65 1.65],[0 100],'--k');
    % xlim([-5 30]);
    set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    changePosition(gca,[0.15 0.13 -0.2 -0.13]);
    put_axes_labels(gca,{'PV Correlation',[0 0 0]},{{'Percentage','of Neurons'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_PV_Corr_%d',CN),600);
return;
end

%%
varC = selC.adj_RR_SP;
varA = selA.adj_RR_SP;
if 1
    CN = 3;
    tcolors = {'k','r'};
    distD(:,1) = varC(:,CN);
    distD(:,2) = varA(:,CN);
    [distDo,allVals,allValsG] = plotDistributions(distD);
    minBin = min(allVals);
    maxBin = max(allVals);

    tcolors = {'k','r'};
    incr = 0.001; %maxBin =
    hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    %    [ha,hb,hca] = plotDistributions(distD,'colors',tcolors,'maxY',maxBin,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    [ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
    % plot([1.65 1.65],[0 100],'--k');
    % xlim([-5 30]);
    set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    changePosition(gca,[0.15 0.13 -0.2 -0.13]);
    put_axes_labels(gca,{'Spatial Correlation',[0 0 0]},{{'Percentage','of Neurons'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_RR_Corr_%d',CN),600);
return;
end

%% pooled
if 1

    mVar = ra.est_marginal_means_wf.Mean;semVar = ra.est_marginal_means_wf.Formula_StdErr;
    combs = ra.mcs_wf.combs; p = ra.mcs_wf.p; h = ra.mcs_wf.p < 0.05;
    xdata = [1 2 3];
    colors = mData.colors;
    hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[3 7 1.5 1],'color','w');
    hold on;
    tcolors ={colors{1};colors{2};colors{3};colors{1};colors{2};colors{3};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 4:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C12','C23','C34','C12','C23','C34'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.06 0.03 -0.3 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Correlation'},[0 0 0]});

    save_pdf(hf,mData.pdf_folder,'cell corr remap pooled',600);
end

%% average correlation of all animals
if 1
    mean_corr_popV_C = selC.mean_PV_corr;
    xs_C = selC.xs;
    ff = makeFigureRowsCols(106,[1 0.5 4 0.5],'RowsCols',[4 4],...
        'spaceRowsCols',[0.1 0.1],'rightUpShifts',[0.13 0.2],'widthHeightAdjustment',...
        [-150 -150]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 3 4 4]);
    ff = show_remapping_corr_plots(mData,ff,mean_corr_popV_C,xs_C,[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('remap_corr_C.pdf'),600);
end
%%
if 1
    mean_corr_popV_C = selA.mean_PV_corr;
    xs_C = selA.xs;
    ff = makeFigureRowsCols(105,[1 0.5 4 0.5],'RowsCols',[4 4],...
        'spaceRowsCols',[0.1 0.1],'rightUpShifts',[0.13 0.2],'widthHeightAdjustment',...
        [-150 -150]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 3 4 4]);
    ff = show_remapping_corr_plots(mData,ff,mean_corr_popV_C,xs_C,[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('remap_corr_C.pdf'),600);
end
%%
