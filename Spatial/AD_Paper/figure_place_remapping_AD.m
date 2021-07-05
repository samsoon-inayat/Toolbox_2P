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

selC = out_C;
selA = out_A;
n = 0;
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

%%
if 1
    for rr = 1:size(varC,1)
        for cc = 1:size(varC,2)
            var_C(rr,cc) = nanmean(varC{rr,cc});
            var_A(rr,cc) = nanmean(varA{rr,cc});
        end
    end
    dataT = array2table([[ones(size(var_C,1),1);2*ones(size(var_A,1),1)] [var_C;var_A]]);
    dataT.Properties.VariableNames = {'Group','C12','C23','C34'};
    dataT.Group = categorical(dataT.Group);
    colVar1 = [1 2 3];    
    within = table(colVar1');
    within.Properties.VariableNames = {'Condition'};
    within.Condition = categorical(within.Condition);
    ra = repeatedMeasuresAnova(dataT,within);
    
    dataTC = dataT(1:3,2:end);
    raC = repeatedMeasuresAnova(dataTC,within);
    
    dataTA = dataT(4:end,2:end);
    raA = repeatedMeasuresAnova(dataTA,within);

    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
    xdata = [1 2 3 5:7];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
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
    changePosition(gca,[0.06 0.03 0.02 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Correlation'},[0 0 0]});

    save_pdf(hf,mData.pdf_folder,'remap bar graph',600);
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
    
    
end
%%
