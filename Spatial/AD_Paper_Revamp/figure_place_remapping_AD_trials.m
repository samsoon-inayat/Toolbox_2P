function figure_place_remapping_AD(fn,allRs,ccs)

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_C = evalin('base','ei10_C'); 
ei_A = evalin('base','ei10_A'); 

selContexts = [1];
rasterNames = {'airD'};
trials = mat2cell([1:10]',ones(10,1));
%%
tic
RsC = Rs_C(:,1); respC = sel_pop_C(:,1);
out_C = find_population_vector_corr_remap_trials(RsC,respC,trials);

RsA = Rs_C(:,1); respA = sel_pop_A(:,1);
out_A = find_population_vector_corr_remap_trials(RsA,respA,trials);
toc
%%

% 
% RsC = get_rasters_data(ei_C,selContexts,rasterNames);
% RsC = repmat(RsC,1,length(trials));
% mRsC = calc_mean_rasters(RsC,1:10);
% for ii = 1:length(trials)
%     [mRsCT(:,ii),RsCT(:,ii)] = calc_mean_rasters(RsC(:,1),trials{ii});
% end
% RsC = find_responsive_rasters(RsC,1:10);
% [CR_C,aCR_C] = find_population_vector_corr(RsC,mRsC,1,0);
% [resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC] = get_responsive_fraction(RsC);
% 
% respC_pop = get_cell_list(resp_valsC,[]);
% respC = get_cell_list(resp_valsC,[1]);
% 
% 
% RsA = get_rasters_data(ei_A,selContexts,rasterNames);
% RsA = repmat(RsA,1,length(trials));
% mRsA = calc_mean_rasters(RsA,1:10);
% for ii = 1:length(trials)
%     [mRsAT(:,ii),RsAT(:,ii)] = calc_mean_rasters(RsA(:,1),trials{ii});
% end
% RsA = find_responsive_rasters(RsA,1:10);
% [CR_A,aCR_A] = find_population_vector_corr(RsA,mRsA,1,0);
% [resp_fractionA,resp_valsA,OIA,mean_OIA,resp_ORA,resp_OR_fractionA,resp_ANDA,resp_AND_fractionA] = get_responsive_fraction(RsA);
% 
% respA_pop = get_cell_list(resp_valsA,[]);
% respA = get_cell_list(resp_valsA,[1]);
% 
% out_C = find_population_vector_corr_remap_trials(RsC,mRsCT,respC);
% out_A = find_population_vector_corr_remap(RsA,mRsAT,respA);
% 
% % out_C_pop = find_population_vector_corr_remap(RsC,mRsCT,respC_pop);
% % out_A_pop = find_population_vector_corr_remap(RsA,mRsAT,respA_pop);

selC = out_C;
selA = out_A;
n = 0;
%%
varC = selC.adj_SP_corr_diag;
varA = selA.adj_SP_corr_diag;
if 1
    CN = 5;
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
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_SP_Corr_%d_%d',CN,selContexts),600);
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
    for cc = 1:size(var_C,2)
        varNames{cc} = sprintf('T%d%d',cc,cc+1);
        xlabels{cc} = sprintf('T%d-T%d',cc,cc+1);
    end
    dataT = array2table([[ones(size(var_C,1),1);2*ones(size(var_A,1),1)] [var_C;var_A]]);
    dataT.Properties.VariableNames = {'Group',varNames{:}};
    dataT.Group = categorical(dataT.Group);
    colVar1 = [1:size(var_C,2)];    
    within = table(colVar1');
    within.Properties.VariableNames = {'TrialPairs'};
    within.TrialPairs = categorical(within.TrialPairs);
    ra = repeatedMeasuresAnova(dataT,within);
%     writetable(dataT,fullfile(mData.pdf_folder,'Remapping_Trials.xls'));

    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
    nbars = length(mVar)/2;
    xdata = [1:nbars ([1:nbars]+nbars+1)];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 5 1.25],'color','w');
    hold on;
    tcolors = colors(1:nbars); tcolors = repmat(tcolors,1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = (nbars+1):length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = xlabels; xticklabels = repmat(xticklabels,1,2);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[-0.06 0.03 0.15 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Correlation'},[0 0 0]});

    save_pdf(hf,mData.pdf_folder,'remap bar graph_trials',600);
return;
end


%% average correlation of all animals
if 1
    ff = makeFigureRowsCols(106,[1 0.5 9 9],'RowsCols',[9 9],...
        'spaceRowsCols',[0.1 0.1],'rightUpShifts',[0.13 0.2],'widthHeightAdjustment',...
        [0 0]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 3 9 9]);
    ff = show_remapping_corr_plots(mData,ff,selC.mean_PV_corr,selC.xs,[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('remap_corr_C.pdf'),600);
end

%% Correlations across trials
typeCorr = {'Spatial Correlation','Pop. Vec. Correlation','\Delta FR Score'};
FF = {'SP','PV','RR'};
ysp = [0.05 0.05 0.1];
ntrials = 9;
for ci = 1%:3;
if 1
    [within,dvn,xlabels] = make_within_table({'Cond'},ntrials);
    switch ci
        case 1
            var_C = arrayfun(@(x) mean(x{1}),selC.adj_SP_corr_diag);var_A = arrayfun(@(x) mean(x{1}),selA.adj_SP_corr_diag);
        case 2
            var_C = arrayfun(@(x) mean(x{1}),selC.adj_PV_corr_diag);var_A = arrayfun(@(x) mean(x{1}),selA.adj_PV_corr_diag);
        case 3
            var_C = arrayfun(@(x) nanmean(x{1}),selC.adj_RR_SP);var_A = arrayfun(@(x) nanmean(x{1}),selA.adj_RR_SP);
    end
    dataT = make_between_table({var_C(:,1:ntrials);var_A(:,1:ntrials)},dvn);
    ra = repeatedMeasuresAnova(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
    colors = mData.colors;
    hf = get_figure(5,[3 7 5 1]);
    tcolors = mData.colors(1:9); tcolors = repmat(tcolors,1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp(ci),'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'T1-T2','T2-T3','T3-T4','T4-T5','T5-T6','T6-T7','T7-T8','T8-T9','T9-T10'};
    xticklabels = repmat(xticklabels,1,2);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    for ii = 10:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    xtickangle(45);
    changePosition(gca,[0.1 0.03 -0.1 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{typeCorr{ci},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_correlation',FF{ci}),600);
end
end