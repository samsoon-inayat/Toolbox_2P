function figure_place_remapping_AD(fn,allRs,ccs)

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_C = evalin('base','ei10_C1'); 
ei_A = evalin('base','ei10_A1'); 

selContexts = [1];
rasterNames = {'airD'};
trials = mat2cell([1:10]',ones(10,1));

RsC = get_rasters_data(ei_C,selContexts,rasterNames);
RsC = find_responsive_rasters(RsC,1:10);
[resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC] = get_responsive_fraction(RsC);
respC = get_cell_list(resp_valsC,[]);
out_C = find_population_vector_corr_remap_trials(RsC,respC,trials);


RsA = get_rasters_data(ei_A,selContexts,rasterNames);
RsA = find_responsive_rasters(RsA,1:10);
[resp_fractionA,resp_valsA,OIA,mean_OIA,resp_ORA,resp_OR_fractionA,resp_ANDA,resp_AND_fractionA] = get_responsive_fraction(RsA);
respA = get_cell_list(resp_valsA,[]);
out_A = find_population_vector_corr_remap_trials(RsA,respA,trials);


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
    ff = makeFigureRowsCols(106,[1 0.5 4 0.5],'RowsCols',[5 5],...
        'spaceRowsCols',[0.1 0.1],'rightUpShifts',[0.13 0.2],'widthHeightAdjustment',...
        [-150 -150]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 3 5 5]);
    ff = show_remapping_corr_plots(mData,ff,selC.mean_PV_corr,selC.xs,[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('remap_corr_C.pdf'),600);
end