function figure_place_remapping_AD_trials_big_ANOVA(fn,allRs,ccs)

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_C = evalin('base','ei10_C1'); 
ei_A = evalin('base','ei10_A1'); 

all_selContexts = [1 2 3 4];
rasterNames = {'airD'};
for ss = 1:length(all_selContexts)
    ss
    trials = mat2cell([1:10]',ones(10,1));
    selContexts = all_selContexts(ss);
    
    RsC = get_rasters_data(ei_C,selContexts,rasterNames);
    RsC = find_responsive_rasters(RsC,1:10);
    [resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC] = get_responsive_fraction(RsC);
    respC = get_cell_list(resp_valsC,[]);
    out_C{ss} = find_population_vector_corr_remap_trials(RsC,respC,trials);


    RsA = get_rasters_data(ei_A,selContexts,rasterNames);
    RsA = find_responsive_rasters(RsA,1:10);
    [resp_fractionA,resp_valsA,OIA,mean_OIA,resp_ORA,resp_OR_fractionA,resp_ANDA,resp_AND_fractionA] = get_responsive_fraction(RsA);
    respA = get_cell_list(resp_valsA,[]);
    out_A{ss} = find_population_vector_corr_remap_trials(RsA,respA,trials);
end

selC = out_C{1};
selA = out_A{1};
n = 0;
%%
all_var_C = [];
all_var_A = [];
all_var_names = [];
all_xlabels = [];
for ii = 1:length(all_selContexts)
    selC = out_C{ii};
    selA = out_A{ii};
    varC = selC.adj_SP_corr_diag;
    varA = selA.adj_SP_corr_diag;
    var_C = []; var_A = [];
    for rr = 1:size(varC,1)
        for cc = 1:size(varC,2)
            var_C(rr,cc) = nanmean(varC{rr,cc});
            var_A(rr,cc) = nanmean(varA{rr,cc});
        end
    end
    for cc = 1:size(var_C,2)
        varNames{cc} = sprintf('C%dT%d%d',ii,cc,cc+1);
        xlabels{cc} = sprintf('C%d-T%d-T%d',ii,cc,cc+1);
    end
    all_var_C = [all_var_C var_C];
    all_var_A = [all_var_A var_A];
    all_var_names = [all_var_names varNames];
    all_xlabels = [all_xlabels xlabels];
end
%%
if 1
    dataT = array2table([[ones(size(var_C,1),1);2*ones(size(var_A,1),1)] [all_var_C;all_var_A]]);
    dataT.Properties.VariableNames = {'Group',all_var_names{:}};
    dataT.Group = categorical(dataT.Group);
    colVar1 = [1:size(var_C,2)]; colVar0 = ones(size(colVar1));
    colVar1 = repmat(colVar1,1,4);    
    colVar0 = [colVar0 2*colVar0 3*colVar0 4*colVar0];
    within = table(colVar0',colVar1');
    within.Properties.VariableNames = {'Condition','TrialPairs'};
    within.TrialPairs = categorical(within.TrialPairs);
    within.Condition = categorical(within.Condition);
    ra = repeatedMeasuresAnova(dataT,within,0.05);
%     writetable(dataT,fullfile(mData.pdf_folder,'Remapping_Trials.xls'));
%%
    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
    nbars = length(mVar)/8;
%     xdata = [1:nbars ([1:nbars]+nbars+1) ([1:nbars]+nbars+1)+nbars+1 ([1:nbars]+nbars+1)+nbars+1+nbars+1];
    xdata1 = [1:9:72]'; xdata1(2:end) = xdata1(2:end) + [1:7]';
    xdata2 = [9:9:72]'; xdata2(2:end) = xdata2(2:end) + [1:7]';
    xdata = []
    for ii = 1:length(xdata1)
        xdata = [xdata xdata1(ii):xdata2(ii)];
    end
    
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[3 7 9 3],'color','w');
    hold on;
    tcolors = colors(1:nbars); tcolors = repmat(tcolors,8,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = (nbars*4+1):length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = all_xlabels;
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(45);
    changePosition(gca,[-0.09 0.03 0.19 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Correlation'},[0 0 0]});

    save_pdf(hf,mData.pdf_folder,'remap bar graph_trials_big_pop',600);
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