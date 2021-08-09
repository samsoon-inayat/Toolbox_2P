function figure_place_remapping_AD_trials_big_ANOVA(fn,allRs,ccs)

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_C = evalin('base','ei10_C1'); 
ei_A = evalin('base','ei10_A1'); 

all_selContexts = [1 2 3 4];
rasterNames = {'airD'};
tic
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
toc
n = 0;
%%
ntrials = 9;
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
    all_var_C = [all_var_C var_C(:,1:ntrials)];
    all_var_A = [all_var_A var_A(:,1:ntrials)];
    all_var_names = [all_var_names varNames(1:ntrials)];
    all_xlabels = [all_xlabels xlabels];
end
%%
if 1
    [within,dvn,xlabels] = make_within_table({'Cond','TrialPairs'},[4 ntrials]);
    dataT = array2table([[ones(size(var_C,1),1);2*ones(size(var_A,1),1)] [all_var_C;all_var_A]]);
    dataT.Properties.VariableNames = {'Group',all_var_names{:}};
    dataT.Group = categorical(dataT.Group);
%     colVar1 = [1:size(var_C,2)]; colVar0 = ones(size(colVar1));
%     colVar1 = repmat(colVar1,1,4);    
%     colVar0 = [colVar0 2*colVar0 3*colVar0 4*colVar0];
%     within = table(colVar0',colVar1');
%     within.Properties.VariableNames = {'Condition','TrialPairs'};
%     within.TrialPairs = categorical(within.TrialPairs);
%     within.Condition = categorical(within.Condition);
    ra = repeatedMeasuresAnova(dataT,within,0.05);
%     writetable(dataT,fullfile(mData.pdf_folder,'Remapping_Trials.xls'));
%%
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
   
    colors = mData.colors;
    hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[3 7 9 3],'color','w');
    hold on;
    tcolors = colors(1:ntrials); tcolors = repmat(tcolors,8,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = (1+length(xdata)/2):length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = all_xlabels;
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(45);
    changePosition(gca,[-0.09 0.03 0.15 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Correlation'},[0 0 0]});

    save_pdf(hf,mData.pdf_folder,'remap bar graph_trials_big_pop',600);
return;
end

%% Correlations across trials
cn = 1;
selC = out_C{cn};
selA = out_A{cn};
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
    ra.ranova
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
    colors = mData.colors;
    hf = get_figure(5,[3 7 5 1]);
    tcolors = mData.colors(1:ntrials); tcolors = repmat(tcolors,1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp(ci),'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'T1-T2','T2-T3','T3-T4','T4-T5','T5-T6','T6-T7','T7-T8','T8-T9','T9-T10'};
    xticklabels = repmat(xticklabels(1:ntrials),1,2);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    for ii = (1+length(hbs)/2):length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    xtickangle(45);
    changePosition(gca,[0.1 0.03 -0.1 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{typeCorr{ci},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_correlation_%d',FF{ci},cn),600);
end
end

%% average correlation of all animals
if 1
    ff = makeFigureRowsCols(106,[1 0.5 9 9],'RowsCols',[10 10],...
        'spaceRowsCols',[0.02 0.02],'rightUpShifts',[0.1 0.1],'widthHeightAdjustment',...
        [-40 -40]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 1 9 9]);
    ff = show_remapping_corr_plots(mData,ff,selC.mean_PV_corr,selC.xs,[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('remap_corr_C.pdf'),600);
end
%% average correlation of all animals
if 1
    ff = makeFigureRowsCols(107,[1 0.5 9 9],'RowsCols',[10 10],...
        'spaceRowsCols',[0.02 0.02],'rightUpShifts',[0.1 0.1],'widthHeightAdjustment',...
        [-40 -40]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 1 9 9]);
    ff = show_remapping_corr_plots(mData,ff,selA.mean_PV_corr,selA.xs,[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('remap_corr_C.pdf'),600);
end
%% testing the HaoRan's correlation of correlation
if 1
    figure(100);clf;
    thisCorrC = squeeze(mean(selC.all_auto_SP_corr_corr,[1 2]));
    thisCorrA = squeeze(mean(selA.all_auto_SP_corr_corr,[1 2]));
    subplot 121;
    imagesc(thisCorr);colorbar;
    thisCorr = mean(selA.all_auto_SP_corr_corr,3);
    subplot 122;
    imagesc(thisCorr);colorbar;
end
%% compare correlations with 5 next trials
if 1
    Nt = length(trials);
%     mask = ones(Nt,Nt); 
%     order = 1; mask = triu(mask,order) & ~triu(mask,order+1);
    trialNum = 1;
    mask = zeros(Nt,Nt); mask(trialNum,:) = 1; mask = mask.*triu(ones(size(mask)),1);

    cn = 1;
    selC = out_C{cn};
    selA = out_A{cn};
    mapvcC = []; mapvcA = [];
    for ii = 1
        [apvc] = get_corr_summary(selC.all_SP_corr,mask);
        tapvc = (arrayfun(@(x) mean(x{1}),apvc(:,1:end)));
        mapvcC = [mapvcC tapvc];
        [apvc] = get_corr_summary(selA.all_SP_corr,mask);
        tapvc = (arrayfun(@(x) mean(x{1}),apvc(:,1:end)));
        mapvcA = [mapvcA tapvc];
    end
    N = size(mapvcC,2);
    [within,dvn,xlabels] = make_within_table({'TrialPairs'},[N]);
    dataT = make_between_table({mapvcC;mapvcA},dvn);
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    ra.ranova
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 2 1]);
    hf = get_figure(5,[8 7 5 1]);
    tcolors = mData.colors(1:N); tcolors = repmat(tcolors,1,2);%[s.m;s.c;s.y];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);
    format_axes(gca);
    xticks = xdata; xticklabels = xlabels;
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.03 0.02 0 -0.011])
    put_axes_labels(gca,{[],[0 0 0]},{'Correlation',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('sp_corr_next_five_trials.pdf'),600);
end
