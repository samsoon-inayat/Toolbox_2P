function light_figure

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei = evalin('base','ei'); 


selContexts = [3 4 5];
rasterNames = {'airT','airT','airT'};
% Rs = get_rasters_data(ei_11_15,selContexts,rasterNames);
RsT = get_rasters_data(ei,selContexts,rasterNames);
RsT = find_responsive_rasters(RsT,1:10);
[resp_fractionCT,resp_valsCT,OICT,mean_OICT,resp_ORCT,resp_OR_fractionCT,resp_ANDCT,resp_AND_fractionCT] = get_responsive_fraction(RsT);
mRT = calc_mean_rasters(RsT,1:10);

selContexts = [3 4 5];
rasterNames = {'airD','airD','airD'};
% Rs = get_rasters_data(ei_11_15,selContexts,rasterNames);
Rs = get_rasters_data(ei,selContexts,rasterNames);
Rs = find_responsive_rasters(Rs,1:10);
[resp_fractionCS,resp_valsCS,OIC,mean_OICS,resp_ORCS,resp_OR_fractionCS,resp_ANDCS,resp_AND_fractionCS] = get_responsive_fraction(Rs);
mR = calc_mean_rasters(Rs,1:10);
n = 0;

%%
resp = resp_ORCS
[CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,1,0);

%%
trials = mat2cell([1:10]',ones(size([1:10]')));
resp = get_cell_list(resp_valsCS,[]);
out1 = find_population_vector_corr_remap_trials(Rs(:,1),resp,trials);
out2 = find_population_vector_corr_remap_trials(Rs(:,2),resp,trials);
out3 = find_population_vector_corr_remap_trials(Rs(:,3),resp,trials);

n = 0;

%% Show sample rasters
% an = 5; cn = 2;
% plotRasters_simplest(Rs{an,cn})
if 1
   an = 3; cn = 1;
    % plotRasters_simplest(Rs{an,cn})
    % find(resp_valsC{an}(:,cn));
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.15 0.03],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
        [-50 -375]);
    gg = 1;
    set(gcf,'color','w');
    set(gcf,'Position',[10 4 4.25 1]);
    ff = sample_rasters(Rs{an,cn},[191 11 96 41],ff);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersD'),600);
    
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.15 0.03],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
        [-50 -375]);
    gg = 1;
    set(gcf,'color','w');
    set(gcf,'Position',[10 4 4.25 1]);
    ff = sample_rasters(RsT{an,cn},[191 11 96 41],ff);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersDT'),600);
    
end
%% Mutual Information Time versus Distance Distributions

if 1
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            R = Rs{rr,cc};
            zMIsC{rr,cc} = R.info_metrics.ShannonMI_Zsh';
            R = RsT{rr,cc};
            zMIsA{rr,cc} = R.info_metrics.ShannonMI_Zsh';
        end
    end

    CN = 3;
    tcolors = {'k','r'};
    distD(:,1) = zMIsC(:,CN);
    distD(:,2) = zMIsA(:,CN);
    [distDo,allVals,allValsG] = plotDistributions(distD);
    minBin = min(allVals);
    maxBin = max(allVals);

    tcolors = {'k','r'};
    incr = 0.001; %maxBin =
    hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    %    [ha,hb,hca] = plotDistributions(distD,'colors',tcolors,'maxY',maxBin,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    [ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
    plot([1.65 1.65],[0 100],'--k');
%     xlim([-5 30]);
    set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    changePosition(gca,[0.15 0.13 -0.2 -0.13]);
    if CN == 1
        put_axes_labels(gca,{'zMI',[0 0 0]},{{'Percentage','of Neurons'},[0 0 0]});
    else
        put_axes_labels(gca,{'zMI',[0 0 0]},{'',[1 0 0]});
    end
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_zMI_%d',CN),600);
end
%% Mutual Information Time versus Distance bar graph
%%
dataT = array2table([[1;1;1;1;1;2;2;2;2;2] [mzMIsC;mzMIsA]]);
dataT.Properties.VariableNames = {'Group','C1','C2','C3','C4'};
within = array2table([1 2 3 4]');
within.Properties.VariableNames = {'Cond'};
within.Cond = categorical(within.Cond);
dataT.Group = categorical(dataT.Group);
ra = repeatedMeasuresAnova(dataT,within,0.05);
% writetable(dataT,fullfile(mData.pdf_folder,'zMI_all_cells.xls'));
%%
mVar = ra.est_marginal_means.Mean;
semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
xdata = [1 2 3 4 6 7 8 9]; 

hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.9 1],'color','w');
hold on;
tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
% tcolors = repmat(tcolors,2,1)';
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);

set(gca,'xlim',[0.25 9.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
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
put_axes_labels(gca,{[],[0 0 0]},{{'Mutual Information','(z-score)'},[0 0 0]});

save_pdf(hf,mData.pdf_folder,'zMI_bar_graph_all_cells.pdf',600);



%% population vector and correlation single animal
if 1
    an = 1;
    ff = makeFigureRowsCols(106,[1 0.5 6 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.07 0.1],'widthHeightAdjustment',...
        [0.01 -60]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 5 3.25 2]);
    [resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC,resp_exc_inh] = get_responsive_fraction(Rs);
    resp = get_cell_list(resp_valsC,[1;2]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    % ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[-0.1 1],[]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[-0.1 1],[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_population_vector_corr.pdf'),600);
end
%%
if 1
    var_1 = []; var_2 = []; var_3 = []; varNames = []; xlabels = []; xdata = [];
    for rr = 1:size(var1,1)
        for cc = 1:size(var1,2)
            var_1(rr,cc) = nanmean(var1{rr,cc});
            var_2(rr,cc) = nanmean(var2{rr,cc});
            var_3(rr,cc) = nanmean(var3{rr,cc});
        end
    end
    ind = 1; ind_val = 1;
    for ii = 1:3
        for cc = 1:size(var_1,2)
            varNames{cc+(size(var_1,2)*(ii-1))} = sprintf('C%d_T%d%d',ii,cc,cc+1);
            xlabels{cc+(size(var_1,2)*(ii-1))} = sprintf('C%d-T%d-T%d',ii,cc,cc+1);
            xdata(ind) = ind_val;
            ind = ind+1; ind_val = ind_val + 1;
        end
        ind_val = ind_val+3;
    end
    
    dataT = array2table([[var_1 var_2 var_3]]);
    dataT.Properties.VariableNames = varNames;
    colVar1 = [1:size(var_1,2)]; colVar2 = [ones(size(colVar1)) 2*ones(size(colVar1)) 3*ones(size(colVar1))];
    colVar1 = repmat(colVar1,1,3);
    
    within = array2table([colVar2' colVar1']);
    within.Properties.VariableNames = {'Condition','TrialPairs'};
    within.TrialPairs = categorical(within.TrialPairs);
    within.Condition = categorical(within.Condition);
    ra = repeatedMeasuresAnova(dataT,within);
%     writetable(dataT,fullfile(mData.pdf_folder,'Remapping_Trials.xls'));

    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
    nbars = length(mVar)/3;
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 5 1.25],'color','w');
    hold on;
    tcolors = colors(1:nbars); tcolors = repmat(tcolors,1,3);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
%     for ii = (nbars+1):length(hbs)
%         set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
%     end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = xlabels; xticklabels = repmat(xticklabels,1,2);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[-0.06 0.03 0.15 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Correlation'},[0 0 0]});

    save_pdf(hf,mData.pdf_folder,'remap bar graph_trials_air_place',600);
return;
end

%%
an = 3;
ff = makeFigureRowsCols(106,[1 0.5 6 1],'RowsCols',[2 3],...
    'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.07 0.1],'widthHeightAdjustment',...
    [0.01 -60]);
set(gcf,'color','w');
set(gcf,'Position',[5 5 3.25 2]);
ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
save_pdf(ff.hf,mData.pdf_folder,sprintf('air_population_vector_corrD.pdf'),600);


ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 3],...
    'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.075 0.2],'widthHeightAdjustment',...
    [0.01 -220]);
set(gcf,'color','w');
set(gcf,'Position',[5 3 3.25 1]);
ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[-0.2 1],[]);
save_pdf(ff.hf,mData.pdf_folder,sprintf('air_average_population_vector_corrD.pdf'),600);
return;% temporary return to just run the code up

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