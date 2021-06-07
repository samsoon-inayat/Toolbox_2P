function light_figure

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_11_15 = evalin('base','ei_11_15'); 
ei_2_3 = evalin('base','ei_2_3'); 

selContexts = [3 4 5];
rasterNames = {'airIT','airIT','airIT'};
Rs = get_rasters_data(ei_11_15,selContexts,rasterNames);
% Rs = get_rasters_data(ei_2_3,selContexts,rasterNames);
Rs = find_responsive_rasters(Rs,1:9);
mR = calc_mean_rasters(Rs,1:9);
[CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,1);
[resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC] = get_responsive_fraction(Rs);
n = 0;
%%
an = 1; cn = 1;
% plotRasters_simplest(Rs{an,cn})
% find(resp_valsC{an}(:,cn));
ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
    'spaceRowsCols',[0.15 0.03],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
    [-50 -375]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[10 4 3.25 1]);
ff = sample_rasters(Rs{an,cn},[558 328 168 55],ff);
save_pdf(ff.hf,mData.pdf_folder,sprintf('airI_rasters'),600);
%% population vector and correlation single animal
an = 1;
ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
    'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.07 0.1],'widthHeightAdjustment',...
    [0.01 -60]);
set(gcf,'color','w');
set(gcf,'Position',[5 5 3.25 2]);
ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:));
save_pdf(ff.hf,mData.pdf_folder,sprintf('airI_population_vector_corr.pdf'),600);

%% average correlation of all animals
ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 3],...
    'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.075 0.2],'widthHeightAdjustment',...
    [0.01 -220]);
set(gcf,'color','w');
set(gcf,'Position',[5 5 3.25 1]);
ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc);
save_pdf(ff.hf,mData.pdf_folder,sprintf('airI_average_population_vector_corr.pdf'),600);
%%
dataT = array2table(resp_fractionC*100);
dataT.Properties.VariableNames = {'L1','L2','L3'};
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
    xticks = xdata(1:end); xticklabels = {'C1','C4','C1'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30)
    changePosition(gca,[0.2 0.03 -0.3 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Light Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('fraction_light_responsive'),600);
%% overlap of cells
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[7 7 1.25 1],'color','w');
    hold on;
    imagesc(mean_OIC);
    axis equal
    colorbar;
    xlim([0.5 3.5]);
    ylim([0.5 3.5]);
    changePosition(gca,[0.1 0.03 0 0]);
    xticklabels = {'C1','C4','C1'''};
    set(gca,'XTick',[1 2 3],'XTickLabels',xticklabels,'YTick',[1 2 3],'YTickLabels',(xticklabels));
    set(gca,'Ydir','reverse','linewidth',0.5,'FontSize',6,'FontWeight','Bold');
    save_pdf(hf,mData.pdf_folder,sprintf('light_overlap_image'),600);
%% overlap RM bar graph
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    for ii = 1:5
        C12(ii) = OIC{ii}(1,2);
        C13(ii) = OIC{ii}(1,3);
        C23(ii) = OIC{ii}(2,3);
    end
    dataT = array2table([C12' C13' C23']);
    dataT.Properties.VariableNames = {'L1','L2','L3'};
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
    xticks = xdata(1:end); xticklabels = {'C1-C4','C1-C1''','C4-C1'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30)
    changePosition(gca,[0.2 0.03 -0.3 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Overlap Index'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('overlap_stats'),600);