function light_figure

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_11_15 = evalin('base','ei'); 

selContexts = [2 7];
rasterNames = {'air55T','air55T'};
Rs = get_rasters_data(ei_11_15,selContexts,rasterNames);
% Rs = get_rasters_data(ei_2_3,selContexts,rasterNames);
Rs = find_responsive_rasters(Rs,1:10);
mR = calc_mean_rasters(Rs,1:10);
[CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,1);
[resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC] = get_responsive_fraction(Rs);
n = 0;
%%
% an = 5; cn = 2;
% plotRasters_simplest(Rs{an,cn})
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
ff = sample_rasters(Rs{an,cn},[140 66 162 181],ff);
save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rasters'),600);
%% population vector and correlation single animal
an = 1;
ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 2],...
    'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.15 0.1],'widthHeightAdjustment',...
    [-70 -60]);
set(gcf,'color','w');
set(gcf,'Position',[5 5 2.2 2]);
ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:));
save_pdf(ff.hf,mData.pdf_folder,sprintf('air_population_vector_corr.pdf'),600);

%% average correlation of all animals
ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 2],...
    'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.15 0.2],'widthHeightAdjustment',...
    [-70 -220]);
set(gcf,'color','w');
set(gcf,'Position',[5 5 2.2 1]);
ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc);
save_pdf(ff.hf,mData.pdf_folder,sprintf('air_average_population_vector_corr.pdf'),600);

%%
dataT = array2table(resp_fractionC*100);
dataT.Properties.VariableNames = {'A1','A2'};
within = array2table([1 2]');
within.Properties.VariableNames = {'Cond'};
within.Cond = categorical(within.Cond);
ra = repeatedMeasuresAnova(dataT,within,0.05);

%%
mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
% row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
% row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

xdata = [1 2]; 
colors = mData.colors;
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
hold on;
%     tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4};};
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.05);
for ii = 5:length(hbs)
set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
end
% plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
%     maxY = maxY + 5 + 2;
set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
xticks = xdata(1:end); xticklabels = {'C2','C2'''};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30)
changePosition(gca,[0.2 0.03 -0.3 -0.05]);
put_axes_labels(gca,{[],[0 0 0]},{{'Air Responsive','Cells (%)'},[0 0 0]});
save_pdf(hf,mData.pdf_folder,sprintf('percentage_air_responsive'),600);
%%
%% overlap RM bar graph
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    for ii = 1:5
        C12(ii) = OIC{ii}(1,2);
    end
%%
    [mVar,semVar] = findMeanAndStandardError(C12);
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
    xticks = xdata(1:end); xticklabels = {'C2-C2'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30)
    changePosition(gca,[0.3 0.03 -0.5 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Overlap Index'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('air_overlap_stats'),600);