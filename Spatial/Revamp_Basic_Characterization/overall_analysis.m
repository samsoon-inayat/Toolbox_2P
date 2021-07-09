function overall_analysis

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei = evalin('base','ei'); 
% selContexts = [1 2 3 4 4 5 6];
% rasterNames = {'light22T','airD','light22T','airD','airD','light22T'};
% Rs = get_rasters_data(ei,selContexts,rasterNames);

selContexts = [1 4 6];
rasterNames = {'light22T','light22T','light22T'};
Rs1 = get_rasters_data(ei,selContexts,rasterNames);
% Rs1 = get_rasters_data(ei,selContexts,rasterNames);

selContexts = [2 7];
rasterNames = {'air55T','air55T'};
% Rs12 = get_rasters_data(ei,selContexts,rasterNames);
Rs12 = get_rasters_data(ei,selContexts,rasterNames);

selContexts = [3 4 5];
rasterNames = {'airD','airD','airD'};
Rs2 = get_rasters_data(ei,selContexts,rasterNames);

selContexts = [1 4 6 2 7 3 4 5];
rasterNames = {'light22T','light22T','light22T','air55T','air55T','airD','airD','airD'};
Rs = [Rs1 Rs12 Rs2];

mRs = calc_mean_rasters(Rs,1:10);
Rs = find_responsive_rasters(Rs,1:10);

[resp_fraction,resp_vals,OI,mean_OI,resp_OR,resp_OR_fraction,resp_AND,resp_AND_fraction] = get_responsive_fraction(Rs);
[CR,aCR] = find_population_vector_corr(Rs,mRs,1,0);
% return;

view_population_vector(Rs,mRs,1,100);
view_population_vector_corr(Rs,mRs,1,200);
return;
%%
if 0
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[7 7 3.25 3.25],'color','w');
hold on;
imagesc(mean_OI);
axis equal
colorbar;
xlim([0.5 8.5]);
ylim([0.5 8.5]);
changePosition(gca,[0 0.03 0 0]);
selContexts = [1 4 6 2 7 3 4 5];
rasterNames = {'light22T','light22T','light22T','air55T','air55T','airD','airD','airD'};
xticklabels = {'C1','C4L','C1''','C2','C2''','C3','C4','C3'''};
set(gca,'XTick',[1:8],'XTickLabels',xticklabels,'YTick',[1:7],'YTickLabels',(xticklabels));
set(gca,'Ydir','reverse','linewidth',0.5,'FontSize',6,'FontWeight','Bold');
save_pdf(hf,mData.pdf_folder,sprintf('overall_overlap_imageD'),600);
return;
end
%%
%% overlap RM bar graph
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    C2 = [];
    for ii = 1:5
        C2(ii,:) = OI{ii}(4,5:8);
    end
    dataT = array2table([C2]);
%     dataT.Properties.VariableNames = {'A1','A2','A3'};
    within = array2table([1 2 3 4]');
    within.Properties.VariableNames = {'Cond'};
    within.Cond = categorical(within.Cond);
    ra = repeatedMeasuresAnova(dataT,within,0.05);

%%
mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
    % row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
    % row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

    xdata = [1 2 3 4]; 
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
%     tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
    tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.05);
    for ii = 5:length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
%     maxY = maxY + 5 + 2;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'C2''','C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30)
    changePosition(gca,[0.2 0.03 -0.3 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Overlap Index'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('air_overlap_statsD'),600);
    
