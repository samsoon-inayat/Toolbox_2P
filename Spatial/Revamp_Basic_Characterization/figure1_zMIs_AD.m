function figure1_number_of_PCs

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei = evalin('base','d15'); 
selContexts = [1 2 3 3 4 4 5 5 6 7];
rasterNames = {'light22T','air55T','air77T','airD','air77T','airD','air77T','airD','light22T','air55T'};
Rs = get_rasters_data(ei,selContexts,rasterNames);

an = 1;
cai_sampling_rate = ei{an}.thorExp.frameRate;
effective_sampling_rate = 1/0.15;
samplingRate = {'Ca','Ef','Ef','Ef','Ef','Ef','Ef','Ef','Ca','Ef'};
timeBefore = [2 5 7 NaN 7 NaN 7 NaN 2 5];

inds = [4 6 8];
for anii = 1:5
    for indsii = 1:length(inds)
        this_raster = Rs{anii,inds(indsii)};
        zMIs = this_raster.info_metrics.ShannonMI_Zsh;
%         zMIs = zMIs(zMIs>1.96);
        var_oi_C(anii,indsii) = nanmean(zMIs);
    end
end
n = 0;
%%
runthis = 1;
if runthis
    dataT = array2table(var_oi_C);
    dataT.Properties.VariableNames = {'C3','C4','C3p'};
    within = table([1 2 3]');
    within.Properties.VariableNames = {'Condition'};
    within.Condition = categorical(within.Condition);
    ra = repeatedMeasuresAnova(dataT,within);
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
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 5:length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
%     maxY = maxY + 5 + 2;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30)
    changePosition(gca,[0.02 0.03 -0.04 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{'zMI',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('zMIs_bargraph'),600);
return;
end

