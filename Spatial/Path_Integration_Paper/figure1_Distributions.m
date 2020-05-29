function figure1_Distributions

% dataAir = evalin('base','data');
% dataBelt = evalin('base','datab');
% dataAOn = evalin('base','dataAOn010');
% dataAOff = evalin('base','dataAOff010');
ei = evalin('base','ei10');
mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;
axes_font_size = mData.axes_font_size;
% allCells = mData.allCells;
selAnimals = [1:4 9];
% selAnimals = 5:8;
% selAnimals = 1:8;
n = 0;
cN = 1;
%%
runthis = 1;
if runthis
%     cN = 1;
    planeNumbers = 'All';
    maxDistTime = [142 15];
    contextNumber = ones(1,4).*cN;
    
    stimMarkers = {'air','air','belt','airI'};
    stimLabel = {'Air','Air','Belt','InterTrial'};
    xstimLabel = {'AirD','AirT','BeltD','AirIT'};
    rasterTypes = {'dist','time','dist','time'};
    rasterLabel = {'Distance','Time','Distance','Time'};
%     stimMarkers = {'air','air','air','air'};
%     rasterTypes = {'dist','dist','dist','dist'};
    selCells = 'areCells';
    allVals = []; gAllVals = [];
    for ss = 1:length(stimMarkers)
        distD = [];
        for jj = 1:length(selAnimals)
            [pcs cns areCells] = getParamValues('placeCells3',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),'air','dist',selCells,maxDistTime);
            [pcs cns areCells] = getParamValues('placeCells3',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),'air','time',selCells,maxDistTime);
            pcs = logical(ones(size(pcs)));
            [tempVals cns ACs] = getParamValues('info_metrics.ShannonMI_Zsh',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
%             [tempVals cns ACs] = getParamValues('gauss_fit_on_mean.rs',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
%             [tempVals cns ACs] = getParamValues('data.rs',ei(selAnimals(jj)),planeNumbers,contextNumber,stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
%             [temppws cns ACs] = getParamValues('data.mean_trial_corr',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
%             tempVals = temppws;
%             tempVals = tempVals(logical(pcs));
            distD = [distD;tempVals(pcs)];
            zMI_mean(jj,ss) = nanmean(tempVals(pcs));
            data{jj,ss} = tempVals(pcs);
        end
        allVals{ss} = distD;
        gAllVals = [gAllVals;distD];
    end
    minBin = min(gAllVals);
    maxBin = 15;%max(gAllVals);
    incr = (maxBin-minBin)/100;
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 6 2 1.5],'color','w');
    hold on;
    [ha,hb,hca,sigR] = plotDistributions(data,'colors',colors,'maxY',90,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    hold on;
    legs = [];
    for ii = 1:length(stimMarkers)
        if ii == 1
            legs{ii} = sprintf('%s-%s (N = %d)',stimLabel{ii},rasterLabel{ii},size(data,1));
        else
            legs{ii} = sprintf('%s-%s',stimLabel{ii},rasterLabel{ii});
        end
    end
    ylim([0 120]);
    xlims = xlim; dx = xlims(2) - xlims(1); ylims = ylim; dy = ylims(2) - ylims(1);
    legs{ii+1} = [xlims(1)+dx/2 dx/30 ylims(1)+dy/3 dy/15];
    if size(data,2) > 2
        putLegend(ha,legs,'colors',colors,'sigR',{sigR,'anova',sigColor,6});
    else
        putLegend(ha,legs,'colors',colors,'sigR',{sigR,'ttest',sigColor,10});
    end
    axes(ha);
    h = xlabel('Mutual Information Z-Score');%changePosition(h,[0 -dy/3 0]);
    h = ylabel('Percentage');changePosition(h,[-0 0 0]);
    set(gca,'FontSize',axes_font_size,'FontWeight','Bold');changePosition(ha,[0.07 0.01 -0.05 -0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution Of zMI %d',cN),600);
    if size(data,2) > 2
        sigR.anova
    else
        sigR.ttest
    end
    
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 1.5 1],'color','w');
    hold on;
    xdata = 1:1.5:6; xdata = xdata(1:4);
    maxY = max(sigR.means + sigR.sems)+3;
    hbs = plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[sigR.anova.multcompare.h sigR.anova.multcompare.p],'colors',colors,'sigColor','k',...
    'maxY',maxY,'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 max(xdata)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata; xticklabels = xstimLabel;
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.13 0.13 -0.1 -0.1]);
    put_axes_labels(gca,{'Raster Type',[0 0 0]},{{'Mutual Information','(z score)'},[-0.1 -0.5 0]});
    xlims = xlim; dx = xlims(2) - xlims(1); ylims = ylim; dy = ylims(2) - ylims(1);
    legs{ii+1} = [xlims(1)+dx/2 dx/30 ylims(1)+dy/3 dy/15];
%     if size(data,2) > 2
%         putLegend(ha,legs,'colors',colors,'sigR',{sigR,'anova',sigColor,6});
%     else
%         putLegend(ha,legs,'colors',colors,'sigR',{sigR,'ttest',sigColor,10});
%     end
    save_pdf(hf,mData.pdf_folder,sprintf('Mean zMI Condition %d',cN),600);
    
return;
end



%%
runthis = 0;
if runthis
%     cN = 1;
    planeNumbers = 'All';
    maxDistTime = [142 15];
    contextNumber = ones(1,4).*cN;
    
    stimMarkers = {'air','air','belt','airI'};
    stimLabel = {'Air','Air','Belt','InterTrial'};
    xstimLabel = {'Air-D','Air-T','Belt-D','AirI-T'};
    rasterTypes = {'dist','time','dist','time'};
    rasterLabel = {'Distance','Time','Distance','Time'};
%     stimMarkers = {'air','air','air','air'};
%     rasterTypes = {'dist','dist','dist','dist'};
    selCells = 'areCells';
    allVals = []; gAllVals = [];
    for ss = 1:length(stimMarkers)
        distD = [];
        for jj = 1:length(selAnimals)
            [pcs cns areCells] = getParamValues('placeCells1',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),'air','dist',selCells,maxDistTime);
            [pcs cns areCells] = getParamValues('placeCells1',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),'air','time',selCells,maxDistTime);
            pcs = logical(ones(size(pcs)));
            [tempVals cns ACs] = getParamValues('gauss_fit_on_mean.rs',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
            distD = [distD;tempVals(pcs)];
            zMI_mean(jj,ss) = nanmean(tempVals(pcs));
            data{jj,ss} = tempVals(pcs);
        end
        allVals{ss} = distD;
        gAllVals = [gAllVals;distD];
    end
    minBin = min(gAllVals);
    maxBin = max(gAllVals);
    incr = (maxBin-minBin)/100;
    %%
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 2 1.5],'color','w');
    hold on;
    [ha,hb,hca,sigR] = plotDistributions(data,'colors',colors,'maxY',90,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    hold on;
    legs = [];
    for ii = 1:length(stimMarkers)
        if ii == 1
            legs{ii} = sprintf('%s-%s (N = %d)',stimLabel{ii},rasterLabel{ii},size(data,1));
        else
            legs{ii} = sprintf('%s-%s',stimLabel{ii},rasterLabel{ii});
        end
    end
    ylim([0 120]);
    xlims = xlim; dx = xlims(2) - xlims(1); ylims = ylim; dy = ylims(2) - ylims(1);
    legs{ii+1} = [xlims(1)+dx/5 dx/30 ylims(1)+dy/1 dy/15];
    if size(data,2) > 2
        putLegend(ha,legs,'colors',colors,'sigR',{sigR,'anova',sigColor,6});
    else
        putLegend(ha,legs,'colors',colors,'sigR',{sigR,'ttest',sigColor,10});
    end
    axes(ha);
    h = xlabel('Rsq');%changePosition(h,[0 -dy/3 0]);
    h = ylabel('Percentage');changePosition(h,[-0 0 0]);
    set(gca,'FontSize',axes_font_size,'FontWeight','Bold');changePosition(ha,[0.07 0.01 -0.05 -0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution Of Rsq %d',cN),600);
    if size(data,2) > 2
        sigR.anova
    else
        sigR.ttest
    end
    %%
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 1.5 1],'color','w');
    hold on;
    xdata = 1:1.5:6; xdata = xdata(1:4);
    maxY = max(sigR.means + sigR.sems)+0.2;
    hbs = plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[sigR.anova.multcompare.h sigR.anova.multcompare.p],'colors',colors,'sigColor','k',...
    'maxY',maxY,'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 max(xdata)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata; xticklabels = xstimLabel;
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.13 0.13 -0.1 -0.1]);
    put_axes_labels(gca,{'Raster Type',[0 0 0]},{{'Rsq'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Mean Rsq Condition %d',cN),600);
    
return;
end

%%
runthis = 0;
if runthis
%     cN = 1;
    planeNumbers = 'All';
    maxDistTime = [142 15];
    contextNumber = ones(1,4).*cN;
    
    stimMarkers = {'air','air','belt','airI'};
    stimLabel = {'Air','Air','Belt','InterTrial'};
    xstimLabel = {'Air-D','Air-T','Belt-D','AirI-T'};
    rasterTypes = {'dist','time','dist','time'};
    rasterLabel = {'Distance','Time','Distance','Time'};
%     stimMarkers = {'air','air','air','air'};
%     rasterTypes = {'dist','dist','dist','dist'};
    selCells = 'areCells';
    allVals = []; gAllVals = [];
    for ss = 1:length(stimMarkers)
        distD = [];
        for jj = 1:length(selAnimals)
            [pcs cns areCells] = getParamValues('placeCells1',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),'air','dist',selCells,maxDistTime);
            [pcs cns areCells] = getParamValues('placeCells1',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),'air','time',selCells,maxDistTime);
            pcs = logical(ones(size(pcs)));
            [tempVals cns ACs] = getParamValues('fractal_dim.HaFD',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
            distD = [distD;tempVals(pcs)];
            zMI_mean(jj,ss) = nanmean(tempVals(pcs));
            data{jj,ss} = tempVals(pcs);
        end
        allVals{ss} = distD;
        gAllVals = [gAllVals;distD];
    end
    minBin = min(gAllVals);
    maxBin = max(gAllVals);
    incr = (maxBin-minBin)/100;
    %%
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 2 1.5],'color','w');
    hold on;
    [ha,hb,hca,sigR] = plotDistributions(data,'colors',colors,'maxY',90,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    hold on;
    legs = [];
    for ii = 1:length(stimMarkers)
        if ii == 1
            legs{ii} = sprintf('%s-%s (N = %d)',stimLabel{ii},rasterLabel{ii},size(data,1));
        else
            legs{ii} = sprintf('%s-%s',stimLabel{ii},rasterLabel{ii});
        end
    end
    ylim([0 120]);
    xlims = xlim; dx = xlims(2) - xlims(1); ylims = ylim; dy = ylims(2) - ylims(1);
    legs{ii+1} = [xlims(1)+dx/5 dx/30 ylims(1)+dy/1 dy/15];
    if size(data,2) > 2
        putLegend(ha,legs,'colors',colors,'sigR',{sigR,'anova',sigColor,6});
    else
        putLegend(ha,legs,'colors',colors,'sigR',{sigR,'ttest',sigColor,10});
    end
    axes(ha);
    h = xlabel('Hausdorff FD');%changePosition(h,[0 -dy/3 0]);
    h = ylabel('Percentage');changePosition(h,[-0 0 0]);
    set(gca,'FontSize',axes_font_size,'FontWeight','Bold');changePosition(ha,[0.07 0.01 -0.05 -0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution Of HaFD %d',cN),600);
    if size(data,2) > 2
        sigR.anova
    else
        sigR.ttest
    end
    %%
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 1.5 1],'color','w');
    hold on;
    xdata = 1:1.5:6; xdata = xdata(1:4);
    maxY = max(sigR.means + sigR.sems)+0.5
    hbs = plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[sigR.anova.multcompare.h sigR.anova.multcompare.p],'colors',colors,'sigColor','k',...
    'maxY',maxY,'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 max(xdata)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata; xticklabels = xstimLabel;
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.13 0.13 -0.1 -0.1]);
    put_axes_labels(gca,{'Raster Type',[0 0 0]},{{'Hausdorff FD'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Mean HaFD Condition %d',cN),600);
return;
end

%%
runthis = 1;
if runthis
%     cN = 1;
    planeNumbers = 'All';
    maxDistTime = [142 15];
    contextNumber = ones(1,4).*cN;
    
    stimMarkers = {'air','air','belt','airI'};
    stimLabel = {'Air','Air','Belt','InterTrial'};
    xstimLabel = {'Air-D','Air-T','Belt-D','AirI-T'};
    rasterTypes = {'dist','time','dist','time'};
    rasterLabel = {'Distance','Time','Distance','Time'};
%     stimMarkers = {'air','air','air','air'};
%     rasterTypes = {'dist','dist','dist','dist'};
    selCells = 'areCells';
    allVals = []; gAllVals = [];
    for ss = 1:length(stimMarkers)
        distD = [];
        for jj = 1:length(selAnimals)
            [pcs cns areCells] = getParamValues('placeCells1',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),'air','dist',selCells,maxDistTime);
            [pcs cns areCells] = getParamValues('placeCells1',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),'air','time',selCells,maxDistTime);
            pcs = logical(ones(size(pcs)));
            [tempVals cns ACs] = getParamValues('fractal_dim.HiFD',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
            distD = [distD;tempVals(pcs)];
            zMI_mean(jj,ss) = nanmean(tempVals(pcs));
            data{jj,ss} = tempVals(pcs);
        end
        allVals{ss} = distD;
        gAllVals = [gAllVals;distD];
    end
    minBin = min(gAllVals);
    maxBin = max(gAllVals);
    incr = (maxBin-minBin)/100;
    %%
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 2 1.5],'color','w');
    hold on;
    [ha,hb,hca,sigR] = plotDistributions(data,'colors',colors,'maxY',90,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    hold on;
    legs = [];
    for ii = 1:length(stimMarkers)
        if ii == 1
            legs{ii} = sprintf('%s-%s (N = %d)',stimLabel{ii},rasterLabel{ii},size(data,1));
        else
            legs{ii} = sprintf('%s-%s',stimLabel{ii},rasterLabel{ii});
        end
    end
    ylim([0 120]);
    xlims = xlim; dx = xlims(2) - xlims(1); ylims = ylim; dy = ylims(2) - ylims(1);
    legs{ii+1} = [xlims(1)+dx/5 dx/30 ylims(1)+dy/1 dy/15];
    if size(data,2) > 2
        putLegend(ha,legs,'colors',colors,'sigR',{sigR,'anova',sigColor,6});
    else
        putLegend(ha,legs,'colors',colors,'sigR',{sigR,'ttest',sigColor,10});
    end
    axes(ha);
    h = xlabel('Higuchi FD');%changePosition(h,[0 -dy/3 0]);
    h = ylabel('Percentage');changePosition(h,[-0 0 0]);
    set(gca,'FontSize',axes_font_size,'FontWeight','Bold');changePosition(ha,[0.07 0.01 -0.05 -0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution Of HiFD %d',cN),600);
    if size(data,2) > 2
        sigR.anova
    else
        sigR.ttest
    end
    %%
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 1.5 1],'color','w');
    hold on;
    xdata = 1:1.5:6; xdata = xdata(1:4);
    maxY = max(sigR.means + sigR.sems)+0.75;
    hbs = plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[sigR.anova.multcompare.h sigR.anova.multcompare.p],'colors',colors,'sigColor','k',...
    'maxY',maxY,'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 max(xdata)+0.75],'ylim',[1.8 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata; xticklabels = xstimLabel;
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.13 0.13 -0.1 -0.1]);
    put_axes_labels(gca,{'Raster Type',[0 0 0]},{{'Higuchi FD'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Mean HiFD Condition %d',cN),600);
return;
end