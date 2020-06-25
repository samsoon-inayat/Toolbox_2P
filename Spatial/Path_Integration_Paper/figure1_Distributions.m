function figure1_Distributions

protocol = '10';
% protocol = '15';
ei = evalin('base',sprintf('ei%s',protocol));
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ET = evalin('base',sprintf('ET%s',protocol));
selAnimals = eval(sprintf('mData.selAnimals%s',protocol));

% in the following variable all the measurements are in the matrices form
% for each variable colums indicate raster and stim marker types specified 
% the rows indicate condition numbers.
paramMs = parameter_matrices('get',protocol);
% after getting all matrics, we can apply selection criteria to select a
% subgroup of cells
% here is the selection criteria in make_selC_structure function
% cellsOrNot = NaN; planeNumber = NaN; zMI_Th = 2; fwids = [0 140]; fcens = [0 140]; rs_th = NaN;
cellsOrNot = NaN; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
conditionsAndRasterTypes = [41 42 43 44 4q5]';% 21 22 23 24 25 31 32 33 34 35 41 42 43 44 45]';
% conditionsAndRasterTypes = [11 13 21 23 31 33 41 43]';
% conditionsAndRasterTypes = [11 21 31 41]';
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th);
[cpMs,pMs] = parameter_matrices('select',protocol,{paramMs,selC});
perc_cells = parameter_matrices('print',protocol,{cpMs,pMs,ET,selAnimals});

gAllVals = [];
for rr = 1:size(pMs,1)
    for cc = 1:size(pMs,2)
        tcond = conditionsAndRasterTypes(rr,cc);
        nds = dec2base(tcond,10) - '0';
        varNames{rr,cc} = sprintf('C%dR%d',nds(1),nds(2));
        for an = 1:length(selAnimals)
            zMIs{an,rr,cc} = (squeeze(pMs{rr,cc}.all_zMIs{selAnimals(an)}(nds(1),nds(2),:)));
            gAllVals = [gAllVals;zMIs{an,rr,cc}];
        end
    end
end

all_conds = []; all_rts = [];
for rr = 1:size(pMs,1)
    for cc = 1:size(pMs,2)
        tcond = conditionsAndRasterTypes(rr,cc);
        nds = dec2base(tcond,10) - '0';
        varNames{rr,cc} = sprintf('C%dR%d',nds(1),nds(2));
        xticklabels{cc,rr} = sprintf('%s-%s',paramMs.stimMarkers{nds(2)},paramMs.rasterTypes{nds(2)}(1));
        all_conds = [all_conds nds(1)]; all_rts = [all_rts nds(2)];
    end
end
all_conds = unique(all_conds); all_rts = unique(all_rts);
n = 0;

%%
runthis = 1;
if runthis
    data = zMIs;
    minBin = min(gAllVals);
    maxBin = 15;%max(gAllVals);
    incr = (maxBin-minBin)/100;
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 6 2 1.5],'color','w');
    hold on;
    [ha,hb,hca,sigR] = plotDistributions(data,'colors',colors,'maxY',90,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    hold on;
    legs = [];
    for ii = 1:length(paramMs.stimMarkers)
        if ii == 1
            legs{ii} = sprintf('%s-%s (N = %d)',paramMs.stimMarkers{ii},paramMs.rasterTypes{ii}(1),size(data,1));
        else
            legs{ii} = sprintf('%s-%s',paramMs.stimMarkers{ii},paramMs.rasterTypes{ii}(1));
        end
    end
    ylim([0 120]);
    xlims = xlim; dx = xlims(2) - xlims(1); ylims = ylim; dy = ylims(2) - ylims(1);
    legs{ii+1} = [xlims(1)+dx/1.5 dx/30 ylims(1)+dy/3 dy/15];
    if size(data,2) > 2
        putLegend(ha,legs,'colors',colors,'sigR',{sigR,'ranova',sigColor,6});
    else
        putLegend(ha,legs,'colors',colors,'sigR',{sigR,'ttest',sigColor,10});
    end
    axes(ha);
    h = xlabel('Mutual Information Z-Score');%changePosition(h,[0 -dy/3 0]);
    h = ylabel('Percentage');changePosition(h,[-0 0 0]);
    set(gca,'FontSize',axes_font_size,'FontWeight','Bold');changePosition(ha,[0.07 0.01 -0.05 -0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution Of zMI %d',all_conds(1)),600);

   
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 6 1.5 1],'color','w');
    hold on;
    xdata = 1:1.5:8; xdata = xdata(1:5);
    means = sigR.ranova.ds.avg; sems = sigR.ranova.ds.sem;
    combs = sigR.ranova.multcompare.combs; hs = sigR.ranova.multcompare.h; ps = sigR.ranova.multcompare.p;
    [hbs,maxYY] = plotBarsWithSigLines(means,sems,combs,[hs ps],'colors',colors,...
        'sigColor','k','ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 max(xdata)+0.75],'ylim',[0 maxYY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata; 
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
    save_pdf(hf,mData.pdf_folder,sprintf('Mean zMI Condition %d',all_conds(1)),600);
    
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