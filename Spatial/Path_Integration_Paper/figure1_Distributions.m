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

%%
runthis = 1;
if runthis
    planeNumbers = 'All';
    maxDistTime = [140 10];
    contextNumber = ones(1,4).*[1 1 1 1];
    
    stimMarkers = {'air','air','belt','airI'};
    stimLabel = {'Air','Air','Belt','InterTrial'};
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
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 1.5 1.5],'color','w');
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
    legs{ii+1} = [xlims(1)+dx/1.5 dx/50 ylims(1)+dy/2 dy/15];
    if size(data,2) > 2
        putLegend(ha,legs,'colors',colors,'sigR',{sigR,'anova',sigColor,6});
    else
        putLegend(ha,legs,'colors',colors,'sigR',{sigR,'ttest',sigColor,10});
    end
    axes(ha);
    h = xlabel('Mutual Information Z-Score');%changePosition(h,[0 -dy/3 0]);
    h = ylabel('Percentage');changePosition(h,[0 0 0]);
    set(gca,'FontSize',axes_font_size,'FontWeight','Bold');changePosition(ha,[-0.04 0.09 0.13 -0.05]);
    save_pdf(hf,mData.pdf_folder,'Distribution Of zMI',600);
    if size(data,2) > 2
        sigR.anova
    else
        sigR.ttest
    end
return;
end

%%
runthis = 1;
if runthis
    varName = 'gauss_fit_on_mean.allgof.rsquare';
    planeNumbers = 'All';
    contextNumber = 1;
%     stimMarkers = {'air','belt','airOnsets010','airOffsets010'};
%     rasterTypes = {'dist','dist','time','time'};
    stimMarkers = {'air','belt','air'};
    rasterTypes = {'dist','dist','time'};
   
    for ss = 1:length(stimMarkers)
        distD = [];
        for jj = 1:length(selAnimals)
            [tempVals cns] = getParamValues(varName,ei(selAnimals(jj)),planeNumbers,contextNumber,stimMarkers{ss},rasterTypes{ss},selCells);
            distD = [distD;tempVals];
            zMI_mean(jj,ss) = nanmean(tempVals);
        end
        data{ss} = distD;
    end
    close all;
    [p,tbl,stats] = anova1(zMI_mean,1:length(stimMarkers),'on');
    figure(101);[c,~,~,gnames] = multcompare(stats,'CType','bonferroni');
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 6.9 2.5],'color','w');
    hold on;
    [ha,hb,hca,sigR] = plotDistributions(data,'colors',colors,'maxY',90,'cumPos',[0.5 0.26 0.25 0.5],'min',0,'incr',0.1,'max',1);
    hold on;
    legs = [];
    for ii = 1:length(stimMarkers)
        legs{ii} = sprintf('%s-%s (N = %d)',stimMarkers{ii},rasterTypes{ii},length(data{ii}));
    end
    legs{ii+1} = [0.3 0.05 70 5];
    putLegend(ha,legs,'colors',colors,'sigR',{sigR,'ks',sigColor,8});
    axes(ha);h = xlabel('Goodness of fit R^2');changePosition(h,[0 -3 0]);
    set(gca,'FontSize',axes_font_size,'FontWeight','Bold');
    changePosition(h,[0 0.7 0]);h = ylabel('Percentage');changePosition(h,[0.0 0 0]);
    axes(hca);set(gca,'FontSize',6);
    set(gca,'FontSize',axes_font_size,'FontWeight','Bold');changePosition(ha,[-0.04 0.09 0.13 -0.05]);
    save_pdf(hf,mData.pdf_folder,'Distribution Of R2',600);
return;
end

%%
% bar graph from anova
% ff = makeFigureWindow__one_axes_only(2,[1 4 1.75 2],[0.2 0.35 0.7 0.6]);
ff = makeFigureWindow__one_axes_only(3,[2 4 2 2],[0.25 0.17 0.7 0.79]);
set(gcf,'color','w');
set(gcf,'Position',[25 4 1.5 2]);
axes(ff.ha);
plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols_all(selGroups),'ySpacingFactor',10);
xlim([0.4 0.6+length(selGroups)]);
%     ylim([0 max(mVals)]);
hyl = ylabel('Average MI (z-score)');
pos = get(hyl,'Position');pos = pos + [+0.1 0 0];set(hyl,'Position',pos);
set(ff.ha,'linewidth',1);
set(ff.ha,'TickDir','out','FontSize',8,'FontWeight','bold');
set(ff.ha,'XTick',[1 2 3 4],'XTickLabel',{'Context 1','Context 2','Context 3','Context4'});
xtickangle(25);
save2pdf('MI_anova_bars.pdf',ff.hf,600);



