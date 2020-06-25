function figure1_Distributions_1

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
cellsOrNot = NaN; planeNumber = NaN; zMI_Th = 3; fwids = [0 140]; fcens = [0 140]; rs_th = 0.4;
cellsOrNot = NaN; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
conditionsAndRasterTypes = [11 12 13 14 15];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th);
[cpMs,pMs] = parameter_matrices('select',protocol,{paramMs,selC});
perc_cells = parameter_matrices('print',protocol,{cpMs,pMs,ET,selAnimals});

for rr = 1:size(pMs,1)
    for cc = 1:size(pMs,2)
        tcond = conditionsAndRasterTypes(rr,cc);
        nds = dec2base(tcond,10) - '0';
        varNames{rr,cc} = sprintf('C%dR%d',nds(1),nds(2));
    end
end
cN = 1;
%%
runthis = 1;
if runthis
    gAllVals = [];
    for si = 1:length(conditionsAndRasterTypes)
        tcond = conditionsAndRasterTypes(si);
        Ndigits = dec2base(tcond,10) - '0';
        for ani = 1:length(selAnimals)
            an = selAnimals(ani);
            data{ani,si} = squeeze(pMs{si}.all_zMIs{an}(Ndigits(1),Ndigits(2),:));
            gAllVals = [gAllVals;data{ani,si}];
        end
        legs{si} = sprintf('C-%d,RT-%d',Ndigits(1),Ndigits(2));
    end
    minBin = min(gAllVals);
    maxBin = 15;%max(gAllVals);
    incr = (maxBin-minBin)/100;
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 6 2 1.5],'color','w');
    hold on;
    [ha,hb,hca,sigR] = plotDistributions(data,'colors',colors,'maxY',90,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    hold on;
    ylim([0 120]);
    xlims = xlim; dx = xlims(2) - xlims(1); ylims = ylim; dy = ylims(2) - ylims(1);
    legs{length(legs)+1} = [xlims(1)+dx/2 dx/30 ylims(1)+dy/3 dy/15];
    putLegend(ha,legs,'colors',colors,'sigR',{sigR,'anova',sigColor,6});
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
    save_pdf(hf,mData.pdf_folder,sprintf('Mean zMI Condition %d',cN),600);
    
return;
end
