function figure_place_remapping(fn,allRs,ccs)

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
cellsOrNot = NaN; planeNumber = NaN; zMI_Th = 3; fwids = [1 120]; fcens = [0 140]; rs_th = 0.4; HaFD_th = NaN; HiFD_th = NaN;
% cellsOrNot = NaN; planeNumber = NaN; zMI_Th = 3; fwids = NaN; fcens = NaN; rs_th = 0.4;
conditionsAndRasterTypes = [21 31];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,HaFD_th,HiFD_th);
[cpMs,pMs] = parameter_matrices('select',protocol,{paramMs,selC});
perc_cells = parameter_matrices('print',protocol,{cpMs,pMs,ET,selAnimals});


%% place field centers scatter plot from one condition to another
runthis = 1;
if runthis
    for si = 1:length(conditionsAndRasterTypes)
        tcond = conditionsAndRasterTypes(si);
        Ndigits = dec2base(tcond,10) - '0';
        pcs = [];
        for ani = 1:length(selAnimals)
            an = selAnimals(ani);
            thisan_pcs = squeeze(cpMs.all_fcenters{an}(Ndigits(1),Ndigits(2),:));
            pcs = [pcs;thisan_pcs];
            apcs_an{ani,si} = thisan_pcs;
        end
        apcs{si} = pcs;
    end
    apcs_an
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 6 4 3],'color','w');
    hold on;
    xvals = apcs{1}; yvals = apcs{2}; ii = 1;
    scatter(xvals,yvals,20,'.','MarkerEdgeColor',colors{ii},'MarkerFaceColor','none');
    mdl = fit_linear(xvals,yvals);
    co = mdl.Coefficients{:,1};
    rsq(ii) = mdl.Rsquared.Ordinary;
    yhvals = feval(mdl,xvals);
%     ft = fittype('poly1');
%     [ftc,gof,output] = fit(xvals,yvals,ft);   rsq(ii) = gof.rsquare;
%     co = coeffvalues(ftc)
%     yhvals = ft(co(1),co(2),xvals);
    hold on;plot(xvals,yhvals,'color',colors{ii},'LineWidth',1)
    text(20,15-3*ii,sprintf('Rsq = %.3f',rsq(ii)),'color',colors{ii});
  

    data{1} = xvals-yvals;
    gAllVals = data{1};
    minBin = min(gAllVals);
    maxBin = max(gAllVals);
    incr = (maxBin-minBin)/50;
    hf = figure(1001);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 8 2 1.5],'color','w');
    hold on;
    [ha,hb,hca] = plotDistributions(data,'colors',colors,'maxY',10,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    ds = descriptiveStatistics(data{1},'decimal_places',1)
    return;
end



%%
runthis = 0;
if runthis
trials = 3:10;
trials10 = 3:9;
% align cells
stimMarkers = paramMs.stimMarkers;
rasterTypes = paramMs.rasterTypes;
CNi = 1; 
conditionNumber = 1; rasterTypeN = 1;

if 0
    for si = 1:length(stimMarkers)
        stimMarker = stimMarkers{si};
        rasterType = rasterTypes{si};
        mRsi = [];
        for ani = 1:length(selAnimals)
            an = selAnimals(ani);
            tei = ei(an);
            selCells = pMs.cellSel{an};
            cns = paramMs.all_cns{an};
            maxDistTime = paramMs.maxDistTime;
            [tempD cnso] = getParamValues('rasters',tei,selC.plane_number,conditionNumber,stimMarker,rasterType,cns(selCells,2:3),maxDistTime);
            if length(tempD) == 0
                continue;
            end
            try
                 mR = findMeanRasters(tempD,trials);
            catch
                 mR = findMeanRasters(tempD,trials10);
            end
            mRsi = [mRsi;mR];
        end
        [temp,~,~] = getParamValues('',ei(1),1,1,stimMarker,rasterType,'areCells',[Inf Inf]);
        dxs = diff(temp.xs); bin_width = dxs(1); xs = 0:bin_width:1000;
        allRs{si} = mRsi;
        time_xs{si} = xs(1:size(mRsi,2));
        raster_labels{si} = sprintf('%s-%s',stimMarker,rasterType);
    end
    [~,~,cellNums] = findPopulationVectorPlot(allRs{CNi},[]);
    for ii = 1:length(stimMarkers)
        mRsi = allRs{ii};
        [allP{ii},allC{ii}] = findPopulationVectorPlot(mRsi,[],cellNums);
    end
else
    for si = 1:4
        stimMarker = stimMarkers{rasterTypeN};
        rasterType = rasterTypes{rasterTypeN};
        mRsi = [];
        for ani = 1:length(selAnimals)
            an = selAnimals(ani);
            tei = ei(an);
            selCells = apMs{si}.cellSel{an};
            cns = paramMs.all_cns{an};
            maxDistTime = paramMs.maxDistTime;
            [tempD cnso] = getParamValues('rasters',tei,selC.plane_number,si,stimMarker,rasterType,cns(selCells,2:3),maxDistTime);
            if length(tempD) == 0
                continue;
            end
            try
                 mR = findMeanRasters(tempD,trials);
            catch
                 mR = findMeanRasters(tempD,trials10);
            end
            mRsi = [mRsi;mR];
        end
        [temp,~,~] = getParamValues('',ei(1),1,1,stimMarker,rasterType,'areCells',[Inf Inf]);
        dxs = diff(temp.xs); bin_width = dxs(1); xs = 0:bin_width:1000;
        allRs{si} = mRsi;
        time_xs{si} = xs(1:size(mRsi,2));
        raster_labels{si} = sprintf('Condition # %d',si);
    end
    [~,~,cellNums] = findPopulationVectorPlot(allRs{CNi},[]);
    for ii = 1:length(stimMarkers)
        mRsi = allRs{ii};
        [allP{ii},allC{ii}] = findPopulationVectorPlot(mRsi,[]);
    end
end





n = 0;

%%
ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 4],...
    'spaceRowsCols',[-0.01 -0.05],'rightUpShifts',[0.1 0.06],'widthHeightAdjustment',...
    [15 -30]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[1 6 3.5 2]);
FS = 4;
for sii = 1:length(stimMarkers)
    P = allP{sii};
    axes(ff.h_axes(1,sii));changePosition(gca,[0 0.05 -0.091 -0.1]);
    imagesc(P);
    box off;
%     axis off
%     if sii == 1
%         text(-5,1,'1','FontSize',FS,'FontWeight','Normal');
%         if size(P,1) < 100
%             text(-7,size(P,1),sprintf('%d',size(P,1)),'FontSize',FS,'FontWeight','Normal');
%         else
%             text(-10,size(P,1),sprintf('%d',size(P,1)),'FontSize',FS,'FontWeight','Normal');
%         end
        if sii == 1
            text(-21,25,sprintf('Cells'),'FontSize',FS+3,'FontWeight','Bold','rotation',90);
        end
%     end
    text(3,size(P,1)+round(size(P,1)/10),sprintf('%s',raster_labels{sii}),'FontSize',FS,'FontWeight','Normal');
    set(gca,'Ydir','Normal','linewidth',0.25,'FontSize',FS,'FontWeight','Bold','YTick',[1 size(P,1)]);
    cols = size(P,2);
    colsHalf = ceil(cols/2);
    ts = round(time_xs{sii}(1:cols));
%     set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    set(gca,'XTick',[]);
    axes(ff.h_axes(2,sii));
    dec = -0.09;
    changePosition(gca,[0.0 0.05 dec dec]);
    imagesc(allC{sii},[-1 1]);
    minC(sii) = min(allC{sii}(:));
    box off;
%     axis equal
%     axis off
%     if sii == 1        
%         text(-8,3,'0','FontSize',FS,'FontWeight','Normal');
%         text(-10,50,num2str(mData.belt_length),'FontSize',FS,'FontWeight','Normal');
%     end
%     text(-1,-3,'0','FontSize',FS,'FontWeight','Normal');
%     text(44,-3,num2str(mData.belt_length),'FontSize',FS,'FontWeight','Normal');
%     if sii == 2
%         text(35,-13,sprintf('Position (cm)'),'FontSize',FS+3,'FontWeight','Bold','rotation',0);
%     end
%     if sii == 1
%         text(-21,-3,sprintf('Position (cm)'),'FontSize',FS+3,'FontWeight','Bold','rotation',90);
%     end
    set(gca,'Ydir','Normal','linewidth',1,'FontSize',FS,'FontWeight','Bold');
    if strcmp(rasterTypes{sii},'dist')
        h = xlabel('Position (cm)');    changePosition(h,[0 0 0]);
    else
        h = xlabel('Time (sec)');    changePosition(h,[0 0 0]);
    end
    if sii == 1
%         h = ylabel('Position (cm)');    changePosition(h,[1 0 0]);
    end
    cols = size(P,2);
    colsHalf = ceil(cols/2);
    ts = round(time_xs{sii}(1:cols));
    set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
end

colormap parula
mI = min(minC);
for ii = 1:4
    axes(ff.h_axes(2,ii));
    caxis([mI 1]);
end

axes(ff.h_axes(2,4));
pos = get(gca,'Position');
h = axes('Position',pos);
changePosition(h,[0.14 0.06 0.0 -0.1]);
hc = colorbar; axis off
set(hc,'YTick',[],'linewidth',0.01);
changePosition(hc,[0 0 0 -0.01]);
ylims = ylim;
xlims = xlim;
text(xlims(1)+0.33,ylims(1)-0.05,sprintf('%.2f',mI),'FontSize',4.5);
text(xlims(1)+0.39,ylims(2)+0.045,'1','FontSize',4.5);

axes(ff.h_axes(1,4));
pos = get(gca,'Position');
h = axes('Position',pos);
changePosition(h,[0.14 0.06 0.0 -0.1]);
hc = colorbar; axis off
set(hc,'YTick',[],'linewidth',0.01);
changePosition(hc,[0 0 0 -0.01]);
ylims = ylim;
xlims = xlim;
text(xlims(1)+0.4,ylims(1)-0.05,sprintf('0'),'FontSize',4.5);
text(xlims(1)+0.39,ylims(2)+0.045,'1','FontSize',4.5);

% delete(ff.h_axes(1,4));
% delete(ff.h_axes(2,4));
fileName = fullfile(mData.pdf_folder,'figure_place_cells_py_10.pdf');
save_pdf(ff.hf,mData.pdf_folder,'figure_place_cells_py_10.pdf',600);
close(ff.hf);
winopen(fileName);
return;
%%
% bar graph from anova
sigR = significanceTesting(allP);
ff = makeFigureWindow__one_axes_only(3,[2 4 2 2],[0.3 0.22 0.68 0.72]);
set(gcf,'color','w');
set(gcf,'Position',[3 4 1.2 1.5]);
axes(ff.ha); hs = sigR.anova.multcompare.h; ps = sigR.anova.multcompare.p;
plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[hs ps],'colors',colors,'sigColor',sigColor,'maxY',0.3,'ySpacing',0.03,'sigTestName','ANOVA');
% plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[hs ps],'colors',colors,'sigColor',sigColor,'maxY',0.57,'ySpacing',0.04,'sigTestName','ANOVA');
xlim([0.4 0.6+length(sigR.means)]);
%     ylim([0 max(mVals)]);
hyl = ylabel('Average Normalized FR');
pos = get(hyl,'Position');pos = pos + [+0.4 0 0];set(hyl,'Position',pos);
% set(ff.ha,'linewidth',1);
set(ff.ha,'TickDir','out','FontSize',7,'FontWeight','Normal');
set(ff.ha,'XTick',[1 2 3 4],'XTickLabel',{'Context 1','Context 2','Context 3','Context4'});
xtickangle(25);
save2pdf('figure_place_cells_bars.pdf',ff.hf,600);
return;
end