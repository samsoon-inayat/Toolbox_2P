function figure_population_vectors(fn,allRs,ccs)

protocol = '10';
% protocol = '15';
ei = evalin('base',sprintf('ei%s',protocol));
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ET = evalin('base',sprintf('ET%s',protocol));
selAnimalsA = eval(sprintf('mData.selAnimals%s',protocol));
selAnimals = selAnimalsA;

% in the following variable all the measurements are in the matrices form
% for each variable colums indicate raster and stim marker types specified 
% the rows indicate condition numbers.
paramMs = parameter_matrices('get',protocol);
% after getting all matrics, we can apply selection criteria to select a
% subgroup of cells
% here is the selection criteria in make_selC_structure function
% cellsOrNot = NaN; planeNumber = NaN; zMI_Th = 3; fwids = [1 120]; fcens = [0 140]; rs_th = 0.4;
cellsOrNot = NaN; planeNumber = NaN; zMI_Th = 5; fwids = NaN; fcens = NaN; rs_th = NaN; HaFD_th = NaN; HiFD_th = NaN;
conditionsAndRasterTypes = [11 12 13 14 15 21 22 23 24 25 31 32 33 34 35 41 42 43 44 45]';
% conditionsAndRasterTypes = [11 13 21 23 31 33 41 43]';
conditionsAndRasterTypes = [11 21 31 41];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,HaFD_th,HiFD_th);
[cpMs,pMs] = parameter_matrices('select',protocol,{paramMs,selC});
% perc_cells = parameter_matrices('print',protocol,{cpMs,pMs,ET,selAnimals});

for rr = 1:size(pMs,1)
    for cc = 1:size(pMs,2)
        tcond = conditionsAndRasterTypes(rr,cc);
        nds = dec2base(tcond,10) - '0';
        varNames{rr,cc} = sprintf('C%dR%d',nds(1),nds(2));
        for an = 1:length(selAnimals)
            zMIs(an,rr,cc) = nanmean(squeeze(pMs{rr,cc}.all_zMIs{selAnimals(an)}(nds(1),nds(2),:)));
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
trials = 3:10;
trials10 = 3:9;
% align cells
stimMarkers = paramMs.stimMarkers;
rasterTypes = paramMs.rasterTypes;
CNi = 3;
rasterTypeN = 1;


for si = 1:length(conditionsAndRasterTypes)
    tcond = abs(conditionsAndRasterTypes(si));
    Ndigits = dec2base(tcond,10) - '0';
    theseRasterTypes{si} = rasterTypes{Ndigits(2)};
    mRsi = [];
    for ani = 1:length(selAnimals)
%         [si ani]
        an = selAnimals(ani);
        tei = ei(an);
        selCells = pMs{si}.cellSel{an};
%         selCells = cpMs.cellSel{an};
        cns = paramMs.all_cns{an};
        maxDistTime = paramMs.maxDistTime;
        [tempD cnso] = getParamValues('rasters',tei,selC.plane_number,Ndigits(1),stimMarkers{Ndigits(2)},rasterTypes{Ndigits(2)},...
            cns(selCells,2:3),maxDistTime);
        if length(tempD) == 0
            continue;
        end
        try
             mR = findMeanRasters(tempD,trials);
        catch
             mR = findMeanRasters(tempD,trials10);
        end
        if ~isempty(mRsi)
            try
                mR = mR(:,1:size(mRsi,2));
            catch
                mRsi = mRsi(:,1:size(mR,2));
            end
        end
        mRsi = [mRsi;mR];
    end
    [temp,~,~] = getParamValues('',ei(1),1,Ndigits(1),stimMarkers{Ndigits(2)},rasterTypes{Ndigits(2)},'areCells',[Inf Inf]);
    dxs = diff(temp.xs); bin_width = dxs(1); xs = 0:bin_width:1000;
    allRs{si} = mRsi;
    time_xs{si} = xs(1:size(mRsi,2));
    raster_labels{si} = sprintf('Cond - %d, Rast - %d',Ndigits(1),Ndigits(2));
end

% [~,~,cellNums] = findPopulationVectorPlot(allRs{CNi},[]);
for ii = 1:length(conditionsAndRasterTypes)
    mRsi = allRs{ii};
    [allP{ii},allC{ii}] = findPopulationVectorPlot(mRsi,[]);
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
for sii = 1:length(conditionsAndRasterTypes)
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
    if size(P,1) > 1
        set(gca,'Ydir','Normal','linewidth',0.25,'FontSize',FS,'FontWeight','Bold','YTick',[1 size(P,1)]);
    else
        set(gca,'Ydir','Normal','linewidth',0.25,'FontSize',FS,'FontWeight','Bold','YTick',[1]);
    end
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
    if strcmp(theseRasterTypes{sii},'dist')
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