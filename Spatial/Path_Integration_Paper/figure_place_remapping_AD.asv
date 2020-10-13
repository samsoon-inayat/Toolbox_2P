function figure_place_remapping_10(fn,allRs,ccs)

protocol_C = '10_C';
protocol_A = '10_A';
ei_C = evalin('base','ei10_C');
ei_A = evalin('base','ei10_A');
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ET_C = evalin('base',sprintf('ET%s',protocol_C));
ET_A = evalin('base',sprintf('ET%s',protocol_A));
selAnimals_C = 1:length(ei_C)
selAnimals_A = 1:length(ei_A)

% in the following variable all the measurements are in the matrices form
% for each variable colums indicate raster and stim marker types specified 
% the rows indicate condition numbers.
paramMs_C = parameter_matrices('get','10_C');
paramMs_A = parameter_matrices('get','10_A');
% after getting all matrics, we can apply selection criteria to select a
% subgroup of cells
% here is the selection criteria in make_selC_structure function
% cellsOrNot = 1; planeNumber = NaN; zMI_Th =3; fwids = [1 120]; fcens = [0 140]; rs_th = 0.4;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
conditionsAndRasterTypes = [11 21 31 41];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN);
[cpMs_C,pMs_C] = parameter_matrices('select','10_C',{paramMs_C,selC});
[cpMs_A,pMs_A] = parameter_matrices('select','10_A',{paramMs_A,selC});
% perc_cells_C = parameter_matrices('print','10_C',{cpMs_C,pMs_C,ET_C,selAnimals_C});
% perc_cells_A = parameter_matrices('print','10_A',{cpMs_A,pMs_A,ET_A,selAnimals_A});

out_C = find_corr_coeff_rasters(pMs_C,paramMs_C,selAnimals_C,ei_C,conditionsAndRasterTypes,selC);
out_A = find_corr_coeff_rasters(pMs_A,paramMs_A,selAnimals_A,ei_A,conditionsAndRasterTypes,selC);


n = 0;
%%
if 1
    data_C_motion = out_C.mean_over_cells;
    data_A_motion = out_A.mean_over_cells;
    dataT_C = array2table(data_C_motion); dataT_C.Properties.VariableNames = {'CC1','CC2','CC3'};
    dataT_A = array2table(data_A_motion); dataT_A.Properties.VariableNames = {'CC1','CC2','CC3'};
    groupT = table([ones(size(dataT_C,1),1);2*ones(size(dataT_A,1),1)]); groupT.Properties.VariableNames = {'Group'};
    between = [groupT [dataT_C;dataT_A]]; between.Group = categorical(between.Group);
    within = table([1 2 3]');
    within.Properties.VariableNames = {'dCond'};
    within.dCond = categorical(within.dCond);
    rmaR = repeatedMeasuresAnova(between,within);
    writetable(between,fullfile(mData.pdf_folder,sprintf('%s_dataT.xlsx',mfilename)));
    writetable(rmaR.ranova,fullfile(mData.pdf_folder,sprintf('%s_dataT_stats_output.xlsx',mfilename)),'WriteRowNames',true);
    
%     findMeanAndStandardError
    mVar = rmaR.est_marginal_means.Mean;
    semVar = rmaR.est_marginal_means.Formula_StdErr;
    combs = rmaR.combs;
    p = rmaR.p; h = p<0.05;
    xdata = [1 2 3 4 5 6]; maxY = 0.3;
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 2.25 1],'color','w');
    hold on;
    tcolors = {colors{1};colors{2};colors{3};colors{1};colors{2};colors{3}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'maxY',maxY,'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 4:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+0],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'Ctrl-12','Ctrl-23','Ctrl-34','AD-12','AD-23','AD-34'};xtickangle(45)
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.02 0.03 0.02 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{'Correlation Coeff.',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bargraph_overall',mfilename),600);
end


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

function out = find_corr_coeff_rasters(pMs_C,paramMs_C,selAnimals_C,ei_C,conditionsAndRasterTypes,selC)

all_conds = []; all_rts = [];
for rr = 1:size(pMs_C,1)
    for cc = 1:size(pMs_C,2)
        tcond = conditionsAndRasterTypes(rr,cc);
        nds = dec2base(tcond,10) - '0';
        varNames{rr,cc} = sprintf('C%dR%d',nds(1),nds(2));
        all_conds = [all_conds nds(1)]; all_rts = [all_rts nds(2)];
        xticklabels{cc,rr} = sprintf('%s-%s',paramMs_C.stimMarkers{nds(2)},paramMs_C.rasterTypes{nds(2)}(1));
        for an = 1:length(selAnimals_C)
            tei = ei_C(an); conditionNumber = nds(1); rasterType = paramMs_C.rasterTypes{nds(2)}; 
            stimMarker = paramMs_C.stimMarkers{nds(2)}; maxDistTime = [140 Inf];%paramMs_C.maxDistTime;
            selCells = pMs_C{conditionNumber}.cellSel{an};
            cns = paramMs_C.all_cns{an};
            [temp_rasters ~] = getParamValues('rasters',tei,selC.plane_number,conditionNumber,stimMarker,rasterType,cns(selCells,2:3),maxDistTime);
            mean_rasters_C{an,cc} = squeeze(nanmean(temp_rasters,1))';
        end
    end
end

combs = [1 2;
        2 3;
        3 4];
all_ccs = cell(1);
for an = 1:length(selAnimals_C)
    ccs = [];
    for ii = 1:size(combs,1)
        cond1 = combs(ii,1); cond2 = combs(ii,2);
        set1 = mean_rasters_C{an,cond1}; set2 = mean_rasters_C{an,cond2};
        cc = [];
        for cn = 1:size(set1,1)
            cell_sig1 = set1(cn,:); cell_sig1 = fillmissing(cell_sig1,'linear',2,'EndValues','nearest');
            cell_sig2 = set2(cn,:); cell_sig2 = fillmissing(cell_sig2,'linear',2,'EndValues','nearest');
            temp = corrcoef(cell_sig1,cell_sig2);
            cc(cn) = temp(1,2);
        end
        ccs(ii,:) = cc;
    end
    all_ccs_C(an) = {ccs};
    mean_over_cells_C(an,:) = nanmean(ccs,2)';
end

out.mean_over_cells = mean_over_cells_C;