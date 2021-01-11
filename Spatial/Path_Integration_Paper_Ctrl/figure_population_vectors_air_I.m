function figure_population_vectors_air_I(fn,allRs,ccs)

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = 1.96; fwids = NaN; fcens = NaN; rs_th = NaN;
% cellsOrNot = 1; planeNumber = NaN; zMI_Th = 1.96; fwids = [1 150]; fcens = [0 150]; rs_th = 0.3;
conditionsAndRasterTypes = [15;25;35;45];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN);
out = read_data_from_base_workspace(selC)

ei_C = out.eis{1}; ei_A = out.eis{2};
pMs_C = out.pMs{1}; pMs_A = out.pMs{2};
paramMs_C = out.paramMs{1}; paramMs_A = out.paramMs{2};
cpMs_C = out.cpMs{1}; cpMs_A = out.cpMs{2};
selAnimals_C = out.selAnimals{1}; selAnimals_A = out.selAnimals{2};
perc_cells_C = out.perc_cells{1}; perc_cells_A = out.perc_cells{2};
% selAnimals_C = [1]; selAnimals_A = [3];
if 0
    cellSel_C = ''; cellSel_A = '';% cellSel_C = cpMs_C; cellSel_A = cpMs_A;
    a_trials = {3:9};
    for ii = 1:length(a_trials)
        out = get_mean_rasters(pMs_C',paramMs_C,selAnimals_C,ei_C,conditionsAndRasterTypes',selC,cellSel_C,a_trials{ii});
        [out.allP_an,out.allC_an,out.avg_C_conds,out.mean_rasters_T] = get_pop_vector_corr(out,conditionsAndRasterTypes,min(out.sz(:)),cellSel_C);
        all_out_C{ii} = out;
        out = get_mean_rasters(pMs_A',paramMs_A,selAnimals_A,ei_A,conditionsAndRasterTypes',selC,cellSel_A,a_trials{ii});
        [out.allP_an,out.allC_an,out.avg_C_conds,out.mean_rasters_T] = get_pop_vector_corr(out,conditionsAndRasterTypes,74,cellSel_A);
        all_out_A{ii} = out;
    end
    save('pop_corr_mat_vector_air_I.mat','all_out_C','all_out_A','a_trials');
    disp('Done');
    return;
else
    temp = load('pop_corr_mat_vector_air_I.mat');
    all_out_C = temp.all_out_C{1};     all_out_A = temp.all_out_A{1}; a_trials = temp.a_trials{1};
end
n = 0;
out_C = all_out_C; out_A = all_out_A;
%%

%%
if 0
    
    all_conds = []; all_rts = [];
    for rr = 1:size(pMs_C,1)
        for cc = 1:size(pMs_C,2)
            tcond = conditionsAndRasterTypes(rr,cc);
            nds = dec2base(tcond,10) - '0';
            varNames{rr,cc} = sprintf('C%dR%d',nds(1),nds(2));
            all_conds = [all_conds nds(1)]; all_rts = [all_rts nds(2)];
            xticklabels{cc,rr} = sprintf('%s-%s',paramMs_C.stimMarkers{nds(2)},paramMs_C.rasterTypes{nds(2)}(1));
        end
    end
    all_conds = unique(all_conds); all_rts = unique(all_rts);
    n = 0;
    
    numCols = length(all_rts);
    data = perc_cells_C;
    cmdTxt = sprintf('dataT_C = table(');
    for ii = 1:(size(data,2)-1)
        cmdTxt = sprintf('%sdata(:,%d),',cmdTxt,ii);
    end
    cmdTxt = sprintf('%sdata(:,size(data,2)));',cmdTxt);
    eval(cmdTxt);
    data = perc_cells_A;
    cmdTxt = sprintf('dataT_A = table(');
    for ii = 1:(size(data,2)-1)
        cmdTxt = sprintf('%sdata(:,%d),',cmdTxt,ii);
    end
    cmdTxt = sprintf('%sdata(:,size(data,2)));',cmdTxt);
    eval(cmdTxt);
    dataT = [dataT_C;dataT_A]
    dataT.Properties.VariableNames = varNames;
    dataT = [table([ones(length(ei_C),1);2*ones(length(ei_A),1)]) dataT];
    dataT.Properties.VariableNames{1} = 'Group';
    dataT.Group = categorical(dataT.Group)
    
    colVar1 = [ones(1,numCols) 2*ones(1,numCols) 3*ones(1,numCols) 4*ones(1,numCols)];    colVar2 = [1:numCols 1:numCols 1:numCols 1:numCols];
%     colVar1 = [1 2 3 4];
    within = table(colVar1');
%     within = table(colVar1',colVar2');
    within.Properties.VariableNames = {'Condition'};
%     within.Properties.VariableNames = {'Condition','Raster'};
    within.Condition = categorical(within.Condition);
%     within.Raster = categorical(within.Raster);
    ra = repeatedMeasuresAnova(dataT,within);

    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
     xdata = [1 2 3 4 6:9]; 
%     xdata = [1 2 3 4];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 2 1],'color','w');
    hold on;
    tcolors ={colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 5:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end

    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+1],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C1','C2','C3','C4','C1','C2','C3','C4'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30);
    changePosition(gca,[0.075 0.0 0 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Percentage of','tuned cells'},[0 -5 0]});
%     rectangle(gca,'Position',[0.75 maxY-20 1 19],'edgecolor','k','facecolor','k');     text(1.85,maxY-14+1,'CRTG','FontSize',6);
%     rectangle(gca,'Position',[6 maxY-20 1 19],'edgecolor','k');     text(7.2,maxY-14+1,'CPTG','FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('Percentage of Cs_I'),600);
return;
end

%%
if 0
    perc_cells_or_C = 100*cpMs_C.numCells./cpMs_C.areCells;
    perc_cells_or_A = 100*cpMs_A.numCells./cpMs_A.areCells;
    [h,p,ci,stats] = ttest2(perc_cells_or_C,perc_cells_or_A)
    
    mVar = [mean(perc_cells_or_C) mean(perc_cells_or_A)]; semVar = [std(perc_cells_or_C)/sqrt(3) std(perc_cells_or_A)/sqrt(5)];
    combs = [1 2]; %p = ra.mcs.p; h = ra.mcs.p < 0.05;
    xdata = [1:2];
%     xdata = [1 2 3 4];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors ={colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 5:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+1],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'CRTG','CPTG'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.17 0 -0.4 -0.09]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Overall percentage','of tuned cells'},[0 -5 0]});
    
    save_pdf(hf,mData.pdf_folder,sprintf('Percentage of Unique Cs_I'),600);
return;
end

%%
cccc = 1;
if cccc == 1
    sel_out = out_C;
    paramMs = paramMs_C;
    an = 2;
    ncols = min(sel_out.sz(:));
else
    sel_out = out_A;
    paramMs = paramMs_A;
    an = 1;
    ncols = min(sel_out.sz(:));
end
% [~,~,cellNums] = findPopulationVectorPlot(sel_out.mean_rasters{an,1},[]);
for ii = 1:length(conditionsAndRasterTypes)
    tcond = abs(conditionsAndRasterTypes(ii));
    Ndigits = dec2base(tcond,10) - '0';
    allP{ii} = sel_out.allP_an{an,ii};
    allC{ii} = sel_out.allC_an{an,ii};
%     mRsi = sel_out.mean_rasters{an,ii};
%     [allP{ii},allC{ii}] = findPopulationVectorPlot(mRsi(:,1:ncols),[]);
    time_xs{ii} = sel_out.xs{an,ii}(1:ncols);
    raster_labels{ii} = sprintf('Cond - %d, Rast - %d',Ndigits(1),Ndigits(2));
    raster_labels{ii} = sprintf('Condition - %d',Ndigits(1));
    theseRasterTypes{ii} = paramMs.rasterTypes{Ndigits(2)};
end
avg_C_conds = sel_out.avg_C_conds;
for ii = 1:length(conditionsAndRasterTypes)
    temp_c_m = nanmean(avg_C_conds{ii},3);
    avg_C_conds{ii} = temp_c_m;
%     avg_C_conds{ii} = imgaussfilt(temp_c_m,2);
end

avg_C_an = avg_C_conds;
n = 0;

%%
ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 4],...
    'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.1 0.1],'widthHeightAdjustment',...
    [0.01 -60]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[1 6 3.45 1.65]);
FS = mData.axes_font_size;
for sii = 1:length(conditionsAndRasterTypes)
    P = allP{sii};
    axes(ff.h_axes(1,sii));changePosition(gca,[0 0.05 -0.091 -0.1]);
    imagesc(P);
    box off;
    if sii == 1
        h = ylabel('Cell No.'); %   changePosition(h,[0 0 0]);
    end
    text(3,size(P,1)+round(size(P,1)/7),sprintf('%s',raster_labels{sii}),'FontSize',FS,'FontWeight','Normal');
    if size(P,1) > 1
        set(gca,'Ydir','Normal','linewidth',0.25,'FontSize',FS,'FontWeight','Bold','YTick',[1 size(P,1)]);
    else
        set(gca,'Ydir','Normal','linewidth',0.25,'FontSize',FS,'FontWeight','Bold','YTick',[1]);
    end
    cols = size(P,2);
    set(gca,'XTick',[]);
    
    %
    axes(ff.h_axes(2,sii));
    dec = -0.09;
    changePosition(gca,[0.0 0.05 dec dec]);
    imagesc(allC{sii},[-1 1]);
    minC(sii) = min(allC{sii}(:));
    maxC(sii) = max(allC{sii}(:));
    box off;
    set(gca,'Ydir','Normal','linewidth',1,'FontSize',FS,'FontWeight','Bold');
    if sii == 1
        h = ylabel('Time (sec)');    changePosition(h,[-5 0 0]);
    end
    cols = size(P,2);
    colsHalf = round(cols/2);
    ts = round(time_xs{sii}(1:cols));
    set(gca,'XTick',[]);
    if sii == 1
    set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    else
        set(gca,'YTick',[]);
    end
    h = xlabel('Time (sec)');%    changePosition(h,[0 0 0]);
    cols = size(P,2);
    colsHalf = ceil(cols/2);
    ts = round(time_xs{sii}(1:cols));
    set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    if sii == 1
        set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    else
        set(gca,'YTick',[]);
    end
end

colormap parula
mI = min(minC);
maxs = [1 1 1];
for ii = 1:4
    axes(ff.h_axes(1,ii)); caxis([0 maxs(1)]);
    axes(ff.h_axes(2,ii)); caxis([mI maxs(2)]);
%     axes(ff.h_axes(3,ii)); caxis([mIa maxs(3)]);
end

hc = putColorBar(ff.h_axes(1,4),[0.0 0.03 0 -0.05],[0 maxs(1)],6,'eastoutside',[0.07 0.07 0.1 0.1]);
hc = putColorBar(ff.h_axes(2,4),[0.0 0.03 0 -0.05],[mI maxs(2)],6,'eastoutside',[0.07 0.07 0.1 0.1]);
% hc = putColorBar(ff.h_axes(3,4),[0.0 0.03 0 -0.05],[mIa maxs(3)],6,'eastoutside',[0.07 0.07 0.1 0.1]);

save_pdf(ff.hf,mData.pdf_folder,sprintf('figure_place_cells_py_10_1_%d_I.pdf',cccc),600);

% return;
%%
ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 4],...
    'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.1 0.23],'widthHeightAdjustment',...
    [0.01 -300]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[1 6 3.45 0.95]);
FS = mData.axes_font_size;
for sii = 1:length(conditionsAndRasterTypes)
%     P = allP{sii};
    axes(ff.h_axes(1,sii));
    dec = -0.09;
    changePosition(gca,[0.0 0.05 dec dec]);
    imagesc(avg_C_an{sii});hold on;
    box off;
    if cccc == 1
        txty = 82;
    else
        txty = 82;
    end
    text(3,txty,sprintf('%s',raster_labels{sii}),'FontSize',FS,'FontWeight','Normal');
    min_avg_C_an(sii) = min(avg_C_an{sii}(:));
    max_avg_C_an(sii) = max(avg_C_an{sii}(:));
    box off;
    set(gca,'Ydir','Normal','linewidth',1,'FontSize',FS,'FontWeight','Bold');
    if sii == 1
        h = ylabel('Time (sec)');    changePosition(h,[-5 0 0]);
    end
    h = xlabel('Time (sec)');    changePosition(h,[0 0 0]);
    cols = size(P,2);
    colsHalf = ceil(cols/2);
    ts = round(time_xs{sii}(1:cols));
    set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    if sii == 1
    set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    else
        set(gca,'YTick',[]);
    end
end

colormap parula
mIa = min(min_avg_C_an);
maxs = [1 1 1];
for ii = 1:4
    axes(ff.h_axes(1,ii)); caxis([mIa maxs(1)]);
end

hc = putColorBar(ff.h_axes(1,4),[0.0 0.07 0 -0.09],[mIa maxs(1)],6,'eastoutside',[0.07 0.07 0.1 0.1]);

save_pdf(ff.hf,mData.pdf_folder,sprintf('figure_place_cells_py_10_2_%d_I.pdf',cccc),600);
