function figure_population_vectors(fn,allRs,ccs)

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN; FR = NaN;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = 1.96; fwids = [1 150]; fcens = [0 150]; rs_th = 0.3; FR = NaN;
conditionsAndRasterTypes = [11 21 31 41];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN,FR);
out = read_data_from_base_workspace(selC)

ei_C = out.eis{1}; ei_A = out.eis{2};
pMs_C = out.pMs{1}; pMs_A = out.pMs{2};
paramMs_C = out.paramMs{1}; paramMs_A = out.paramMs{2};
cpMs_C = out.cpMs{1}; cpMs_A = out.cpMs{2};
selAnimals_C = out.selAnimals{1}; selAnimals_A = out.selAnimals{2};
perc_cells_C = out.perc_cells{1}; perc_cells_A = out.perc_cells{2};

if 0
    cellSel_C = ''; cellSel_A = ''; %cellSel_C = cpMs_C; cellSel_A = cpMs_A;
    a_trials = {3:10};
    for ii = 1:length(a_trials)
        out = get_mean_rasters(pMs_C',paramMs_C,selAnimals_C,ei_C,conditionsAndRasterTypes',selC,cellSel_C,a_trials{ii});
        [out.allP_an,out.allC_an,out.avg_C_conds,out.mean_rasters_T] = get_pop_vector_corr(out,conditionsAndRasterTypes,min(out.sz(:)));
        all_out_C{ii} = out;
        out = get_mean_rasters(pMs_A',paramMs_A,selAnimals_A,ei_A,conditionsAndRasterTypes',selC,cellSel_A,a_trials{ii});
        [out.allP_an,out.allC_an,out.avg_C_conds,out.mean_rasters_T] = get_pop_vector_corr(out,conditionsAndRasterTypes,49);
        all_out_A{ii} = out;
    end
    save('pop_corr_mat_vector.mat','all_out_C','all_out_A','a_trials');
    disp('Done');
    return;
else
    temp = load('pop_corr_mat_vector.mat');
    all_out_C = temp.all_out_C;     all_out_A = temp.all_out_A; a_trials = temp.a_trials{1};
end
n = 0;
out_C = all_out_C{1}; out_A = all_out_A{1};
%%
cccc = 1;
if cccc == 1
    sel_out = out_C;
    paramMs = paramMs_C;
    an = 1;
    ncols = min(sel_out.sz(:));
else
    sel_out = out_A;
    paramMs = paramMs_A;
    an = 4;
    ncols = 49;
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
        h = ylabel('Position (cm)');    changePosition(h,[-5 0 0]);
    end
    cols = size(P,2);
    colsHalf = round(cols/2);
    ts = round(time_xs{sii}(1:cols));
    set(gca,'XTick',[]);
    if sii == 1
    set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1)-2 ts(colsHalf) ts(cols)+2]);
    else
        set(gca,'YTick',[]);
    end
    h = xlabel('Position (cm)');%    changePosition(h,[0 0 0]);
    cols = size(P,2);
    colsHalf = ceil(cols/2);
    ts = round(time_xs{sii}(1:cols));
    set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1)-2 ts(colsHalf) ts(cols)+2]);
    if sii == 1
        set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1)-2 ts(colsHalf) ts(cols)+2]);
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

save_pdf(ff.hf,mData.pdf_folder,sprintf('figure_place_cells_py_10_1_%d.pdf',cccc),600);


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
        txty = 53;
    else
        txty = 55;
    end
    text(3,txty,sprintf('%s',raster_labels{sii}),'FontSize',FS,'FontWeight','Normal');
    min_avg_C_an(sii) = min(avg_C_an{sii}(:));
    max_avg_C_an(sii) = max(avg_C_an{sii}(:));
    box off;
    set(gca,'Ydir','Normal','linewidth',1,'FontSize',FS,'FontWeight','Bold');
    if sii == 1
        h = ylabel('Position (cm)');    changePosition(h,[-5 0 0]);
    end
    h = xlabel('Position (cm)');    changePosition(h,[0 0 0]);
    cols = size(P,2);
    colsHalf = ceil(cols/2);
    ts = round(time_xs{sii}(1:cols));
    set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1)-2 ts(colsHalf) ts(cols)+2]);
    if sii == 1
    set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1)-2 ts(colsHalf) ts(cols)+2]);
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

save_pdf(ff.hf,mData.pdf_folder,sprintf('figure_place_cells_py_10_2_%d.pdf',cccc),600);
