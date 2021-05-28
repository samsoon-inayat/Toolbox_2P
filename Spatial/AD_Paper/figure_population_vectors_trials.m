function figure_population_vectors_trials(fn,allRs,ccs)

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
% cellsOrNot = 1; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = 1.96; fwids = [1 150]; fcens = [0 150]; rs_th = 0.3; FR = NaN;
conditionsAndRasterTypes = [11;21;31;41];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN,FR);
out = read_data_from_base_workspace_AD(selC)

ei_C = out.eis{1}; ei_A = out.eis{2};
pMs_C = out.pMs{1}; pMs_A = out.pMs{2};
paramMs_C = out.paramMs{1}; paramMs_A = out.paramMs{2};
cpMs_C = out.cpMs{1}; cpMs_A = out.cpMs{2};
selAnimals_C = out.selAnimals{1}; selAnimals_A = out.selAnimals{2};
perc_cells_C = out.perc_cells{1}; perc_cells_A = out.perc_cells{2};

% trials = 7:10;
% out_C = get_mean_rasters(pMs_C',paramMs_C,selAnimals_C,ei_C,conditionsAndRasterTypes',selC,'',trials);
% out_A = get_mean_rasters(pMs_A',paramMs_A,selAnimals_A,ei_A,conditionsAndRasterTypes',selC,'',trials);
filename = fullfile(mData.pd_folder,'pop_corr_mat_trials_AD.mat');
if 0
    cellSel_C = ''; cellSel_A = ''; cellSel_C = cpMs_C; cellSel_A = cpMs_A;
    a_trials = {1:2,3:4,5:6,7:8,9:10};
    for ii = 1:length(a_trials)
        out = get_mean_rasters(pMs_C',paramMs_C,selAnimals_C,ei_C,conditionsAndRasterTypes',selC,cellSel_C,a_trials{ii});
        [out.allP_an,out.allC_an,out.avg_C_conds,out.mean_rasters_T,out.all_corr_an,out.all_corr_cell_an] = get_pop_vector_corr(out,conditionsAndRasterTypes,min(out.sz(:)),cellSel_C);
        all_out_C{ii} = out;
        out = get_mean_rasters(pMs_A',paramMs_A,selAnimals_A,ei_A,conditionsAndRasterTypes',selC,cellSel_A,a_trials{ii});
        [out.allP_an,out.allC_an,out.avg_C_conds,out.mean_rasters_T,out.all_corr_an,out.all_corr_cell_an] = get_pop_vector_corr(out,conditionsAndRasterTypes,49,cellSel_A);
        all_out_A{ii} = out;
    end
    save(filename,'all_out_C','all_out_A','a_trials');
    disp('Done');
    n= 0;
    return;
end
temp = load(filename);
all_out_C = temp.all_out_C;     all_out_A = temp.all_out_A; a_trials = temp.a_trials;

% filename = fullfile(mData.pd_folder,'pop_corr_mat_trials.mat');
if 0
pos_cell_corr_C = get_pop_vector_corr_trials(all_out_C);
pos_cell_corr_A = get_pop_vector_corr_trials(all_out_A);
save(filename,'all_out_C','all_out_A','a_trials','pos_cell_corr_C','pos_cell_corr_A');
disp('Done');
return;
end

pos_cell_corr_C = temp.pos_cell_corr_C;
pos_cell_corr_A = temp.pos_cell_corr_A;

% pause here then execute below
n = 0;
%%
all_CC_C = cell(5,5);
all_CC_A = cell(5,5);
mean_cell_corr_C = [];
mean_cell_corr_A = [];
group = 1; cn = 1;
for rr = 1:5
    for cc = 1:5
        thisCC = [];
        for an = 1:size(pos_cell_corr_C.pos_corr{1,1},1)
            thisCC(:,:,an) = pos_cell_corr_C.pos_corr{rr,cc}{an,cn};
            mean_cell_corr_C(rr,cc,an) = mean(diag(pos_cell_corr_C.cell_corr{rr,cc}{an,cn}));
        end
        all_CC_C{rr,cc} = nanmean(thisCC,3);
        
        thisCC = [];
        for an = 1:size(pos_cell_corr_A.pos_corr{1,1},1)
            thisCC(:,:,an) = pos_cell_corr_A.pos_corr{rr,cc}{an,cn};
            mean_cell_corr_A(rr,cc,an) = mean(diag(pos_cell_corr_A.cell_corr{rr,cc}{an,cn}));
        end
        all_CC_A{rr,cc} = nanmean(thisCC,3);
    end
end

if group == 1
    sel_out = all_CC_C;
    sel_p = all_out_C;
%     an = 1; cn = 3;
else
    sel_out = all_CC_A;
    sel_p = all_out_A;
end

ff = makeFigureRowsCols(107,[1 0.5 2 2],'RowsCols',[5 5],...
    'spaceRowsCols',[0.01 0.01],'rightUpShifts',[0.15 0.12],'widthHeightAdjustment',...
    [-45 -45]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[1 6 1.7 1.7]);
FS = mData.axes_font_size;
maskdisp = triu(ones(5,5),0);

for rr = 1:5
    for cc = 1:5
        if ~maskdisp(rr,cc)
            delete(ff.h_axes(rr,cc));
            continue;
        end
        axes(ff.h_axes(rr,cc));
        this_pos_corr = sel_out{rr,cc};
        imagesc(this_pos_corr);
        box off;
        set(gca,'Ydir','Normal','linewidth',1,'FontSize',FS-1,'FontWeight','Bold');
        if rr == 1 && cc == 1
            cols = size(this_pos_corr,2);
            colsHalf = round(cols/2);
            ts = round(sel_p{cc}.xs{3,1}(1:cols));
            set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1)-2 ts(colsHalf) ts(cols)+2]);
            h = ylabel('Pos (cm)');%    changePosition(h,[0 0 0]);
            set(gca,'XTick',[]);
        elseif rr == 5 && cc == 5
            h = xlabel('Pos (cm)');%    changePosition(h,[0 0 0]);
            cols = size(this_pos_corr,2);
            colsHalf = round(cols/2);
            ts = round(sel_p{cc}.xs{3,1}(1:cols));
            set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1)-2 ts(colsHalf) ts(cols)+2]);
            set(gca,'YTick',[]);
        else
            axis off
        end
        minC = min(this_pos_corr(:));
        maxC = max(this_pos_corr(:));
%         text(5,cols-5,sprintf('(%.1f, %.1f)',minC,maxC),'FontSize',5,'Color','w');
%         hc = putColorBar(ff.h_axes(rr,cc),[0.0 0.03 0 -0.05],[minC maxC],5,'eastoutside',[0.07 0.11 0.1 0.16]);
    end
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('pos corr trials_group_%d_cond_%d',group,cn),600);



    temp = mean_cell_corr_C(:,:,1);
    mask = ones(size(temp)); mask = triu(mask,1) & ~triu(mask,2);
    for ii = 1:size(mean_cell_corr_C,3)
        temp = mean_cell_corr_C(:,:,ii);
        var_C(ii,:) = temp(mask)';
    end
    for ii = 1:size(mean_cell_corr_A,3)
        temp = mean_cell_corr_A(:,:,ii);
        var_A(ii,:) = temp(mask)';
    end
    dataT = array2table([[ones(size(var_C,1),1);2*ones(size(var_A,1),1)] [var_C;var_A]]);
    dataT.Properties.VariableNames = {'Group','T1234','T3456','T5678','T78910'};
    dataT.Group = categorical(dataT.Group);
    colVar1 = [1 2 3 4];    
    within = table(colVar1');
    within.Properties.VariableNames = {'Condition'};
    within.Condition = categorical(within.Condition);
    ra = repeatedMeasuresAnova(dataT,within);
    
    dataTC = dataT(1:3,2:end);
    raC = repeatedMeasuresAnova(dataTC,within);
    
    dataTA = dataT(4:end,2:end);
    raA = repeatedMeasuresAnova(dataTA,within);

    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
    xdata = [1 2 3 4 6:9];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.75 1],'color','w');
    hold on;
    tcolors ={colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 5:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'TC1','TC2','TC3','TC4','TC1','TC2','TC3','TC4'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.06 0.03 0.02 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Correlation'},[0 0 0]});

    save_pdf(hf,mData.pdf_folder,sprintf('cell corr remap trials _ group_%d_cond_%d',group,cn),600);

% return;
%%
n = 0;
out_C = all_out_C{4};
out_A = all_out_A{3};

cccc = 1;
if cccc == 1
    sel_out = out_C;
    paramMs = paramMs_C;
    an = 3;
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
    mRsi = sel_out.mean_rasters{an,ii};
    [allP{ii},allC{ii}] = findPopulationVectorPlot(mRsi(:,1:ncols),[]);
    time_xs{ii} = sel_out.xs{an,ii}(1:ncols);
    raster_labels{ii} = sprintf('Cond - %d, Rast - %d',Ndigits(1),Ndigits(2));
    raster_labels{ii} = sprintf('Condition - %d',Ndigits(1));
    theseRasterTypes{ii} = paramMs.rasterTypes{Ndigits(2)};
end

commonZ = zeros(ncols,ncols);
avg_C_an = repmat(commonZ,1,1,size(sel_out.sz,1));
for ii = 1:length(conditionsAndRasterTypes)
    avg_C_conds{ii} = avg_C_an;
end


for an = 1:size(sel_out.sz,1)
    for ii = 1:length(conditionsAndRasterTypes)
        tcond = abs(conditionsAndRasterTypes(ii));
        Ndigits = dec2base(tcond,10) - '0';
        mRsi = sel_out.mean_rasters{an,ii};
        if size(mRsi,2) < ncols
            cncols = size(mRsi,2);
            mRsi(:,(cncols+1:ncols)) = nan(size(mRsi,1),length(cncols+1:ncols));
        end
        [allP_an{an,ii},allC_an{an,ii}] = findPopulationVectorPlot(mRsi(:,1:ncols),[]);
        avg_C_conds{ii}(:,:,an) = allC_an{an,ii};
    end
end
% save('pop_corr_A.mat','avg_C_conds');
% save('pop_corr_C.mat','avg_C_conds');
for ii = 1:length(conditionsAndRasterTypes)
    temp_c_m = nanmean(avg_C_conds{ii},3);
    avg_C_conds{ii} = temp_c_m;
%     avg_C_conds{ii} = imgaussfilt(temp_c_m,2);
end
avg_C_an = avg_C_conds;
n = 0;

%%
if 1
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
return
end

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
%     text(3,txty,sprintf('%s',raster_labels{sii}),'FontSize',FS,'FontWeight','Normal');
    min_avg_C_an(sii) = min(avg_C_an{sii}(:));
    max_avg_C_an(sii) = max(avg_C_an{sii}(:));
    box off;
    set(gca,'Ydir','Normal','linewidth',1,'FontSize',FS,'FontWeight','Bold');
    if sii == 1
        h = ylabel('Position (cm)');    changePosition(h,[-5 0 0]);
    end
    h = xlabel('Position (cm)');    changePosition(h,[0 0 0]);
    cols = size(avg_C_an{sii},2);
    colsHalf = ceil(cols/2);
    ts = round(time_xs{sii}(1:cols));
    set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1)-2 ts(colsHalf) ts(cols)+2]);
    if sii == 1
        set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1)-2 ts(colsHalf) ts(cols)+2]);
    else
        set(gca,'YTick',[]);
    end
%     set(gca,'XTick',[]);
end

colormap parula
mIa = min(min_avg_C_an);
maxs = [1 1 1];
for ii = 1:4
    axes(ff.h_axes(1,ii)); caxis([mIa maxs(1)]);
end

hc = putColorBar(ff.h_axes(1,4),[0.0 0.07 0 -0.09],[mIa maxs(1)],6,'eastoutside',[0.07 0.07 0.1 0.1]);

save_pdf(ff.hf,mData.pdf_folder,sprintf('figure_place_cells_py_10_2_%d.pdf',cccc),600);
