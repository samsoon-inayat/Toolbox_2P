function figure_place_remapping_AD(fn,allRs,ccs)

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = 1.96; fwids = [1 150]; fcens = [0 150]; rs_th = 0.3;
conditionsAndRasterTypes = [11 21 31 41];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN);
out = read_data_from_base_workspace(selC)

ei_C = out.eis{1}; ei_A = out.eis{2};
pMs_C = out.pMs{1}; pMs_A = out.pMs{2};
paramMs_C = out.paramMs{1}; paramMs_A = out.paramMs{2};
cpMs_C = out.cpMs{1}; cpMs_A = out.cpMs{2};
selAnimals_C = out.selAnimals{1}; selAnimals_A = out.selAnimals{2};
perc_cells_C = out.perc_cells{1}; perc_cells_A = out.perc_cells{2};

andor = 1;
if andor == 1
    selAnimals_A = [3 4 5];
end

group = 1; 


if 0
    cellSel_C = ''; cellSel_A = ''; cellSel_C = cpMs_C; cellSel_A = cpMs_A;
    a_trials = {3:10};
    for ii = 1:length(a_trials)
        out = get_mean_rasters(pMs_C',paramMs_C,selAnimals_C,ei_C,conditionsAndRasterTypes',selC,cellSel_C,a_trials{ii});
        [out.allP_an,out.allC_an,out.avg_C_conds,out.mean_rasters_T,out.all_corr_an,out.all_corr_cell_an] = get_pop_vector_corr(out,conditionsAndRasterTypes,min(out.sz(:)));
        all_out_C{ii} = out;
        out = get_mean_rasters(pMs_A',paramMs_A,selAnimals_A,ei_A,conditionsAndRasterTypes',selC,cellSel_A,a_trials{ii});
        [out.allP_an,out.allC_an,out.avg_C_conds,out.mean_rasters_T,out.all_corr_an,out.all_corr_cell_an] = get_pop_vector_corr(out,conditionsAndRasterTypes,49);
        all_out_A{ii} = out;
    end
    if andor == 1
        save('pop_corr_mat_remapping_and.mat','all_out_C','all_out_A','a_trials');
    else
        save('pop_corr_mat_remapping_or.mat','all_out_C','all_out_A','a_trials');
    end
    disp('Done');
    return;
else
    temp = load('pop_corr_mat_remapping_and.mat');
    all_out_C_and = temp.all_out_C{1};     all_out_A_and = temp.all_out_A{1}; a_trials = temp.a_trials{1};
    temp = load('pop_corr_mat_remapping_or.mat');
    all_out_C_or = temp.all_out_C{1};     all_out_A_or = temp.all_out_A{1}; a_trials = temp.a_trials{1};
end

if andor == 1
    all_out_C = all_out_C_and;
    all_out_A = all_out_A_and;
else
    all_out_C = all_out_C_or;
    all_out_A = all_out_A_or;
end

all_CC_C = cell(4,4);
all_CC_A = cell(4,4);
for rr = 1:4
    for cc = 1:4
        thisCC = [];
        for an = 1:length(all_out_C.all_corr_an)
            thisCC(:,:,an) = all_out_C.all_corr_an{an}{rr,cc};
            mean_cell_corr_C(rr,cc,an) = mean(diag(all_out_C.all_corr_cell_an{an}{rr,cc}));
        end
        all_CC_C{rr,cc} = nanmean(thisCC,3);
        
        thisCC = [];
        for an = 1:length(all_out_A.all_corr_an)
            thisCC(:,:,an) = all_out_A.all_corr_an{an}{rr,cc};
            mean_cell_corr_A(rr,cc,an) = mean(diag(all_out_A.all_corr_cell_an{an}{rr,cc}));
        end
        all_CC_A{rr,cc} = nanmean(thisCC,3);
    end
end


n = 0;
%%
if 0
    n = 0;
    perc_cells_common_C = 100*(cpMs_C.numCells./cpMs_C.areCells);
    perc_cells_common_A = 100*(cpMs_A.numCells./cpMs_A.areCells);
    [h,p,ci,stats] = ttest2(perc_cells_common_C,perc_cells_common_A);
    mVar = [mean(perc_cells_common_C) mean(perc_cells_common_A)]; 
    semVar = [std(perc_cells_common_C)./sqrt(length(perc_cells_common_C)) std(perc_cells_common_A)./sqrt(length(perc_cells_common_A))];
    combs = [1 2];
    xdata = [1 2];
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
    changePosition(gca,[0.17 0.03 -0.3 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Common Spatially','Tuned Cells (%)'},[0 0 0]});

    save_pdf(hf,mData.pdf_folder,'Perc Of Place Cells_common',600);
    
    return;
end


%%
if 1
% out_CC = all_out_C.all_corr_an{3};
% out_AA = all_out_A.all_corr_an{1};
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
    dataT.Properties.VariableNames = {'Group','C12','C23','C34'};
    dataT.Group = categorical(dataT.Group);
    colVar1 = [1 2 3];    
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
    xdata = [1 2 3 5:7];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    tcolors ={colors{1};colors{2};colors{3};colors{1};colors{2};colors{3};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = 4:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C12','C23','C34','C12','C23','C34'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.06 0.03 0.02 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Correlation'},[0 0 0]});

    save_pdf(hf,mData.pdf_folder,'cell corr remap',600);
return;
end


%%
% out_CC = all_out_C.all_corr_an{3};
% out_AA = all_out_A.all_corr_an{1};

if group == 1
    sel_out = all_CC_C;
    sel_p = all_out_C;
else
    sel_out = all_CC_A;
    sel_p = all_out_A;
end

ff = makeFigureRowsCols(107,[1 0.5 2 2],'RowsCols',[4 4],...
    'spaceRowsCols',[0.05 0.07],'rightUpShifts',[0.1 0.1],'widthHeightAdjustment',...
    [-95 -95]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[1 6 3.4 3.4]);
FS = mData.axes_font_size;
maskdisp = triu(ones(4,4),0);

for rr = 1:4
    for cc = 1:4
        if ~maskdisp(rr,cc)
            delete(ff.h_axes(rr,cc));
            continue;
        end
        axes(ff.h_axes(rr,cc));
        imagesc(sel_out{rr,cc});
        box off;
        set(gca,'Ydir','Normal','linewidth',1,'FontSize',FS-1,'FontWeight','Bold');
        if rr == cc
            cols = size(sel_out{rr,cc},2);
            colsHalf = round(cols/2);
            ts = round(sel_p.xs{3,1}(1:cols));
            set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1)-2 ts(colsHalf) ts(cols)+2]);
            h = xlabel('Position (cm)');%    changePosition(h,[0 0 0]);
            h = ylabel('Position (cm)');%    changePosition(h,[0 0 0]);
            cols = size(sel_out{rr,cc},2);
            colsHalf = round(cols/2);
            ts = round(sel_p.xs{3,1}(1:cols));
            set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1)-2 ts(colsHalf) ts(cols)+2]);

        else
            axis off
        end
        minC = min(sel_out{rr,cc}(:));
        maxC = max(sel_out{rr,cc}(:));
        hc = putColorBar(ff.h_axes(rr,cc),[0.0 0.03 0 -0.05],[minC maxC],5,'eastoutside',[0.07 0.08 0.1 0.13]);
    end
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('pos corr remap_%d',group),600);
return;
