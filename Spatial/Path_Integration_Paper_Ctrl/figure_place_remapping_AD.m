function figure_place_remapping_AD(fn,allRs,ccs)

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = 1.96; fwids = [1 150]; fcens = [0 150]; rs_th = 0.3;
conditionsAndRasterTypes = [11;21;31;41];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN);
out = read_data_from_base_workspace(selC)

ei_C = out.eis{1}; ei_A = out.eis{2};
pMs_C = out.pMs{1}; pMs_A = out.pMs{2};
paramMs_C = out.paramMs{1}; paramMs_A = out.paramMs{2};
cpMs_C = out.cpMs{1}; cpMs_A = out.cpMs{2};
selAnimals_C = out.selAnimals{1}; selAnimals_A = out.selAnimals{2};
perc_cells_C = out.perc_cells{1}; perc_cells_A = out.perc_cells{2};

if 0
    cellSel_C = ''; cellSel_A = ''; cellSel_C = cpMs_C; cellSel_A = cpMs_A;
    a_trials = {3:10};
    for ii = 1:length(a_trials)
        out = get_mean_rasters(pMs_C',paramMs_C,selAnimals_C,ei_C,conditionsAndRasterTypes',selC,cellSel_C,a_trials{ii});
        [out.allP_an,out.allC_an,out.avg_C_conds,out.mean_rasters_T,out.all_corr_an] = get_pop_vector_corr(out,conditionsAndRasterTypes,min(out.sz(:)));
        all_out_C{ii} = out;
        out = get_mean_rasters(pMs_A',paramMs_A,selAnimals_A,ei_A,conditionsAndRasterTypes',selC,cellSel_A,a_trials{ii});
        [out.allP_an,out.allC_an,out.avg_C_conds,out.mean_rasters_T,out.all_corr_an] = get_pop_vector_corr(out,conditionsAndRasterTypes,49);
        all_out_A{ii} = out;
    end
    save('pop_corr_mat_remapping_or.mat','all_out_C','all_out_A','a_trials');
    disp('Done');
    return;
else
    temp = load('pop_corr_mat_remapping_and.mat');
    all_out_C_and = temp.all_out_C{1};     all_out_A_and = temp.all_out_A{1}; a_trials = temp.a_trials{1};
    temp = load('pop_corr_mat_remapping_or.mat');
    all_out_C_or = temp.all_out_C{1};     all_out_A_or = temp.all_out_A{1}; a_trials = temp.a_trials{1};
end

group = 2; andor = 1;

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
        for ii = 1:length(all_out_C.all_corr_an)
            thisCC(:,:,ii) = all_out_C.all_corr_an{ii}{rr,cc};
        end
        all_CC_C{rr,cc} = nanmean(thisCC,3);
        
        thisCC = [];
        for ii = 1:length(all_out_A.all_corr_an)
            thisCC(:,:,ii) = all_out_A.all_corr_an{ii}{rr,cc};
        end
        all_CC_A{rr,cc} = nanmean(thisCC,3);
    end
end


n = 0;
%%
% out_CC = all_out_C.all_corr_an{3};
% out_AA = all_out_A.all_corr_an{1};

if group == 1
    sel_out = all_CC_C;
else
    sel_out = all_CC_A;
end

ff = makeFigureRowsCols(107,[1 0.5 2 2],'RowsCols',[4 4],...
    'spaceRowsCols',[0.01 0.01],'rightUpShifts',[0.03 0.1],'widthHeightAdjustment',...
    [-20 -60]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[1 6 4 4]);
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
%         if rr ~= 4 & cc ~=1
            axis off
%         end
    end
end
return;
%%
for sii = 1:length(conditionsAndRasterTypes)

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


