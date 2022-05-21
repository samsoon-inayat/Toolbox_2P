function air_place_cells_figure
%%
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei = evalin('base','ei'); 
    selContexts = [3 4 5]; rasterNames = {'airT','airT','airT'};
    oT = get_data(ei,selContexts,rasterNames);
    selContexts = [3 4 5]; rasterNames = {'airD','airD','airD'};
    oD = get_data(ei,selContexts,rasterNames);
    selContexts = [3 4 5 3 4 5]; rasterNames = {'airD','airD','airD','airT','airT','airT'};
    oDT = get_data(ei,selContexts,rasterNames);
    break;
end
n = 0;
%% find correlations across conditions

try
    Rs = oD.Rs; mR = oD.mR;
catch
    si = si_seq(setdiff(1:11,[1 11 9 2 10]));
    si = si([1 3 5]);
    Rs = o.Rs(:,si);mR = o.mR(:,si);
end
props1 = get_props_Rs(Rs,50);
gFR = props1.good_FR;
gFR_OR = cell_list_op(gFR,[],'or');

% resp = resp_ORCS;
% [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,1,0);

outRemap = find_population_vector_corr_remap(Rs,mR,gFR_OR);
n = 0;

%% find trial by trial comparison - correlations
trials = mat2cell([1:10]',ones(size([1:10]')));
% trials = mat2cell([1:3]',ones(size([1:3]')));
resp = gFR_OR;
out1 = find_population_vector_corr_remap_trials(Rs(:,1),resp,trials);
out2 = find_population_vector_corr_remap_trials(Rs(:,2),resp,trials);
out3 = find_population_vector_corr_remap_trials(Rs(:,3),resp,trials);
% save(fullfile(mData.pd_folder,'air_place_cells_figure.mat'),'out1','out2','out3');
n = 0;
%% look at PV correlation across trials
% average correlation of all animals
while 1
    selORemap = out2;
    ff = makeFigureRowsCols(106,[1 0.5 9 9],'RowsCols',[10 10],...
        'spaceRowsCols',[0.013 0.013],'rightUpShifts',[0.06 0.051],'widthHeightAdjustment',...
        [-20 -20]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 1 4 4]);
    ff = show_remapping_corr_plots(mData,ff,selORemap.mean_PV_corr,selORemap.xs,[]);
    for rr = 1:size(ff.h_axes,1)
        for cc = 1:size(ff.h_axes,2)
            if rr == cc
                if rr > 1 & rr < 10
                    ff.h_axes(rr,cc).YLabel.String = '';                    ff.h_axes(rr,cc).XTick = [];                    ff.h_axes(rr,cc).YTick = []; ff.h_axes(rr,cc).XLabel.String = '';
                end
                if rr == 1
                    ff.h_axes(rr,cc).XLabel.String = '';                    ff.h_axes(rr,cc).XTick = [];
                end
                if rr == 10
                    ff.h_axes(rr,cc).YLabel.String = '';                    ff.h_axes(rr,cc).YTick = [];
                end
            end
        end
    end
    
    colormap_ig
    save_pdf(ff.hf,mData.pdf_folder,sprintf('remap_corr_C_trials.pdf'),600);
    break;
end

%% Spatial correlation across conditions
while 1
    props1 = get_props_Rs(Rs,50);
    gFR = props1.good_FR;
    for an = 1:length(outRemap.all_SP_corr_diag)
        this_sp_corr_diag = outRemap.all_SP_corr_diag{an};
        var_CE_mat{an,1} = this_sp_corr_diag{1,2};
        var_CE_mat{an,2} = this_sp_corr_diag{2,3};
        var_CE_mat{an,3} = this_sp_corr_diag{1,3};
    end
    var_CE1 = exec_fun_on_cell_mat(var_CE_mat,'mean'); 
%     var_CE = exec_fun_on_cell_mat(outRemap.adj_SP_corr_diag,'mean');
    [within,dvn,xlabels] = make_within_table({'Cond'},3);
    dataT = make_between_table({var_CE1},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 0 1]);
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors = mData.dcolors(6:8);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.02,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xlabels = {'Ar-ArL','ArL-Ar*','Ar-Ar*'};
    xticks = xdata(1:end)+0; xticklabels = xlabels;
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     set(hbs(2),'facecolor','none','edgecolor',tcolors{2});
    xtickangle(45);
%     changePosition(gca,[0.2 0.03 -0.55 -0.05]);
    changePosition(gca,[0.1 0.03 -0.4 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatial Correlation'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'air_place_cell corr remap',600);
    break;
end
