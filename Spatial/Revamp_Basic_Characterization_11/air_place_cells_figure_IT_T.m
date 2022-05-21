function air_place_cells_figure_IT_T

%%
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei = evalin('base','ei'); 
    selContexts = [3 4 5]; rasterNames = {'airIT','airIT','airIT'};
    oT = get_data(ei,selContexts,rasterNames);
    selContexts = [3 4 5]; rasterNames = {'airID','airID','airID'};
    oD = get_data(ei,selContexts,rasterNames);
    selContexts = [3 4 5 3 4 5]; rasterNames = {'airIT','airIT','airIT','airID','airID','airID'};
    oDT = get_data(ei,selContexts,rasterNames);
    break;
end
n = 0;

%% Show sample rasters for presentation
% an = 5; cn = 2;
% plotRasters_simplest(Rs{an,cn})
si = [Ar_i_T ArL_i_T Ars_i_T];
Rs = o.Rs(:,si);mR = o.mR(:,si);
props1 = get_props_Rs(Rs,50);
while 1
   an = 4; cn = 2;
    % plotRasters_simplest(Rs{an,cn})
    % find(resp_valsC{an}(:,cn));
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
        [-75 -475]);
    set(gcf,'color','w'); set(gcf,'Position',[10 4 3.25 1]);
%     ff = sample_rasters(Rs{an,cn},[75 278 145 181],ff);
%     ff = sample_rasters(Rs{an,cn},[75 71 115 6],ff); % tuned
    ff = sample_rasters(Rs{an,cn},[3 33 7 37],ff); % untuned
    for ii = 1:length(ff.h_axes(1,:)) set_obj(ff.h_axes(1,ii),{'xtick',[1 68 136],'xticklabels',[0 7.5 15]});  end 
    colormap parula
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersIT'),600);
    
%     ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
%         'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
%         [-75 -475]);
%     set(gcf,'color','w'); set(gcf,'Position',[10 4 3.25 1]);
%     ff = sample_rasters(RsT{an,cn},[191 11 96 41],ff);
%     save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersDT'),600);
    break;
end


%% find correlations across conditions


Rs = oT.Rs; mR = oT.mR;
props1 = get_props_Rs(Rs,50);
gFR = props1.good_FR;
gFR_OR = cell_list_op(gFR,[],'or');

% resp = resp_ORCS;
% [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,1,0);

outRemap = find_population_vector_corr_remap(Rs,mR,gFR_OR);
n = 0;

%% look at PV correlation across conditions
% average correlation of all animals
while 1
    ff = makeFigureRowsCols(106,[1 0.5 9 9],'RowsCols',[3 3],...
        'spaceRowsCols',[0.1 0.1],'rightUpShifts',[0.1 0.1],'widthHeightAdjustment',...
        [-130 -130]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 3 3 3]);
    ff = show_remapping_corr_plots_time(mData,ff,outRemap.mean_PV_corr,outRemap.xs,[]);
    colormap_ig
    save_pdf(ff.hf,mData.pdf_folder,sprintf('remap_corr_C_IT.pdf'),600);
    break;
end

%% look at percentage of spatially tuned cells, place field widths on the belt based on location of their peak in one of the three bins
% In this code make changes manually for plotting spatially tuned cells or
% place field widths
while 1
    try
        Rs = oD.Rs;
    catch
        si = si_seq(setdiff(1:11,[1 11 9 2 10]));
        si = si([1 3 5]);
        Rs = o.Rs(:,si);mR = o.mR(:,si);
    end
    props1 = get_props_Rs(Rs,50);
    gFR = props1.good_FR;
    pL = props1.peak_locations;
    rs = props1.centers;
    propI = props1.PWs;
    minBin = 0;     maxBin = 150;     incr = 50; % choosing more than 3 bins, give significant anova but not significant multcompare
    % three bins also make sense because LED in condition 2 comes ON at
    % 110cm
    bins = minBin:incr:maxBin;
    vals_rs = []; vals = []; vals_propI = [];
    for rr = 1:size(pL,1)
        one_row = []; one_row_rs = []; one_row_propI = [];
        for cc = 1:size(pL,2)
            tpL = pL{rr,cc};
            tgFR = gFR{rr,cc};
            tpL = tpL(tgFR);
            trs = rs{rr,cc}; tpropI = propI{rr,cc};
            tBinV = []; tBinV_rs = []; tBinV_propI = [];
            for bb = 2:length(bins)
                tBinV(bb-1) = 100*length(find(tpL > bins(bb-1) & tpL <= bins(bb)))/length(tgFR);
                inds_cells = find(tpL > bins(bb-1) & tpL <= bins(bb));
                tBinV_rs(bb-1) = nanmean(trs(inds_cells)-tpL(inds_cells));
                tBinV_propI(bb-1) = nanmean(tpropI(inds_cells));
            end
            one_row = [one_row tBinV]; one_row_rs = [one_row_rs tBinV_rs]; one_row_propI = [one_row_propI tBinV_propI];
        end
        vals(rr,:) = one_row; vals_rs(rr,:) = one_row_rs; vals_propI(rr,:) = one_row_propI;
    end
    [within,dvn,xlabels] = make_within_table({'Cond','Bins'},[3,3]);
    % see here which variable is being set up for the anova
    dataT = make_between_table({vals_propI},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_Bins','hsd'},[1 1 1]);
    h(h==1) = 0;
%     s = generate_shades(3); tcolors = s.g;
    tcolors = [mData.colors(6);mData.colors(6);mData.colors(6);...
                mData.colors(7);mData.colors(7);mData.colors(7);...
                mData.colors(8);mData.colors(8);mData.colors(8);];
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.75 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
    format_axes(gca);
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[0 maxY])
    xticks = [xdata(1:end)]; xticklabels = {'B1','B2','B3'};
%     set_title(gca,'Air-Belt Mismatch',[-0.35,95],5)
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.08 0.02 -0.005 -0.011])
    
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
%     put_axes_labels(gca,{[],[0 0 0]},{{'Place Field','Width (cm)'},[0 0 0]});
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Spatially Tuned Cells_IT'),600);
    
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);%h(h==1) = 0;
%     s = generate_shades(3); tcolors = s.g;
    tcolors = mData.dcolors(8:10);
     hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 3 1.5 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
    format_axes(gca);
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[0 maxY])
    xticks = [xdata(1:end)]; xticklabels = {'B1','B2','B3'};
%     set_title(gca,'Air-Belt Mismatch',[-0.35,95],5)
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.12 0.02 -0.35 -0.091])
    
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Spatially Tuned Cells_pooled_IT'),600);
    
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
    save_pdf(hf,mData.pdf_folder,'air_place_cell corr remap _IT',600);
    break;
end

%% Delta FR score (RAte remap) between conditions
while 1
    props1 = get_props_Rs(Rs,50);
    gFR = props1.good_FR;
    for an = 1:length(outRemap.all_RR_SP)
        this_sp_corr_diag = outRemap.all_RR_SP{an};
        var_CE_mat{an,1} = this_sp_corr_diag{1,2};
        var_CE_mat{an,2} = this_sp_corr_diag{2,3};
        var_CE_mat{an,3} = this_sp_corr_diag{1,3};
    end
    var_CE1 = exec_fun_on_cell_mat(var_CE_mat,'nanmean'); 
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
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0.7 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xlabels = {'Ar-ArL','ArL-Ar*','Ar-Ar*'};
    xticks = xdata(1:end)+0; xticklabels = xlabels;
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     set(hbs(2),'facecolor','none','edgecolor',tcolors{2});
    xtickangle(45);
%     changePosition(gca,[0.2 0.03 -0.55 -0.05]);
    changePosition(gca,[0.1 0.03 -0.4 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'\Delta FR Score'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'air_place_cell RR corr remap_IT',600);
    break;
end

%% Popuation Vector correlation between conditions
while 1
    props1 = get_props_Rs(Rs,50);
    gFR = props1.good_FR;
    for an = 1:length(outRemap.all_PV_corr_diag)
        this_sp_corr_diag = outRemap.all_PV_corr_diag{an};
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
    put_axes_labels(gca,{[],[0 0 0]},{{'PV Correlation'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'air_place_cell PV corr remapIT',600);
    break;
end


%% Place cell emergence disruption stability
props1 = get_props_Rs(Rs,50);
respAnB = props1.good_FR;

for tC = 1:4
while 1
    txtT = {'Unique','New','Disrupted','Common'};
    or_cells = get_cell_list(respAnB,[1;2;3]);
    percCells = []; 
    switch tC
        case 1 % unique place cells
            respCells1 = get_cell_list(respAnB,[1 -2 -3]);    respCells2 = get_cell_list(respAnB,[-1 2 -3]);    respCells3 = get_cell_list(respAnB,[-1 -2 3]);
        case 2 % new place cells
            respCells1 = get_cell_list(respAnB,[-1 2]);    respCells2 = get_cell_list(respAnB,[-2 3]);    respCells3 = get_cell_list(respAnB,[-1 3]);
        case 3 % disrupted place cells
            respCells1 = get_cell_list(respAnB,[1 -2]);    respCells2 = get_cell_list(respAnB,[2 -3]);    respCells3 = get_cell_list(respAnB,[1 -3]);
        case 4 % remained place cells
            respCells1 = get_cell_list(respAnB,[1 2]);    respCells2 = get_cell_list(respAnB,[2 3]);    respCells3 = get_cell_list(respAnB,[1 3]);
    end
    
    respCells = [respCells1(:,1) respCells2(:,1) respCells3(:,1)];
    percCells = 100*exec_fun_on_cell_mat(respCells,'sum')./exec_fun_on_cell_mat(or_cells,'length');

    [within,dvn,xlabels] = make_within_table({'Conds'},[3]);
    if tC > 1
       xlabels = {'Ar-ArL','ArL-Ar*','Ar-Ar*'};
    else
        xlabels = {'Ar','ArL','Ar*'};
    end
    dataT = make_between_table({percCells},dvn);
    ra = RMA(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Conds','hsd'},[1 0.25 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    if tC > 1
        tcolors = mData.dcolors(6:8);%[s.m;s.c;s.y];
    else
        tcolors = mData.colors(6:8);%[s.m;s.c;s.y];
    end
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);
    format_axes(gca);
    xticks = xdata; xticklabels = xlabels;
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.1 0.02 -0.45 -0.05])
    yshifts = [0 0 -2 0]
    put_axes_labels(gca,{[],[0 0 0]},{sprintf('%s Cells (%%)',txtT{tC}),[0 yshifts(tC) 0]});
%     text(0.75,maxY+3,txtT{tC},'FontSize',6)
    save_pdf(hf,mData.pdf_folder,sprintf('%s_cells_IT.pdf',txtT{tC}),600);
    ra.ranova
%     ra.mauchly
    break;
end
end
