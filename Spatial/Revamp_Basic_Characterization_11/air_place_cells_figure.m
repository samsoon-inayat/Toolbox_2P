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
    si = [Ar_t_D ArL_t_D Ars_t_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
end
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
    ff = show_remapping_corr_plots(mData,ff,outRemap.mean_PV_corr,outRemap.xs,[]);
    colormap_ig
    save_pdf(ff.hf,mData.pdf_folder,sprintf('remap_corr_C.pdf'),600);
    break;
end

%% find trial by trial comparison - correlations
trials = mat2cell([1:10]',ones(size([1:10]')));
resp = get_cell_list(resp_valsCS,[]);
out1 = find_population_vector_corr_remap_trials(Rs(:,1),resp,trials);
out2 = find_population_vector_corr_remap_trials(Rs(:,2),resp,trials);
out3 = find_population_vector_corr_remap_trials(Rs(:,3),resp,trials);
% save(fullfile(mData.pd_folder,'air_place_cells_figure.mat'),'out1','out2','out3');
n = 0;
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
    dataT = make_between_table({vals},dvn);
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
    save_pdf(hf,mData.pdf_folder,sprintf('Spatially Tuned Cells'),600);
    
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
    save_pdf(hf,mData.pdf_folder,sprintf('Spatially Tuned Cells_pooled'),600);
    
    break;
end
%% Air Belt mismatch
while 1
    sRs = Rs(:,1:3);
%     sRs = Rs(:,4:6);
    for rr = 1:size(sRs,1)
        for cc = 1:size(sRs,2)
            R = sRs{rr,cc};
            sp = R.space;
            [msp,mspi] = min(sp,[],2);
            air_belt_distance(rr,cc) = mean((mspi-1) * R.bin_width);
%             figure(1000);clf;imagesc(R.space);
        end
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},[3]);
    dataT = make_between_table({air_belt_distance},dvn);
    ra = RMA(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','bonferroni'},[1 1 1]);
%     s = generate_shades(3); tcolors = s.g;
    tcolors = mData.colors(1:3);
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
    format_axes(gca);
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[0 maxY])
    xticks = [xdata(1:end)]; xticklabels = {'C3','C4','C3'''};
    set_title(gca,'Air-Belt Mismatch',[-0.35,95],5)
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.17 0.02 -0.45 -0.011])
    
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Avg. Air-Onset to Belt', 'Marker Dist (cm)'},[0 -5 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Air_Belt_Mismatch'),600);
    break;
end
%% Show sample rasters
% an = 5; cn = 2;
% plotRasters_simplest(Rs{an,cn})
si = [Ar_t_D ArL_t_D Ars_t_D];
Rs = o.Rs(:,si);mR = o.mR(:,si);
while 1
   an = 3; cn = 1;
    % plotRasters_simplest(Rs{an,cn})
    % find(resp_valsC{an}(:,cn));
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
        [-75 -475]);
    set(gcf,'color','w'); set(gcf,'Position',[10 4 3.25 1]);
    ff = sample_rasters(Rs{an,cn},[191 11 96 41],ff);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersD'),600);
    
%     ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
%         'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
%         [-75 -475]);
%     set(gcf,'color','w'); set(gcf,'Position',[10 4 3.25 1]);
%     ff = sample_rasters(RsT{an,cn},[191 11 96 41],ff);
%     save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersDT'),600);
    break;
end

%% Show sample rasters for presentation
% an = 5; cn = 2;
% plotRasters_simplest(Rs{an,cn})
si = [Ar_t_D ArL_t_D Ars_t_D];
Rs = o.Rs(:,si);mR = o.mR(:,si);
while 1
   an = 3; cn = 1;
    % plotRasters_simplest(Rs{an,cn})
    % find(resp_valsC{an}(:,cn));
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[2 4],...
        'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
        [-75 -475]);
    set(gcf,'color','w'); set(gcf,'Position',[10 4 3.25 1]);
    ff = sample_rasters(Rs{an,cn},[191 11 96 41],ff);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersD'),600);
    colormap parula
%     ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
%         'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
%         [-75 -475]);
%     set(gcf,'color','w'); set(gcf,'Position',[10 4 3.25 1]);
%     ff = sample_rasters(RsT{an,cn},[191 11 96 41],ff);
%     save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersDT'),600);
    break;
end


%% Mutual Information Time versus Distance Distributions
while 1
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            R = Rs{rr,cc};
            zMIsC{rr,cc} = R.info_metrics.ShannonMI_Zsh';
            R = RsT{rr,cc};
            zMIsA{rr,cc} = R.info_metrics.ShannonMI_Zsh';
        end
    end

    CN = 1;
    tcolors = {'m','c'};
    distD(:,1) = zMIsC(:,CN);
    distD(:,2) = zMIsA(:,CN);
    [distDo,allVals,allValsG] = plotDistributions(distD);
    minBin = min(allVals);
    maxBin = max(allVals);

    
    incr = 0.001; %maxBin =
    hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    %    [ha,hb,hca] = plotDistributions(distD,'colors',tcolors,'maxY',maxBin,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    [ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
    plot([1.65 1.65],[0 100],'--k');
    myxlims2 = [7 9 13];
    xlims = xlim; xlim([xlims(1) myxlims2(CN)]);
    set(gca,'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
    changePosition(gca,[0.15 0.13 -0.4 -0.13]);
    if CN == 1
        put_axes_labels(gca,{'zMI',[0 0 0]},{{'Percentage','of Neurons'},[0 0 0]});
    else
        put_axes_labels(gca,{'zMI',[0 0 0]},{'',[1 0 0]});
    end
    if CN == 3
        legs = {'Dist (D)','Time (T)'};
        legs{end+1} = [4 2 60 5];
        putLegend(gca,legs,'colors',tcolors)
    end
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_zMI_%d',CN),600);
    break;
end
%% Mutual Information Time versus Distance bar graph
mzMIsC = arrayfun(@(x) nanmean(x{1}),zMIsC);
mzMIsA = arrayfun(@(x) nanmean(x{1}),zMIsA);
[within,dvn,xlabels] = make_within_table({'DT','Cond'},[2,3]);
dataT = make_between_table({mzMIsC,mzMIsA},dvn);
ra = repeatedMeasuresAnova(dataT,within,0.05);
[xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0.25 1]);
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
hold on;
tcolors = [mData.shades.m;mData.shades.c];
% tcolors = repmat(tcolors,2,1)';
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
set(gca,'xlim',[0.35 xdata(end)+.65],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
xticks = xdata; xticklabels = {'C3-D','C4-D','C3''-D','C3-T','C4-T','C3''-T'}; xticklabels = repmat(xticklabels,1,2);
set(gca,'xtick',xticks,'xticklabels',xticklabels);
xtickangle(45)
changePosition(gca,[0.05 0.02 -0.05 -0.011])
put_axes_labels(gca,{[],[0 0 0]},{{'zMI'},[0 0 0]});

save_pdf(hf,mData.pdf_folder,'zMI_bar_graph_all_cells.pdf',600);



%% spatial and population vector correlation single animal
while 1
    an = 1; cn = 2;
    ff = makeFigureRowsCols(106,[1 0.5 6 1],'RowsCols',[2 2],...
        'spaceRowsCols',[0.09 -0.03],'rightUpShifts',[0.07 0.1],'widthHeightAdjustment',...
        [-100 -150]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 5 3.25 2]);
%     [resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC,resp_exc_inh] = get_responsive_fraction(Rs);
%     resp = get_cell_list(resp_valsC,[]);
    [CRc,aCRc,mRR,spCR] = find_population_correlations(Rs,mR,1,0);
    axes(ff.h_axes(1,1)); imagesc(mRR{an,cn}); colorbar; set(gca,'Ydir','Normal');
    axes(ff.h_axes(2,1)); imagesc(CRc{an,cn}); colorbar; set(gca,'Ydir','Normal');
    axes(ff.h_axes(1,2)); imagesc(spCR{an,cn}); colorbar; set(gca,'Ydir','Normal');
    set(ff.h_axes(2,2),'Visible','Off');
    set_colormap;
    set_obj(ff,{'FontWeight','Normal','FontSize',6,'LineWidth',0.25});
    ht = get_obj(ff,'title'); hyl = get_obj(ff,'ylabel'); changePosition(hyl(1,1),[4 0 0]);
%     set_obj(ht,'String',{'Pop. Activity','Pop. Activity','Pop. Activity';'Pop. Correlation','Pop.Correlation','Pop.Correlation'});
    set_obj(ht,{'FontSize',5,'FontWeight','Normal'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_correlations.pdf'),600);
    break;
end

%% Percentage of Responsive Cells Unique in Conditions
while 1
    percCells = [];
    percCells(:,1) = get_cell_list(resp_valsCS,[1 -2 -3],1)'; percCells(:,2) = get_cell_list(resp_valsCS,[-1 2 -3],1)'; percCells(:,3) = get_cell_list(resp_valsCS,[-1 -2 3],1)';
    within = make_within_table({'Cond'},3);
    dataT = make_between_table({percCells*100},{'C31','C4','C32'});
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
%     s = generate_shades(3); tcolors = s.g;
    tcolors = mData.colors;
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
    format_axes(gca);
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[0 25])
    xticks = [xdata(1:end)]; xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    any_mean = mean(100*resp_OR_fractionCS);    any_sem = std(100*resp_OR_fractionCS)/sqrt(5);
%     pmchar=char(177); any_text = sprintf('%.0f%c%.0f%%',any_mean,pmchar,any_sem); text(0.75,20,any_text,'FontSize',6);
    changePosition(gca,[0.17 0.02 -0.4 -0.011])
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('percentage_PCs_responsive_unique'),600);
    break;
end
%% Spatial correlation between adjacent trails
while 1
    [within,dvn,xlabels] = make_within_table({'Condition','TrialPairs'},[3,9]);
    var1 = arrayfun(@(x) mean(x{1}),out1.adj_SP_corr_diag); 
    var2 = arrayfun(@(x) mean(x{1}),out2.adj_SP_corr_diag);
    var3 = arrayfun(@(x) mean(x{1}),out3.adj_SP_corr_diag);
%     var1 = arrayfun(@(x) mean(x{1}),out1.adj_PV_corr_diag); 
%     var2 = arrayfun(@(x) mean(x{1}),out2.adj_PV_corr_diag);
%     var3 = arrayfun(@(x) mean(x{1}),out3.adj_PV_corr_diag);
    between = make_between_table({var1,var2,var3},dvn);
    ra = RMA(between,within);
    ra.ranova
    %%
    hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_line_graph(mData,ra,0);
    tcolors = mData.colors;
    ii = 1; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 2; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linestyle','-.','linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 3; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    legs = {'C3','C4','C3'''};
    legs{end+1} = [1 0.2 0.25 0.25];
    putLegendH(gca,legs,tcolors)
    ylims = ylim;
    set(gca,'xlim',[0.5 xdata(end)+0.5],'ylim',[ylims(1) ylims(2)],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'T1-T2','T2-T3','T3-T4','T4-T5','T5-T6','T6-T7','T7-T8','T8-T9','T9-T10'};
    xticklabels = repmat(xticklabels,1,2);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(45);
    changePosition(gca,[0.05 0.12 0.05 -0.08]);
    put_axes_labels(gca,{'Trial-Pairs',[0 0 0]},{{'Spatial Correlation'},[0 -0.05 0]});
    save_pdf(hf,mData.pdf_folder,'spatial_correlation_trials',600);
    %%
    break;
end
%% Spatial correlation between adjacent trails (taking mean) 
while 1
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Condition_by_TrialPairs','hsd'},[1 0.5 1]);
    hollowsep = 19;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.75 1],'color','w');
    hold on;
    tcolors = colors(1:9); tcolors = repmat(tcolors,[1 3]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
            'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
            'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    ylims = ylim;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[ylims(1) maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'T1-T2','T2-T3','T3-T4','T4-T5','T5-T6','T6-T7','T7-T8','T8-T9','T9-T10'};
    xticklabels = repmat(xticklabels,1,3);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(45);
%     changePosition(gca,[0.1 0.11 -0.2 -0.05]);
    changePosition(gca,[0.0 0.03 0.08 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatial','Correlation'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'spatial_correlation_trials_bar_graph',600);
    break;
end

%% Spatial correlation trails with mean of trials
while 1
    [within,dvn,xlabels] = make_within_table({'Condition','TrialPairs'},[3,10]);
    var1 = []; var2 = []; var3 = [];
    for ii = 1:length(out1.all_SP_corr_diag_with_mean)
        var1(ii,:) = (arrayfun(@(x) mean(x{1}),out1.all_SP_corr_diag_with_mean{ii}))'; 
        var2(ii,:) = (arrayfun(@(x) mean(x{1}),out2.all_SP_corr_diag_with_mean{ii}))'; 
        var3(ii,:) = (arrayfun(@(x) mean(x{1}),out3.all_SP_corr_diag_with_mean{ii}))'; 
    end
    between = make_between_table({var1,var2,var3},dvn);
    ra = repeatedMeasuresAnova(between,within);
    %%
    hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_line_graph(mData,ra,0);
    tcolors = mData.colors;
    ii = 1; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 2; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linestyle','-.','linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 3; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    legs = {'C3','C4','C3'''};
    legs{end+1} = [1 0.2 0.25 0.25];
    putLegendH(gca,legs,tcolors)
    ylims = ylim;
    set(gca,'xlim',[0.5 xdata(end)+0.5],'ylim',[ylims(1) ylims(2)],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'T1-T2','T2-T3','T3-T4','T4-T5','T5-T6','T6-T7','T7-T8','T8-T9','T9-T10'};
    xticklabels = repmat(xticklabels,1,2);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(45);
    changePosition(gca,[0.05 0.12 0.05 -0.08]);
    put_axes_labels(gca,{'Trial-Pairs',[0 0 0]},{{'Spatial Correlation'},[0 -0.05 0]});
    save_pdf(hf,mData.pdf_folder,'spatial_correlation_trials',600);
end
%% Spatial correlation trails with mean of trials (mean bar graph)
while 1
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0.5 1]);
    hollowsep = 19;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 3.25 1],'color','w');
    hold on;
    tcolors = colors(1:10); tcolors = repmat(tcolors,[1 3]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
            'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
            'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    ylims = ylim;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[ylims(1) maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'T1','T2','T3','T4-T5','T5-T6','T6-T7','T7-T8','T8-T9','T9-T10'};
    xticklabels = repmat(xticklabels,1,3);
    set(gca,'xtick',xticks,'xticklabels',xlabels);
    xtickangle(45);
%     changePosition(gca,[0.1 0.11 -0.2 -0.05]);
    changePosition(gca,[-0.02 0.03 0.1 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Correlation'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'spatial_correlation_trials_bar_graph',600);
end
%% Rate remap (Delta FR score) between adjacent trails
while 1
    [within,dvn,xlabels] = make_within_table({'Condition','TrialPairs'},[3,9]);
    var1 = arrayfun(@(x) nanmean(x{1}),out1.adj_RR_SP); 
    var2 = arrayfun(@(x) nanmean(x{1}),out2.adj_RR_SP);
    var3 = arrayfun(@(x) nanmean(x{1}),out3.adj_RR_SP);
    between = make_between_table({var1,var2,var3},dvn);
    ra = repeatedMeasuresAnova(between,within);
    %%
    hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_line_graph(mData,ra,0);
    tcolors = mData.colors;
    ii = 1; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 2; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linestyle','-.','linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 3; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    legs = {'C3','C4','C3'''};
    legs{end+1} = [1.5 0.3 1.02 0.25];
    putLegendH(gca,legs,tcolors)
    ylims = ylim;
    set(gca,'xlim',[0.5 xdata(end)+0.5],'ylim',[ylims(1)-0.01 ylims(2)],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'T1-T2','T2-T3','T3-T4','T4-T5','T5-T6','T6-T7','T7-T8','T8-T9','T9-T10'};
    xticklabels = repmat(xticklabels,1,2);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(45);
    changePosition(gca,[0.065 0.12 0.05 -0.08]);
    put_axes_labels(gca,{'Trial-Pairs',[0 0 0]},{{'\Delta FR Score'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'Delta_FR_score_trials',600);
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
    save_pdf(hf,mData.pdf_folder,'air_place_cell RR corr remap',600);
    break;
end

% while 1
%     [within,dvn,xlabels] = make_within_table({'Cond'},2);
%     var_CE = arrayfun(@(x) nanmean(x{1}),outRemap.adj_RR_SP);
%     dataT = make_between_table({var_CE},dvn);
%     ra = RMA(dataT,within);
%     ra.ranova
%     [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 0 1]);
%     colors = mData.colors;
%     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
%     hold on;
%     tcolors ={colors{7};colors{8};};
%     [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
%         'ySpacing',0.02,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
%         'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+0.1],'FontSize',6,'FontWeight','Normal','TickDir','out');
%     xticks = xdata(1:end)+0; xticklabels = {'C3-C4','C4-C3'''};
%     set(gca,'xtick',xticks,'xticklabels',xticklabels);
% %     set(hbs(2),'facecolor','none','edgecolor',tcolors{2});
%     xtickangle(45);
% %     changePosition(gca,[0.2 0.03 -0.55 -0.05]);
%     changePosition(gca,[0.1 0.03 -0.4 -0.1]);
%     put_axes_labels(gca,{[],[0 0 0]},{{'\Delta FR Score'},[0 0 0]});
%     save_pdf(hf,mData.pdf_folder,'air_place_cell rate remap (DFR)',600);
% end

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
    save_pdf(hf,mData.pdf_folder,'air_place_cell PV corr remap',600);
    break;
end

%%
%% overlap RM bar graph
while 1
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    resp_valsCS = gFR;
    OI_E = get_overlap_index(resp_valsCS);
    for ii = 1:5
        C12(ii,1) = OI_E{ii}(1,2);
        C12(ii,2) = OI_E{ii}(2,3);
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},2);
    var_CE = C12;
    dataT = make_between_table({var_CE},dvn);
    ra = repeatedMeasuresAnova(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0 1]);
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors ={colors{7};colors{8};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+0.1],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C3-C4','C4-C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     set(hbs(2),'facecolor','none','edgecolor',tcolors{2});
    xtickangle(45);
%     changePosition(gca,[0.2 0.03 -0.55 -0.05]);
    changePosition(gca,[0.1 0.03 -0.4 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Overlap Index'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('air_overlap_stats_place_cells'),600);
end

%% Place Field Properties
% Rs = oD.Rs;
props1 = get_props_Rs(Rs,50);
respAnB = cell_list_op(props1.good_FR,props1.good_Gauss,'and');
% good_FR_OR = cell_list_op(good_FR,[],'or');
props = {'Field Width (cm)','Field Center (cm)','Field FR (AU)'};
fileNames = {'PFW','PFC','PFFR'};
all_incr = [1 1 1];
% legs_v = 
for pri = 3%:3;
while 1
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            R = Rs{rr,cc};
            [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
            switch pri
                case 1
                    zMIsC{rr,cc} = PWs(respAnB{rr,cc})';
                case 2
                    zMIsC{rr,cc} = centers(respAnB{rr,cc})';
                case 3
                    zMIsC{rr,cc} = MFR(respAnB{rr,cc})';
            end
        end
    end
%%
%     tcolors = mData.colors;
%     distD = zMIsC;
%     [distDo,allVals] = getAveragesAndAllValues(distD);
%     minBin = min(allVals);
%     maxBin = max(allVals);
%     incr = all_incr(pri); %maxBin =
%     hf = get_figure(8,[5 7 1.25 1]);
%     [ha,hb,~,bins] = plotAverageDistributions(distD,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
%     format_axes(gca);
%     if pri == 1
%         legs = {'C3','C4','C3''',[50 3 70 5]}; putLegend(gca,legs,tcolors);
%     end
%     changePosition(gca,[0.1 0.13 -0.25 -0.13]);
%     put_axes_labels(gca,{props{pri},[0 0 0]},{{'Neurons (%)'},[0 0 0]});
%     save_pdf(hf,mData.pdf_folder,sprintf('Distribution_%s',fileNames{pri}),600);
end
%%
while 1
    mzMIsC = arrayfun(@(x) nanmean(x{1}),zMIsC);
    [within,dvn,xlabels] = make_within_table({'Cond'},[3]);
    dataT = make_between_table({mzMIsC},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 0.25 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'Ar','ArL','Ar*'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.07 0.02 -0.4 -0.011])
    put_axes_labels(gca,{[],[0 0 0]},{props{pri},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells.pdf',fileNames{pri}),600);
end
end

%% Place Field Properties widths vs centers scatter plot and correlation
Rs = oD.Rs;
props1 = get_props_Rs(Rs,50);
respAnB = cell_list_op(props1.good_FR,props1.good_Gauss,'and');
props = {'Field Width (cm)','Field Center (cm)','Field FR (AU)'};
fileNames = {'PFW','PFC','PFFR'};
all_incr = [1 1 1];

while 1
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            R = Rs{rr,cc};
            [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
            W{rr,cc} = PWs(respAnB{rr,cc})';
            C{rr,cc} = centers(respAnB{rr,cc})';
            CW_corr(rr,cc) = corr(C{rr,cc},W{rr,cc});
        end
    end
    tcolors = mData.colors;
    hf = get_figure(8,[5 7 1.25 1]);hold on;
    an = 3; 
    for cn = 1:3
        scatter(C{an,cn},W{an,cn},2,tcolors{cn});
        [x,p] = polyfit(C{an,cn},W{an,cn},1);
        y_hat = polyval(x,C{an,cn});
        plot(C{an,cn},y_hat,'color',tcolors{cn});
    end
    format_axes(gca);
    changePosition(gca,[0.1 0.13 -0.25 -0.13]);
    put_axes_labels(gca,{'Field Center (cm)',[0 0 0]},{'Field Width (cm)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('scatter_%s','PFW_PFC'),600);

    [within,dvn,xlabels] = make_within_table({'Cond'},[3]);
    dataT = make_between_table({CW_corr},dvn);
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0.25 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1)-0.1 maxY]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.23 0.02 -0.4 -0.011])
    put_axes_labels(gca,{[],[0 0 0]},{{'Corr. Field Width','Vs Field Center'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells.pdf','CW_corr'),600);
end

%%
while 1
    Rs = oD.Rs;
    props1 = get_props_Rs(Rs,50);
    respAnB = cell_list_op(props1.good_FR,props1.good_Gauss,'and');
%     respAnB = props1.good_FR;
    var1 = 100*exec_fun_on_cell_mat(respAnB,'sum')./exec_fun_on_cell_mat(respAnB,'length');
    [within,dvn,xlabels] = make_within_table({'Cond'},[size(Rs,2)]);
    dataT = make_between_table({var1},dvn);
    ra = RMA(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    h(h==1) = 0;
    hf = get_figure(5,[8 7 1.75 1]);
    
    tcolors = colors;%[s.m;s.c;s.y];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
        set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);

    format_axes(gca);
    xticks = xdata; xticklabels = {'Ar-B1','Ar-B2','Ar-B3','ArL-B1','ArL-B2','ArL-B3','Ar*-B1','Ar*-B2','Ar*-B3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)

        changePosition(gca,[0.01 0.02 -0.03 -0.011])

    put_axes_labels(gca,{[],[0 0 0]},{'%',[0 -2 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('bar_graph_place_cells_belt.pdf'),600);
    ra.ranova
    ra.mauchly
    break
end


%% Place Field Properties on the belt
% Rs = oD.Rs;
si = [Ar_t_D ArL_t_D Ars_t_D];
Rs = o.Rs(:,si);mR = o.mR(:,si);
props1 = get_props_Rs(Rs,50);
% respAnB = cell_list_op(props1.good_FR,props1.good_Gauss,'and');
respAnB = props1.good_FR;
respCells1 = get_cell_list(respAnB,[1 -2 -3]);     respCells2 = get_cell_list(respAnB,[-1 2 -3]);     respCells3 = get_cell_list(respAnB,[-1 -2 3]);
all_respCells = {respCells1,respCells2,respCells3};
props = {'Field Width (cm)',{'Spatially Tuned','Cells (%)'},'Field FR (AU)'};
fileNames = {'PFW','PFC','PFFR'}; yspacing = [10 5 5]; 
pri = 2;
%%
while 1
    sRs = Rs(:,1:3);
    minBin = 0;     maxBin = 150;     incr = 50; % choosing more than 3 bins, give significant anova but not significant multcompare
    % three bins also make sense because LED in condition 2 comes ON at
    % 110cm
    bins = minBin:incr:maxBin;
    num_cells = []; PFW_belt = []; MFR_belt = [];
    for rr = 1:size(sRs,1)
        for cc = 1:size(sRs,2)
            R = sRs{rr,cc};
            [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
            resp = respAnB{rr,cc}; 
%             resp = all_respCells{cc}{rr};
            centers = centers(resp)'; PWs = PWs(resp)'; MFR = MFR(resp)';
            for ii = 1:(length(bins)-1)
                binS = bins(ii); binE = bins(ii+1);
                inds = find(centers >= binS & centers < binE);
                num_cells(rr,cc,ii) = 100*length(inds)/sum(resp);
                PFW_belt(rr,cc,ii) = mean(PWs(inds));
                MFR_belt(rr,cc,ii) = mean(MFR(inds));
            end
        end
    end
    mzMIsC = [];
    for rr = 1:size(sRs,1)
        switch pri
            case 1
                mzMIsC(rr,:) = reshape(squeeze(PFW_belt(rr,:,:))',1,(length(bins)-1)*size(sRs,2));                
            case 2
                mzMIsC(rr,:) = reshape(squeeze(num_cells(rr,:,:))',1,(length(bins)-1)*size(sRs,2));
            case 3
                mzMIsC(rr,:) = reshape(squeeze(MFR_belt(rr,:,:))',1,(length(bins)-1)*size(sRs,2));
        end
    end
    [within,dvn,xlabels] = make_within_table({'Cond','Bins'},[size(sRs,2),(length(bins)-1)]);
    dataT = make_between_table({mzMIsC},dvn);
    ra = RMA(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_Bins','hsd'},[1 1 1]);
    h(h==1) = 0;
    hf = get_figure(5,[8 7 1.75 1]);
    s = generate_shades(length(bins)-1);
    tcolors = colors;%[s.m;s.c;s.y];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',yspacing(pri),'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    if pri == 2
        set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);
    else
        set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);
    end
    format_axes(gca);
    xticks = xdata; xticklabels = {'Ar-B1','Ar-B2','Ar-B3','ArL-B1','ArL-B2','ArL-B3','Ar*-B1','Ar*-B2','Ar*-B3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    if pri == 2
        changePosition(gca,[0.11 0.02 -0.03 -0.011])
    else
        changePosition(gca,[0.01 0.02 -0.03 -0.011])
    end
    put_axes_labels(gca,{[],[0 0 0]},{props{pri},[0 -2 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells_belt.pdf',fileNames{pri}),600);
    ra.ranova
    %%
    break;
end


%% percent of spatially tuned cells on the belt
while 1
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Bins','hsd'},[1 0.25 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    s = generate_shades(length(bins)-1);
    tcolors = colors;%[s.m;s.c;s.y];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',15,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);
    format_axes(gca);
    xticks = xdata; xticklabels = xlabels;
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0 0.02 -0.35 -0.05])
    put_axes_labels(gca,{[],[0 0 0]},{props{pri},[0 -2 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells_belt_pooled.pdf',fileNames{pri}),600);
    ra.ranova
    ra.mauchly
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
    save_pdf(hf,mData.pdf_folder,sprintf('%s_cells.pdf',txtT{tC}),600);
    ra.ranova
%     ra.mauchly
    break;
end
end


%% Place cell emergence disruption stability - big combined test
while 1
    respCells1 = get_cell_list(respAnB,[1]); respCells2 = get_cell_list(respAnB,[2]);    respCells3 = get_cell_list(respAnB,[3]);
    respCells1o2 = get_cell_list(respAnB,[1;2]); respCells2o3 = get_cell_list(respAnB,[2;3]);    respCells1o3 = get_cell_list(respAnB,[1;3]);
    or_cells = get_cell_list(respAnB,[1;2;3])
    txtT = {'New','Disrupted'};
    percCells = [];
    % new place cells
    respCells12 = get_cell_list(respAnB,[-1 2]);    respCells23 = get_cell_list(respAnB,[-2 3]);    respCells13 = get_cell_list(respAnB,[-1 3]);
    for ii = 1:length(respCells1)
        percCells(ii,1) = sum(respCells12{ii})/length(or_cells{ii});
        percCells(ii,2) = sum(respCells23{ii})/length(or_cells{ii});
        percCells(ii,3) = sum(respCells13{ii})/length(or_cells{ii});
    end
    % disrupted place cells
    respCells12 = get_cell_list(respAnB,[1 -2]);    respCells23 = get_cell_list(respAnB,[2 -3]);    respCells13 = get_cell_list(respAnB,[1 -3]);
    for ii = 1:length(respCells1)
        percCells(ii,4) = sum(respCells12{ii})/length(or_cells{ii});%length(respCells1{ii});
        percCells(ii,5) = sum(respCells23{ii})/length(or_cells{ii});%length(respCells2{ii});
        percCells(ii,6) = sum(respCells13{ii})/length(or_cells{ii});%length(respCells3{ii});
    end
     % remained place cells
    respCells12 = get_cell_list(respAnB,[1 2]);    respCells23 = get_cell_list(respAnB,[2 3]);    respCells13 = get_cell_list(respAnB,[1 3]);
    for ii = 1:length(respCells1)
        percCells(ii,7) = sum(respCells12{ii})/length(or_cells{ii});%length(respCells1{ii});
        percCells(ii,8) = sum(respCells23{ii})/length(or_cells{ii});%length(respCells2{ii});
        percCells(ii,9) = sum(respCells13{ii})/length(or_cells{ii});%length(respCells3{ii});
    end
    

    [within,dvn,xlabels] = make_within_table({'Conds','Type'},[3 3]);
    dataT = make_between_table({percCells},dvn);
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
    hf = get_figure(5,[8 7 2.5 1]);
    tcolors = colors;%[s.m;s.c;s.y];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);
    format_axes(gca);
    xticks = xdata; xticklabels = xlabels;
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0 0.02 -0.15 -0.05])
    put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('dynamics_cells.pdf'),600);
    ra.ranova
    ra.mauchly
end

%% trial by trial comparison - difference in peak locations
while 1
    Rs = oD.Rs;
    props1 = get_props_Rs(Rs,50);
    gFR = props1.good_FR;
    pL_trials = props1.peak_locations_trials;
    vals = [];
    for rr = 1:size(gFR,1)
        one_row = [];
        for cc = 1:size(gFR,2)
            tgFR = gFR{rr,cc};
            tplT = pL_trials{rr,cc}; % peak locations of all trials (of all cells)
            mean_diff_tplT = nanmean(diff(tplT(tgFR,:),[],2));
            one_row = [one_row mean_diff_tplT];
        end
        vals(rr,:) = one_row;
    end
    [within,dvn,xlabels] = make_within_table({'Cond','TP'},[3,9]);
    dataT = make_between_table({vals},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_TP','hsd'},[1 1 1]);
    h(h==1) = 0;
%     s = generate_shades(3); tcolors = s.g;
    tcolors = colors;
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 3.25 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
    format_axes(gca);
    ylims = ylim;
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[ylims(1) maxY])
    xticks = [xdata(1:end)]; xticklabels = {'B1','B2','B3'};
%     set_title(gca,'Air-Belt Mismatch',[-0.35,95],5)
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.0 -0.011])
    
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Difference Peak','Location (cm)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Difference_Peak_Locations'),600);
    
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'TP','hsd'},[1 1 1]);%h(h==1) = 0;
%     s = generate_shades(3); tcolors = s.g;
    tcolors = colors;
     hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 3 2 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
    format_axes(gca);
    ylims = ylim;
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[ylims(1) maxY])
    xticks = [xdata(1:end)]; xticklabels = {'B1','B2','B3'};
%     set_title(gca,'Air-Belt Mismatch',[-0.35,95],5)
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.12 0.02 -0.15 -0.091])
    
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Difference Peak','Location (cm)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Difference_Peak_Locations_pooled'),600);
    break;
end

%% trial by trial comparison - place location gaussian fitting
% [resp_fractionCS,respAnB,OIC,mean_OICS,resp_ORCS,resp_OR_fractionCS,resp_ANDCS,resp_AND_fractionCS] = get_responsive_fraction(Rs);
% resp = get_cell_list(respAnB,[1;2;3]);
props1 = get_props_Rs(Rs,50);
gFR = props1.good_FR;
pL_trials = props1.peak_locations_trials;
out_tt_PC = find_trial_by_trial_Gauss_Fit_Results(Rs,mR,gFR);
[within,dvn,xlabels] = make_within_table({'Conds','DC'},[3,9]);
dataT1 = [];
for rr = 1:size(out_tt_PC.diff_centers_from_mean,1)
    thisD = [];
    for cc = 1:size(out_tt_PC.diff_centers_from_mean,2)
        thisDT = nanmean(out_tt_PC.diff_centers{rr,cc});
        thisD = [thisD thisDT];
    end
    dataT1(rr,:) = thisD;
end

dataT = make_between_table({dataT1},dvn);
ra = RMA(dataT,within,0.05);
[xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'DC','bonferroni'},[1 1 1]);
hf = get_figure(5,[8 7 9 1]);
% s = generate_shades(length(bins)-1);
tcolors = colors;%[s.m;s.c;s.y];
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
ylims = ylim;
set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]);
format_axes(gca);
xticks = xdata; xticklabels = xlabels;
set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
changePosition(gca,[0.05 0.02 -0.35 -0.05])
put_axes_labels(gca,{[],[0 0 0]},{'Jitter Centers (cm)',[0 0 0]});
save_pdf(hf,mData.pdf_folder,sprintf('jitter_centers.pdf'),600);
ra.ranova
ra.mauchly