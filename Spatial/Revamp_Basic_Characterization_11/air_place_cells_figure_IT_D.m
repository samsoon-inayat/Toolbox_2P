function air_place_cells_figure

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei = evalin('base','ei'); 


selContexts = [3 4 5];
rasterNames = {'airIT','airIT','airIT'};
% Rs = get_rasters_data(ei_11_15,selContexts,rasterNames);
RsT = get_rasters_data(ei,selContexts,rasterNames);
% RsT = filterRasters(RsT);
RsT = find_responsive_rasters(RsT,1:9);
respT = get_responsive_cells(RsT);
resp_FRT = get_responsive_fraction_FR(RsT);
mRT = calc_mean_rasters(RsT,1:9);

selContexts = [3 4 5];
rasterNames = {'airID','airID','airID'};
% Rs = get_rasters_data(ei_11_15,selContexts,rasterNames);
Rs = get_rasters_data(ei,selContexts,rasterNames);
% Rs = do_dtw(Rs);
% Rs = filterRasters(Rs);
Rs = find_responsive_rasters(Rs,1:9);
respA = get_responsive_cells(Rs);
resp_FR = get_responsive_fraction_FR(Rs);
mR = calc_mean_rasters(Rs,1:9);

selContexts = [3 4 5];
% rasterNames = {'beltD','beltD','beltD'};
% % Rs = get_rasters_data(ei_11_15,selContexts,rasterNames);
% RsB = get_rasters_data(ei,selContexts,rasterNames);
% % Rs = filterRasters(Rs);
% RsB = find_responsive_rasters(RsB,1:10);
% respB = get_responsive_cells(RsB);
% resp_FRB = get_responsive_fraction_FR(RsB);
% mRB = calc_mean_rasters(RsB,1:10);
% 
% respAnB = sep_cell_list(respA.vals,respB.vals);
n = 0;
%% find correlations across conditions
resp = resp_ORCS;
[CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,1,0);

outRemap = find_population_vector_corr_remap(Rs,mR,resp_ORCS);
n = 0;

%% find trial by trial comparison - correlations
trials = mat2cell([1:10]',ones(size([1:10]')));
resp = get_cell_list(resp_valsCS,[]);
out1 = find_population_vector_corr_remap_trials(Rs(:,1),resp,trials);
out2 = find_population_vector_corr_remap_trials(Rs(:,2),resp,trials);
out3 = find_population_vector_corr_remap_trials(Rs(:,3),resp,trials);
% save(fullfile(mData.pd_folder,'air_place_cells_figure.mat'),'out1','out2','out3');
n = 0;

%% Speed Figure
if 1
    for ii = 1:length(ei)
        b1 = ei{ii}.b;
        for jj = 1:10
            alds(ii,jj) = b1.dist(b1.stim_r(jj+10)) - b1.dist(b1.air_puff_r(jj+20));
        end
    end
    ald = round(mean(alds(:)));
     ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[-0.01 0.13],'widthHeightAdjustment',...
        [-45 -250]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 5 3.25 1]);
    for cn = 1:3
        for an = 1:size(Rs,1)
            R = Rs{an,cn};
            szSp(an,cn) = size(R.dist,2);
        end
    end
    
    [Y,E] = discretize(1:49,3);
    all_speeds = [];
    for cn = 1:3
        mean_speed_over_trials = [];
        aThisSpeed = [];
        for an = 1:size(Rs,1)
            thisSpeed = nanmean(Rs{an,cn}.speed);
            mean_speed_over_trials(an,:) = thisSpeed;
            for bb = 1:3
                try
                    aThisSpeed(an,bb) = nanmean(thisSpeed(Y==bb));
                catch
                    aThisSpeed(an,bb) = NaN;
                end
            end
        end
        all_speeds = [all_speeds aThisSpeed];
        axes(ff.h_axes(1,cn));
        hold on;
        xs = Rs{1,2}.xs; N = length(xs);
        mspeed = mean(mean_speed_over_trials(:,1:N)); semspeed = std(mean_speed_over_trials(:,1:N))/sqrt(5);
        plot(xs,mspeed);
        shadedErrorBar(xs,mspeed,semspeed);
        changePosition(gca,[0.1 0.15 -0.05 -0.15]);
        if cn == 1
            put_axes_labels(gca,{'',[0 0 0]},{{'Speed (cm/sec)'},[0 -1 0]});
        end
        plot([ald ald],[0 30],'m:','linewidth',0.5);
        plot([50 50],[0 30],'k--','linewidth',0.25);
        plot([100 100],[0 30],'k--','linewidth',0.25);
         bTxt = {'B1','B2','B3'}; xbTxt = [25 75 125]-7; ybTxt = 31;
        for ii = 1:length(bTxt)
            text(xbTxt(ii),ybTxt,bTxt{ii},'FontSize',5);
        end
        xlabel('Position (cm)');
        ylim([0 30]);
        box off;
        set(gca,'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('speeds_345'),600);
%%
    
    [within,dvn,xlabels] = make_within_table({'Conds','Bins'},[3,3]);
    dataT = make_between_table({all_speeds},dvn);
    ra = RMA(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Bins','bonferroni'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors(5:9);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'B1','B2','B3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.35 -0.05])
    put_axes_labels(gca,{[],[0 0 0]},{'Avg. Speed (cm/sec)',[0 0 0]});
    set_title(gca,'AFoR',[1.3,39],5)
    save_pdf(hf,mData.pdf_folder,sprintf('avg_speed_anova.pdf'),600);
    ra.ranova
    ra.mauchly
end
%% Air Belt mismatch
if 1
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
end
%% Show sample rasters
% an = 5; cn = 2;
% plotRasters_simplest(Rs{an,cn})
if 1
   an = 3; cn = 1;
    % plotRasters_simplest(Rs{an,cn})
    % find(resp_valsC{an}(:,cn));
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
        [-75 -475]);
    set(gcf,'color','w'); set(gcf,'Position',[10 4 3.25 1]);
    ff = sample_rasters(Rs{an,cn},[191 11 96 41],ff);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersD'),600);
    
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
        [-75 -475]);
    set(gcf,'color','w'); set(gcf,'Position',[10 4 3.25 1]);
    ff = sample_rasters(RsT{an,cn},[191 11 96 41],ff);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersDT'),600);
    
end
%% Mutual Information Time versus Distance Distributions

if 1
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

%% population vector and correlation single animal
if 1
    an = 1;
    ff = makeFigureRowsCols(106,[1 0.5 6 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.07 0.1],'widthHeightAdjustment',...
        [0.01 -60]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 5 3.25 2]);
%     [resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC,resp_exc_inh] = get_responsive_fraction(Rs);
%     resp = get_cell_list(resp_valsC,[]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,respAnB,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[-0.1 1],[]);
    set_obj(ff,{'FontWeight','Normal','FontSize',6,'LineWidth',0.25});
    ht = get_obj(ff,'title'); hyl = get_obj(ff,'ylabel'); changePosition(hyl(1,1),[4 0 0]);
    set_obj(ht,'String',{'Pop. Activity','Pop. Activity','Pop. Activity';'Pop. Correlation','Pop.Correlation','Pop.Correlation'});
    set_obj(ht,{'FontSize',5,'FontWeight','Normal'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_population_vector_corr.pdf'),600);
end
%% average correlation of all animals
if 1
    ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.075 0.2],'widthHeightAdjustment',...
        [0.01 -250]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 5 3.25 1]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[-0.1 1],[]);
    set_obj(ff,{'FontWeight','Normal','FontSize',6,'LineWidth',0.25});
    ht = get_obj(ff,'title'); 
    set_obj(ht,'String',{'Avg. Pop. Correlation','Avg. Pop. Correlation','Avg. Pop. Correlation'});
    set_obj(ht,{'FontSize',5,'FontWeight','Normal'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_average_population_vector_corrD.pdf'),600);
end

%% population vector and correlation single animal TRIAL WISE trial wise
% For belt FoR, exclude cells that are also PCs for Air FoR
if 1
    an = 5;
    ff = makeFigureRowsCols(106,[1 0.5 6 1],'RowsCols',[2 9],...
        'spaceRowsCols',[0 -0.07],'rightUpShifts',[0.03 0.1],'widthHeightAdjustment',...
        [+60 -90]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 5 13 5]);
    tR = Rs(an,:)'; [tR,tmR] = separate_into_trials(tR,{1;2;3;4;5;6;7;8;9});
    respAsB_T = respA.vals(an,:)';respAsB_T = repmat(respAsB_T,1,9);
    cn = 1;
    [CRc,aCRc,mRR] = find_population_vector_corr(tR,tmR,respAsB_T,3);
    ff = show_population_vector_and_corr(mData,ff,tR(cn,:),mRR(cn,:),CRc(cn,:),[],[]);
    
    set_obj(ff,{'FontWeight','Normal','FontSize',6,'LineWidth',0.25});
    ht = get_obj(ff,'title'); hyl = get_obj(ff,'ylabel'); changePosition(hyl(1,1),[4 0 0]);
%     colormap jet
%     set_obj(ht,'String',{'Pop. Activity','Pop. Activity','Pop. Activity','Pop. Activity','Pop. Activity','Pop. Activity';'Pop. Correlation','Pop.Correlation','Pop.Correlation','Pop. Correlation','Pop.Correlation','Pop.Correlation'});
    set_obj(ht,{'FontSize',5,'FontWeight','Normal'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_population_vector_corr_trial_wise.pdf'),600);
end

%% Percentage of Responsive Cells
if 1
    [within,dvn,xlabels] = make_within_table({'FR','Cond'},[2,3]);
    dataT = make_between_table({resp_fractionCS*100,resp_fractionCSB*100},dvn);
    ra = RMA(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR_by_Cond','bonferroni'},[1 1 1]);
%     s = generate_shades(3); tcolors = s.g;
    tcolors = colors;%mData.colors(1:3);
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
    pmchar=char(177); any_text = sprintf('%.0f%c%.0f%%',any_mean,pmchar,any_sem); text(0.75,20,any_text,'FontSize',6);
    any_mean = mean(100*resp_OR_fractionCSB);    any_sem = std(100*resp_OR_fractionCSB)/sqrt(5);
    pmchar=char(177); any_text = sprintf('%.0f%c%.0f%%',any_mean,pmchar,any_sem); text(4.75,20,any_text,'FontSize',6);
    changePosition(gca,[0.17 0.02 -0.1 -0.011])
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('percentage_PCs_responsive'),600);
end

%% Percentage of Responsive Cells Unique in Conditions
if 1
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
end
%% Spatial correlation between adjacent trails
if 1
    [within,dvn,xlabels] = make_within_table({'Condition','TrialPairs'},[3,9]);
    var1 = arrayfun(@(x) mean(x{1}),out1.adj_SP_corr_diag); 
    var2 = arrayfun(@(x) mean(x{1}),out2.adj_SP_corr_diag);
    var3 = arrayfun(@(x) mean(x{1}),out3.adj_SP_corr_diag);
    between = make_between_table({var1,var2,var3},dvn);
    ra = RMA(between,within);
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
%% Spatial correlation between adjacent trails (taking mean) 
if 1
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'TrialPairs','bonferroni'},[1 0.5 1]);
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
end

%% Spatial correlation trails with mean of trials
if 1
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
if 1
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
if 1
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
%% Spatial correlation between conditions
if 1
    [within,dvn,xlabels] = make_within_table({'Cond'},2);
    var_CE = arrayfun(@(x) mean(x{1}),outRemap.adj_SP_corr_diag);
    dataT = make_between_table({var_CE},dvn);
    ra = repeatedMeasuresAnova(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0 1]);
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors ={colors{7};colors{8};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.02,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C3-C4','C4-C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     set(hbs(2),'facecolor','none','edgecolor',tcolors{2});
    xtickangle(45);
%     changePosition(gca,[0.2 0.03 -0.55 -0.05]);
    changePosition(gca,[0.1 0.03 -0.4 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatial Correlation'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'air_place_cell corr remap',600);
end

%% Delta FR score (RAte remap) between conditions
if 1
    [within,dvn,xlabels] = make_within_table({'Cond'},2);
    var_CE = arrayfun(@(x) nanmean(x{1}),outRemap.adj_RR_SP);
    dataT = make_between_table({var_CE},dvn);
    ra = repeatedMeasuresAnova(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0 1]);
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors ={colors{7};colors{8};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.02,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+0.1],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C3-C4','C4-C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     set(hbs(2),'facecolor','none','edgecolor',tcolors{2});
    xtickangle(45);
%     changePosition(gca,[0.2 0.03 -0.55 -0.05]);
    changePosition(gca,[0.1 0.03 -0.4 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'\Delta FR Score'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'air_place_cell rate remap (DFR)',600);
end

%% Popuation Vector correlation between conditions
if 1
    [within,dvn,xlabels] = make_within_table({'Cond'},2);
    var_CE = arrayfun(@(x) nanmean(x{1}),outRemap.adj_PV_corr_diag);
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
    changePosition(gca,[0.2 0.03 -0.4 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Pop. Vec.','Correlation'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'air_place_cell PV corr',600);
end

%%
%% overlap RM bar graph
if 1
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    
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
props = {'Field Width (cm)','Field Center (cm)','Field FR (AU)'};
fileNames = {'PFW','PFC','PFFR'};
all_incr = [1 1 1];
% legs_v = 
for pri = 1%:3;
if 1
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

    tcolors = mData.colors;
    distD = zMIsC;
    [distDo,allVals] = getAveragesAndAllValues(distD);
    minBin = min(allVals);
    maxBin = max(allVals);
    incr = all_incr(pri); %maxBin =
    hf = get_figure(8,[5 7 1.25 1]);
    [ha,hb,~,bins] = plotAverageDistributions(distD,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
    format_axes(gca);
    if pri == 1
        legs = {'C3','C4','C3''',[50 3 70 5]}; putLegend(gca,legs,tcolors);
    end
    changePosition(gca,[0.1 0.13 -0.25 -0.13]);
    put_axes_labels(gca,{props{pri},[0 0 0]},{{'Neurons (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_%s',fileNames{pri}),600);
end

if 1
    mzMIsC = arrayfun(@(x) nanmean(x{1}),zMIsC);
    [within,dvn,xlabels] = make_within_table({'Cond'},[3]);
    dataT = make_between_table({mzMIsC},dvn);
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0.25 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.07 0.02 -0.4 -0.011])
    put_axes_labels(gca,{[],[0 0 0]},{props{pri},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells.pdf',fileNames{pri}),600);
end
end

%% Place Field Properties widths vs centers scatter plot and correlation
props = {'Field Width (cm)','Field Center (cm)','Field FR (AU)'};
fileNames = {'PFW','PFC','PFFR'};
all_incr = [1 1 1];

if 1
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

%% Place Field Properties on the belt
respCells1 = get_cell_list(respAnB,[1 -2 -3]);     respCells2 = get_cell_list(respAnB,[-1 2 -3]);     respCells3 = get_cell_list(respAnB,[-1 -2 3]);
all_respCells = {respCells1,respCells2,respCells3};
props = {'Field Width (cm)',{'Spatially Tuned','Cells (%)'},'Field FR (AU)'};
fileNames = {'PFW','PFC','PFFR'}; yspacing = [10 5 5]; 
pri = 2;
if 1
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
                num_cells(rr,cc,ii) = 100*length(inds)/length(resp);
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
    ra = RMA(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_Bins','bonferroni'},[1 0.25 1]);
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
    xticks = xdata; xticklabels = {'C3-B1','C3-B2','C3-B3','C4-B1','C4-B2','C4-B3','C3''-B1','C3''-B2','C3''-B3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    if pri == 2
        changePosition(gca,[0.11 0.02 -0.03 -0.011])
    else
        changePosition(gca,[0.01 0.02 -0.03 -0.011])
    end
    put_axes_labels(gca,{[],[0 0 0]},{props{pri},[0 -2 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells_belt.pdf',fileNames{pri}),600);
    ra.ranova
    ra.mauchly
end


%% percent of spatially tuned cells on the belt
if 1
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Bins','bonferroni'},[1 0.25 1]);
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
for tC = 1:4
if 1
    txtT = {'Unique','New','Disrupted','Stable'};
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

    for ii = 1:length(respCells1)
        percCells(ii,1) = 100*sum(respCells1{ii})/length(or_cells{ii});%length(respCells1{ii});
        percCells(ii,2) = 100*sum(respCells2{ii})/length(or_cells{ii});%length(respCells2{ii});
        percCells(ii,3) = 100*sum(respCells3{ii})/length(or_cells{ii});%length(respCells3{ii});
    end

    [within,dvn,xlabels] = make_within_table({'Conds'},[3]);
    if tC > 1
       xlabels = {'C3-C4','C4-C3''','C3-C3'''};
    end
    dataT = make_between_table({percCells},dvn);
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0.25 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    tcolors = colors;%[s.m;s.c;s.y];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',15,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);
    format_axes(gca);
    xticks = xdata; xticklabels = xlabels;
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.1 0.02 -0.35 -0.05])
    put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    text(0.75,maxY,txtT{tC},'FontSize',6)
    save_pdf(hf,mData.pdf_folder,sprintf('%s_cells.pdf',txtT{tC}),600);
    ra.ranova
    ra.mauchly
end
end


%% Place cell emergence disruption stability - big combined test
if 1
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

%% trial by trial comparison - place location gaussian fitting
[resp_fractionCS,respAnB,OIC,mean_OICS,resp_ORCS,resp_OR_fractionCS,resp_ANDCS,resp_AND_fractionCS] = get_responsive_fraction(Rs);
resp = get_cell_list(respAnB,[1;2;3]);
out_tt_PC = find_trial_by_trial_Gauss_Fit_Results(Rs,mR,resp);
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