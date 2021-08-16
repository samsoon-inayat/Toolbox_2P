function belt_place_cells_figure

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei = evalin('base','ei'); 

selContexts = [3 4 5];
rasterNames = {'airD','airD','airD'};
% Rs = get_rasters_data(ei_11_15,selContexts,rasterNames);
RsA = get_rasters_data(ei,selContexts,rasterNames);
% Rs = filterRasters(Rs);
RsA = find_responsive_rasters(RsA,1:10);
% [resp_fractionCSA,resp_valsCSA,OICA,mean_OICSA,resp_ORCSA,resp_OR_fractionCSA,resp_ANDCSA,resp_AND_fractionCSA] = get_responsive_fraction(RsA);
respA = get_responsive_cells(RsA);
resp_FRA = get_responsive_fraction_FR(RsA);
mRA = calc_mean_rasters(RsA,1:10);

selContexts = [3 4 5];
rasterNames = {'beltD','beltD','beltD'};
% Rs = get_rasters_data(ei_11_15,selContexts,rasterNames);
Rs = get_rasters_data(ei,selContexts,rasterNames);
% Rs = filterRasters(Rs);
Rs = find_responsive_rasters(Rs,1:10);
% [resp_fractionCS,resp_valsCS,OIC,mean_OICS,resp_ORCS,resp_OR_fractionCS,resp_ANDCS,resp_AND_fractionCS] = get_responsive_fraction(Rs);
respB = get_responsive_cells(Rs);
resp_FR = get_responsive_fraction_FR(Rs);
mR = calc_mean_rasters(Rs,1:10);

respBnA = sep_cell_list(respB.vals,respA.vals);
respBnAandFR = cell_list_op(respBnA,resp_FR,'and');
sum_respBnAandFR = exec_fun_on_cell_mat(respBnAandFR,'sum');
length_respBnAandFR = exec_fun_on_cell_mat(respBnAandFR,'length');
frac_respBnAandFR = sum_respBnAandFR./length_respBnAandFR;
n = 0;
%% find correlations across conditions
resp = get_cell_list(respBnA,[1;2;3],0);
[CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,1,0);

outRemap = find_population_vector_corr_remap(Rs,mR,respBnA);
n = 0;

%% find trial by trial comparison - correlations
trials = mat2cell([1:10]',ones(size([1:10]')));
resp = get_cell_list(respBnA,[]);
out1 = find_population_vector_corr_remap_trials(Rs(:,1),resp,trials);
out2 = find_population_vector_corr_remap_trials(Rs(:,2),resp,trials);
out3 = find_population_vector_corr_remap_trials(Rs(:,3),resp,trials);
% save(fullfile(mData.pd_folder,'air_place_cells_figure.mat'),'out1','out2','out3');

n = 0;
%%
% an = 1; cn = 1;
% plotRasters_simplest(Rs{an,cn})
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
    [Y,E] = discretize(1:49,3);
    all_speeds = [];
    for cn = 1:3
        mean_speed_over_trials = [];
        aThisSpeed = [];
        for an = 1:size(Rs,1)
            thisSpeed = nanmean(Rs{an,cn}.speed);
            mean_speed_over_trials(an,:) = thisSpeed;
            for bb = 1:3
                aThisSpeed(an,bb) = nanmean(thisSpeed(Y==bb));
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
save_pdf(ff.hf,mData.pdf_folder,sprintf('speeds_345_B'),600);
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
    set_title(gca,'BFoR',[1.3,35],5)
    save_pdf(hf,mData.pdf_folder,sprintf('avg_speed_anova_B.pdf'),600);
    ra.ranova
    ra.mauchly
end
%% Show sample rasters
% an = 5; cn = 2;
% plotRasters_simplest(Rs{an,cn})
if 1
   an = 3; cn = 1;
    % plotRasters_simplest(Rs{an,cn})
    % find(resp_valsC{an}(:,cn));
    ff = makeFigureRowsCols(2021,[0.5 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
        [-75 -475]);
    set(gcf,'color','w'); set(gcf,'Position',[10 4 3.25 1]);
    ff = sample_rasters(Rs{an,cn},[191 11 96 41],ff);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersD'),600);
   
end
%% population vector and correlation single animal
if 1
    an = 1;
    ff = makeFigureRowsCols(106,[1 0.5 6 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.07 0.1],'widthHeightAdjustment',...
        [0.01 -60]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 5 3.25 2]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,respBnA,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[-0.1 1],[]);
    set_obj(ff,{'FontWeight','Normal','FontSize',6,'LineWidth',0.25});
    ht = get_obj(ff,'title'); hyl = get_obj(ff,'ylabel'); changePosition(hyl(1,1),[4 0 0]);
    set_obj(ht,'String',{'Pop. Activity','Pop. Activity','Pop. Activity';'Pop. Correlation','Pop.Correlation','Pop.Correlation'});
    set_obj(ht,{'FontSize',5,'FontWeight','Normal'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('population_vector_corr_B.pdf'),600);
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
    save_pdf(ff.hf,mData.pdf_folder,sprintf('average_population_vector_corrD_B.pdf'),600);
end

%% Percentage of Responsive Cells
if 1
    [within,dvn,xlabels] = make_within_table({'FR','Cond'},[2,3]);
    dataT = make_between_table({resp_fractionCSA*100,resp_fractionCSB*100},dvn);
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
    any_mean = mean(100*resp_OR_fractionCSA);    any_sem = std(100*resp_OR_fractionCSA)/sqrt(5);
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
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Condition_by_TrialPairs','bonferroni'},[1 0.5 1]);
    hollowsep = 19;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 3.75 1],'color','w');
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
%% Spatial correlation trails with mean of trials (mean bar graph)
if 1
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Condition_by_TrialPairs','bonferroni'},[1 0.5 1]);
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
    
    OI_E = get_overlap_index(respBnA);
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
                    zMIsC{rr,cc} = PWs(respBnA{rr,cc})';
                case 2
                    zMIsC{rr,cc} = centers(respBnA{rr,cc})';
                case 3
                    zMIsC{rr,cc} = MFR(respBnA{rr,cc})';
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
            W{rr,cc} = PWs(respBnA{rr,cc})';
            C{rr,cc} = centers(respBnA{rr,cc})';
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
props = {'Field Width (cm)',{'Spatially Tuned','Cells (%)'},'Field FR (AU)'};
fileNames = {'PFW','PFC','PFFR'}; yspacing = [10 5 5]; 
pri = 1;
Air = 0;
if 1
    sRs = Rs;
    minBin = 0;     maxBin = 150;     incr = 50; % choosing more than 3 bins, give significant anova but not significant multcompare
    % three bins also make sense because LED in condition 2 comes ON at
    % 110cm
    bins = minBin:incr:maxBin;
    num_cells = []; PFW_belt = []; MFR_belt = [];
    for rr = 1:size(sRs,1)
        for cc = 1:size(sRs,2)
            R = sRs{rr,cc};
            [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
            resp = respBnA{rr,cc};
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
    [within,dvn,xlabels] = make_within_table({'Cond','Bins'},[3,3]);
    dataT = make_between_table({mzMIsC},dvn);
    ra = RMA(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_Bins','bonferroni'},[1 0.25 1]);
    hf = get_figure(5,[8 7 3.49 1]);
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
    put_axes_labels(gca,{[],[0 0 0]},{props{pri},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells_belt.pdf',fileNames{pri}),600);
    ra.ranova
    ra.mauchly
end


%% percent of spatially tuned cells on the belt
if 1
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','bonferroni'},[1 0.25 1]);
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
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells_belt_pooled_B.pdf',fileNames{pri}),600);
    ra.ranova
    ra.mauchly
end

%% Place cell emergence disruption stability
for tC = 1:4
if 1
    txtT = {'Unique','New','Disrupted','Stable'};
    or_cells = get_cell_list(respBnA,[1;2;3]);
    percCells = []; 
    switch tC
        case 1 % unique place cells
            respCells1 = get_cell_list(respBnA,[1 -2 -3]);    respCells2 = get_cell_list(respBnA,[-1 2 -3]);    respCells3 = get_cell_list(respBnA,[-1 -2 3]);
        case 2 % new place cells
            respCells1 = get_cell_list(respBnA,[-1 2]);    respCells2 = get_cell_list(respBnA,[-2 3]);    respCells3 = get_cell_list(respBnA,[-1 3]);
        case 3 % disrupted place cells
            respCells1 = get_cell_list(respBnA,[1 -2]);    respCells2 = get_cell_list(respBnA,[2 -3]);    respCells3 = get_cell_list(respBnA,[1 -3]);
        case 4 % remained place cells
            respCells1 = get_cell_list(respBnA,[1 2]);    respCells2 = get_cell_list(respBnA,[2 3]);    respCells3 = get_cell_list(respBnA,[1 3]);
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
    respCells1 = get_cell_list(respBnA,[1]); respCells2 = get_cell_list(respBnA,[2]);    respCells3 = get_cell_list(respBnA,[3]);
    respCells1o2 = get_cell_list(respBnA,[1;2]); respCells2o3 = get_cell_list(respBnA,[2;3]);    respCells1o3 = get_cell_list(respBnA,[1;3]);
    or_cells = get_cell_list(respBnA,[1;2;3])
    txtT = {'New','Disrupted'};
    percCells = [];
    % new place cells
    respCells12 = get_cell_list(respBnA,[-1 2]);    respCells23 = get_cell_list(respBnA,[-2 3]);    respCells13 = get_cell_list(respBnA,[-1 3]);
    for ii = 1:length(respCells1)
        percCells(ii,1) = sum(respCells12{ii})/length(or_cells{ii});
        percCells(ii,2) = sum(respCells23{ii})/length(or_cells{ii});
        percCells(ii,3) = sum(respCells13{ii})/length(or_cells{ii});
    end
    % disrupted place cells
    respCells12 = get_cell_list(respBnA,[1 -2]);    respCells23 = get_cell_list(respBnA,[2 -3]);    respCells13 = get_cell_list(respBnA,[1 -3]);
    for ii = 1:length(respCells1)
        percCells(ii,4) = sum(respCells12{ii})/length(or_cells{ii});%length(respCells1{ii});
        percCells(ii,5) = sum(respCells23{ii})/length(or_cells{ii});%length(respCells2{ii});
        percCells(ii,6) = sum(respCells13{ii})/length(or_cells{ii});%length(respCells3{ii});
    end
     % remained place cells
    respCells12 = get_cell_list(respBnA,[1 2]);    respCells23 = get_cell_list(respBnA,[2 3]);    respCells13 = get_cell_list(respBnA,[1 3]);
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
[resp_fractionCS,respBnA,OIC,mean_OICS,resp_ORCS,resp_OR_fractionCS,resp_ANDCS,resp_AND_fractionCS] = get_responsive_fraction(Rs);
resp = get_cell_list(respBnA,[1;2;3]);
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