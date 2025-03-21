function conjunctive_representation

%% Load Data
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei = evalin('base','eiA'); 
    
    selContexts = [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4];
    rasterNames = {'airD','airD','airD','airD','airIT','airIT','airIT','airIT','airT','airT','airT','airT','airID','airID','airID','airID','beltD','beltD','beltD','beltT','beltT','beltT','beltT','beltT'};
    rasterNamesTxt = {'Ar-t-D','ArC-t-D','ArCB-t-D','ArB-t-D','Ar-i-T','ArC-i-T','ArCB-i-T','ArB-i-T','Ar-t-T','ArC-t-T','ArCB-t-T','ArB-t-T','Ar-i-D','ArC-i-D','ArCB-i-D','ArB-i-D'};

    o = get_data(ei,selContexts,rasterNames);
    for ii = 1:length(selContexts)
        all_xl{ii} = sprintf('%c%c-%d',rasterNames{ii}(1),rasterNames{ii}(end),selContexts(ii));
    end
    Ar_t_D = 1; ArC_t_D = 2; ArCB_t_D = 3; ArB_t_D = 4;     Ar_i_T = 5; ArC_i_T = 6; ArCB_i_T = 7; ArB_i_T = 8;
    Ar_t_T = 9; ArC_t_T = 10; ArCB_t_T = 11; ArB_t_T = 12;     Ar_i_D = 13; ArC_i_D = 14; ArCB_i_D = 15; ArB_i_D = 16;
    break
end
n = 0;
%% trial formation
while 1
    filename = fullfile(mData.pd_folder,sprintf('%s_trials_formation',mfilename));
    if 1
        si = [Lb Ab Ar_t_D Ar_i_T];        Rs = o.Rs(:,si);
        trials = mat2cell([1:10]',ones(size([1:10]')));
        props1 = get_props_Rs(Rs,50);
        parfor ii = 1:size(Rs,2)
            outTrials{ii} = find_population_vector_corr_remap_trials(Rs(:,ii),props1.good_FR(:,ii),trials);
%             outTrials_C{ii} = find_population_vector_trial_to_trial_corr(Rs(:,ii),props1.good_FR(:,ii));
        end
        outTrials_tuned = [];
%         parfor ii = 1:size(Rs,2)
%             outTrials_tuned{ii} = find_population_vector_corr_remap_trials(Rs(:,ii),props1.good_FR_and_tuned(:,ii),trials);
%         end
        n = 0;
        save(filename,'outTrials','trials');
    else
        si = [Lb Ab Ar_t_D Ar_i_T];        Rs = o.Rs(:,si);
        trials = mat2cell([1:10]',ones(size([1:10]'))); props1 = get_props_Rs(Rs,50);
        load(filename);
    end
    break;
end
    n = 0;
    %% avg correlation of across trials
while 1
    meancorr_trials = [];
    for ii = 1:11
        toutTrials = outTrials{ii};
        meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(toutTrials.adj_SP_corr_diag,'nanmean')];
    end
    [within,dvn,xlabels] = make_within_table({'Cond','TP'},[11,9]);
    dataT = make_between_table({meancorr_trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
%     xdata = xdataG;
    ptab = 0;
    if ptab h(h==1) = 0; end
    hf = get_figure(5,[8 7 1.75 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    
    if ptab
    [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
        'BaseValue',0.01,'xdata',xdatag,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.5);
    else
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'yspacing',0.1);
    end
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    txl = rasterNamesTxt(si); 
    xticks = xdata; xticklabels = txl;
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 0.1]); xtickangle(45);
    changePosition(gca,[0.02 -0.01 0.06 0.01]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. Correlation',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('avg correlation'),600);
    %%
break;
end
%% Show sample rasters
while 1
    Rs = o.Rs;
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
    break;
end

%% compare the zMIs
while 1
    ntrials = 50;
    si = si_seqG;%si_no_brake_dist; 
    Rs = o.Rs(:,si)
    props1 = get_props_Rs(Rs,ntrials);
    good_FR = props1.good_FR;
    all_zMIs = props1.zMI;
    zMIs = [];
    for rr = 1:size(all_zMIs,1)
        for cc = 1:size(all_zMIs,2)
            resp = good_FR{rr,cc};
            tzmis = all_zMIs{rr,cc};
            zMIs(rr,cc) = nanmean(tzmis(resp));
        end
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},[11]);
    dataT = make_between_table({zMIs},dvn);
    ra = RMA(dataT,within);
    ra.ranova
%     ra.mauchly
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = xdataG;
    ptab = 0;
    if ptab h(h==1) = 0; end
    hf = get_figure(5,[8 7 1.75 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    

    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'yspacing',0.5);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rasterNamesTxt(si);
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 round(maxY/2) round(maxY)]); xtickangle(45);
    if ptab
    changePosition(gca,[0.08 0.01 0.0 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. zMI',[-1 3 0]});
    ha = gca; ptable = extras.pvalsTable;
    display_p_table(ha,hbs,[0 -0.07 0 0.9],ptable);
    else
    changePosition(gca,[0.01 -0.01 0.06 0.01]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. zMI',[0 0 0]});
    end
    save_pdf(hf,mData.pdf_folder,sprintf('zMIs_good_FR_all_Conditions.pdf'),600);
    %%
    break;
end

%% compare the peak locations of dist MIs versus time MIs i.e., which cells have larger MI for dist or time
while 1
    props1 = get_props_Rs(o.Rs,3);
    si = si_no_brake_dist;
    dMI = prop_op(props1.zMI(:,si_no_brake_dist),props1.zMI(:,si_no_brake_time),0);
    good_FR = props1.good_FR(:,si);
    DgT = dMI.resp_D_g_T;
    TgD = dMI.resp_T_g_D;
    all_zMIs = props1.peak_locations(:,si);
    zMIs = []; zMIsT = [];
    for rr = 1:size(all_zMIs,1)
        for cc = 1:size(all_zMIs,2)
            resp = good_FR{rr,cc} & DgT{rr,cc};
            respT = good_FR{rr,cc} & TgD{rr,cc};
            tzmis = all_zMIs{rr,cc};
            zMIs(rr,cc) = nanmean(tzmis(resp));
            zMIsT(rr,cc) = nanmean(tzmis(respT));
        end
    end
    azMIs = [zMIs zMIsT];
    [within,dvn,xlabels] = make_within_table({'DT','Cond','T'},[2,3,2]);
    dataT = make_between_table({azMIs},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    ra.mauchly
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'DT_by_T','bonferroni'},[1 1 1]);
%     h(h==1) = 0;
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rasterNames(si);
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[-0.01 -0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Active Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('zMIs_good_FR_all_Conditions.pdf'),600);
    %%
    break;
end

%% distribution of trials in which cells fired
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    si = [Lb_T ArL_L_T Lbs_T Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
    props1 = get_props_Rs(o.Rs,ntrials);
    good_FR = props1.N_Resp_Trials(:,si);
    %%
    hf = get_figure(8,[5 7 2.25 1.5]);hold on;
    mData.colors = [mData.colors;mData.colors];
    
   for cn = 1:(size(good_FR,2))
        distD = good_FR(:,cn);
        [distDo,allVals] = getAveragesAndAllValues(distD);
        minBin = min(allVals);
        maxBin = max(allVals);
        incr = 20;
        tcolors = mData.colors(cn,:);
        [ha,hb,~,bins] = plotAverageDistributions(distD,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'pdf_or_cdf','pdf');
    end
    format_axes(gca);
%     changePosition(gca,[0.1 0.13 -0.25 -0.13]);
%     put_axes_labels(gca,{props{pri},[0 0 0]},{{'Neurons (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_N_trials_resp_%d',cn),600);
    
    %%
    mean_N_trials_resp = exec_fun_on_cell_mat(good_FR,'mean');
    [within,dvn,xlabels] = make_within_table({'Cond'},[length(si)]);
    dataT = make_between_table({mean_N_trials_resp},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    
     %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = make_xdata([3 2 3 3],[1 1.5]);
    h(h==1) = 0;
    hf = get_figure(5,[10 7 2.2 1.5]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.5);
    
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rNaG;
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 20 40]); xtickangle(30);
%     changePosition(gca,[0.08 0.01 0.0 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. Responsive Trials (%)',[-0.7 +25 0]});
    changePosition(gca,[0.06 0.01 0.03 -0.55]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. Responsive Trials (%)',[-1.25 25 0]});
    ha = gca; ptable = extras.pvalsTable;
    display_p_table_img(ha,hbs,[0 0.2 0 0.4],ptable);
%     changePosition(gca,[0.06 0.01 0.03 -0.55]); put_axes_labels(gca,{[],[0 0 0]},{'Silent Cells (%)',[-1.25 -10 0]});
%     ha = gca; ptable = extras.pvalsTable;
%     display_p_table_img(ha,hbs,[0 0.23 0 0.37],ptable);%ytickangle(10);
    
%     [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
%     hf = get_figure(5,[8 7 1.75 1.5]);
%     % s = generate_shades(length(bins)-1);
%     tcolors = colors;
%     [hbs,maxY] = plot_bars_sig_lines(mVar,semVar,combs,[p],'colors',tcolors,'sigColor','k',...
%         'ySpacing',15,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
%         'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
%     ylims = ylim;
%     format_axes(gca);
%     set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 260]); format_axes(gca);
%     xticks = xdata; xticklabels = rasterNamesTxt(si);
%     set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 50 100]); xtickangle(45)
%     changePosition(gca,[0.03 0.01 0.05 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. Responsive Trials (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('average_responsive_trials_%d.pdf',ntrials),600);
    %%
%     ff = makeFigureRowsCols(2020,[10 4 6.99 3],'RowsCols',[length(mVar) length(mVar)],'spaceRowsCols',[0 0],'rightUpShifts',[0 0],'widthHeightAdjustment',[0 0]);
    hf = get_figure(2020,[10 3 5.5 3.5]);%axis off;
    [xs,ys] = display_p_table_independent(gca,extras.pvalsTableTxt,[3.5 2.5 0.5 0.7],'');
    xticklabels = rasterNamesTxt(si);
    set(gca,'xtick',xs,'xticklabels',xticklabels,'ytick',ys,'yticklabels',xticklabels); xtickangle(45);ytickangle(45);
    format_axes(gca);
    changePosition(gca,[-0.05 0.01 0.1 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{[],[0 0 0]});
    format_axes(gca);
    save_pdf(hf,mData.pdf_folder,sprintf('average_responsive_trials_%d_p_value_table.pdf',ntrials),600);
    %%
    break;
end

%% compare percent responsive cells
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    si = [Ar_t_D,ArC_t_D,ArCB_t_D,ArB_t_D,Ar_i_T,ArC_i_T,ArCB_i_T,ArB_i_T];
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    good_FR = props1.good_FR;
    good_FR_any = cell_list_op(good_FR,[],'or');
    good_FR_all = cell_list_op(good_FR,[],'and');
    per_active =[]; per_active = [];
    for rr = 1:size(good_FR,1)
        for cc = 1:size(good_FR,2)
            tts = good_FR{rr,cc};
            per_active(rr,cc) = 100*sum(tts)/length(tts);
        end
        tts = good_FR_any{rr,1};        per_active_any(rr) = 100*sum(tts)/length(tts);
        tts = good_FR_all{rr,1};        per_active_all(rr) = 100*sum(tts)/length(tts);
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},[length(si)]);
    dataT = make_between_table({per_active},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = make_xdata([4 4],[1 1.5]);
    h(h==1) = 0;
    hf = get_figure(5,[8 7 2.25 2]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors;
    [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.5);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rasterNamesTxt(si);
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 50 100]); xtickangle(45);
%     changePosition(gca,[0.08 0.01 0.0 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{'Responsive Cells (%)',[-0.7 +15 0]});
    changePosition(gca,[0.06 0.01 0.02 -0.55]); put_axes_labels(gca,{[],[0 0 0]},{'Responsive Cells (%)',[-1.25 -5 0]});
    ha = gca; ptable = extras.pvalsTable;
    ha = display_p_table_img(ha,hbs,[0 0.23 0 0.37],ptable); %ytickangle(20)
%     htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    
%     [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
% %     h(h==1) = 0;
%     hf = get_figure(5,[8 7 2 1.5]);
%     % s = generate_shades(length(bins)-1);
%     tcolors = colors;
%     [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
%         'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
%         'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
%     maxY = maxY + 5;
%     ylims = ylim;
%     format_axes(gca);
%     htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
% %     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
%     set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
%     xticks = xdata; xticklabels = rasterNamesTxt(si);
%     set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 50 100]); xtickangle(45)
%     changePosition(gca,[0.01 0.01 0.05 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Responsive Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end

%% compare percent responsive cells (considering both distance and time rasters of trials and intertrials)
while 1
    ntrials = 50; 
    si = [Lb_T ArL_L_T Lbs_T];          props_light = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ab_T Abs_T];                  props_air_rest = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ar_t_D Ar_t_T];               props_Ar_trials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [ArL_t_D ArL_t_T];             props_ArL_trials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ars_t_D Ars_t_T];             props_Ars_trials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ar_i_D Ar_i_T];               props_Ar_itrials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [ArL_i_D ArL_i_T];             props_ArL_itrials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ars_i_D Ars_i_T];             props_Ars_itrials = get_props_Rs(o.Rs(:,si),ntrials);
    
    gFR_Lb = props_light.good_FR; gFR_Ab = props_air_rest.good_FR;
    gFR_tAr = cell_list_op(props_Ar_trials.good_FR,[],'or');
    gFR_tArL = cell_list_op(props_ArL_trials.good_FR,[],'or');
    gFR_tArs = cell_list_op(props_Ars_trials.good_FR,[],'or');
    gFR_iAr = cell_list_op(props_Ar_itrials.good_FR,[],'or');
    gFR_iArL = cell_list_op(props_ArL_itrials.good_FR,[],'or');
    gFR_iArs = cell_list_op(props_Ars_itrials.good_FR,[],'or');
    
    gFR_Ar = cell_list_op(gFR_tAr,gFR_iAr,'or');
    gFR_ArL = cell_list_op(gFR_tArL,gFR_iArL,'or'); gFR_ArL = cell_list_op(gFR_ArL(:,1),gFR_Lb(:,2),'or');
    gFR_Ars = cell_list_op(gFR_tArs,gFR_iArs,'or');
    
    good_FR = [];
    good_FR = [gFR_Lb(:,[1 3]) gFR_Ab gFR_Ar(:,1) gFR_ArL(:,1) gFR_Ars(:,1)];
    good_FR_any = cell_list_op(good_FR,[],'or');
    good_FR_all = cell_list_op(good_FR,[],'and');
    per_active= [];
    for rr = 1:size(good_FR,1)
        for cc = 1:size(good_FR,2)
            tts = good_FR{rr,cc};
            per_active(rr,cc) = 100*sum(tts)/length(tts);
        end
        tts = good_FR_any{rr,1};        per_active_any(rr) = 100*sum(tts)/length(tts);
        tts = good_FR_all{rr,1};        per_active_all(rr) = 100*sum(tts)/length(tts);
    end
    per_active1 = per_active(:,[1 3 5 6 7 2 4]);
    [within,dvn,xlabels] = make_within_table({'Cond'},[7]);
    dataT = make_between_table({per_active1},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = [1 2 [4 5 6]-0.5 [8 9]-1];
    h(h==1) = 0;
    hf = get_figure(5,[8 7 1.55 2]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.5);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Lb','Ab','Ar','ArL','Ar*','Lb*','Ab*'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 50 100]); xtickangle(45);
%     changePosition(gca,[0.08 0.01 0.0 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{'Responsive Cells (%)',[-0.7 +15 0]});
    changePosition(gca,[0.05 0.01 0.03 -0.55]); put_axes_labels(gca,{[],[0 0 0]},{'Responsive Cells (%)',[-0.25 15 0]});
    ha = gca; ptable = extras.pvalsTable;
    ha = display_p_table_img(ha,hbs,[0 0.25 0 0.35],ptable); %ytickangle(20)
    htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    
%     [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
% %     h(h==1) = 0;
%     hf = get_figure(5,[8 7 2 1.5]);
%     % s = generate_shades(length(bins)-1);
%     tcolors = colors;
%     [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
%         'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
%         'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
%     maxY = maxY + 5;
%     ylims = ylim;
%     format_axes(gca);
%     htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
% %     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
%     set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
%     xticks = xdata; xticklabels = rasterNamesTxt(si);
%     set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 50 100]); xtickangle(45)
%     changePosition(gca,[0.01 0.01 0.05 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Responsive Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d.pdf',ntrials),600);
    %%
    break;
end

%% compare percent silent cells (combined)
while 1
    ntrials = 50; 
    si = [Lb_T ArL_L_T Lbs_T];      props_light = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ab_T Abs_T];              props_air_rest = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ar_t_D Ar_t_T];           props_Ar_trials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [ArL_t_D ArL_t_T];         props_ArL_trials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ars_t_D Ars_t_T];         props_Ars_trials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ar_i_D Ar_i_T];           props_Ar_itrials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [ArL_i_D ArL_i_T];         props_ArL_itrials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ars_i_D Ars_i_T];         props_Ars_itrials = get_props_Rs(o.Rs(:,si),ntrials);
    
        
    gFR_Lb = props_light.silent_cells; gFR_Ab = props_air_rest.silent_cells;
    gFR_tAr = cell_list_op(props_Ar_trials.silent_cells,[],'and');
    gFR_tArL = cell_list_op(props_ArL_trials.silent_cells,[],'and');
    gFR_tArs = cell_list_op(props_Ars_trials.silent_cells,[],'and');
    gFR_iAr = cell_list_op(props_Ar_itrials.silent_cells,[],'and');
    gFR_iArL = cell_list_op(props_ArL_itrials.silent_cells,[],'and');
    gFR_iArs = cell_list_op(props_Ars_itrials.silent_cells,[],'and');
    
    gFR_Ar = cell_list_op(gFR_tAr,gFR_iAr,'and');
    gFR_ArL = cell_list_op(gFR_tArL,gFR_iArL,'and'); gFR_ArL = cell_list_op(gFR_ArL(:,1),gFR_Lb(:,2),'and');
    gFR_Ars = cell_list_op(gFR_tArs,gFR_iArs,'and');
    
    
    good_FR = [gFR_Lb(:,[1 3]) gFR_Ab gFR_Ar(:,1) gFR_ArL(:,1) gFR_Ars(:,1)];
    per_silent = [];
    for rr = 1:size(good_FR,1)
        for cc = 1:size(good_FR,2)
            tts = good_FR{rr,cc};
            per_silent(rr,cc) = 100*sum(tts)/length(tts);
        end
    end
    per_silent1 = per_silent(:,[1 3 5 6 7 2 4]);
    [within,dvn,xlabels] = make_within_table({'Cond'},[7]);
    dataT = make_between_table({per_silent1},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = [1 2 [4 5 6]-0.5 [8 9]-1];
    h(h==1) = 0;
    hf = get_figure(5,[8 7 1.5 1.5]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.5);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Lb','Ab','Ar','ArL','Ar*','Lb*','Ab*'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 50 100]); xtickangle(45);
%     changePosition(gca,[0.08 0.01 0.0 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{'Responsive Cells (%)',[-0.7 +15 0]});
    changePosition(gca,[0.06 0.01 0.03 -0.55]); put_axes_labels(gca,{[],[0 0 0]},{'Silent Cells (%)',[-0.25 15 0]});
    ha = gca; ptable = extras.pvalsTable;
    ha = display_p_table_img(ha,hbs,[0 0.25 0 0.35],ptable); %ytickangle(20)
    save_pdf(hf,mData.pdf_folder,sprintf('silent_cells_across_conditions_%d.pdf',ntrials),600);
    %%
    break;
end

%% compare percent average percentage of responsive trials cells (combined)
while 1
    ntrials = 50; 
    si = [Lb_T ArL_L_T Lbs_T];                          props_light = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ab_t_T Abs_t_T Ab_i_T Abs_i_T];                       props_air_rest = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ar_t_D Ar_t_T];              props_Ar_trials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [ArL_t_D ArL_t_T];             props_ArL_trials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ars_t_D Ars_t_T];             props_Ars_trials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ar_i_D Ar_i_T];             props_Ar_itrials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [ArL_i_D ArL_i_T];            props_ArL_itrials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ars_i_D Ars_i_T];            props_Ars_itrials = get_props_Rs(o.Rs(:,si),ntrials);
    
    gFR_Lb = cell2mat(props_light.N_Resp_Trials); gFR_Ab = cell2mat(props_air_rest.N_Resp_Trials);
    gFR_tAr = mean(cell2mat(props_Ar_trials.N_Resp_Trials),2);
    gFR_tArL = mean(cell2mat(props_ArL_trials.N_Resp_Trials),2);
    gFR_tArs = mean(cell2mat(props_Ars_trials.N_Resp_Trials),2);
    gFR_iAr = mean(cell2mat(props_Ar_itrials.N_Resp_Trials),2);
    gFR_iArL = mean(cell2mat(props_ArL_itrials.N_Resp_Trials),2);
    gFR_iArs = mean(cell2mat(props_Ars_itrials.N_Resp_Trials),2);
    
    gFR_Ar = mean([gFR_tAr,gFR_iAr],2);
    gFR_ArL = mean([gFR_tArL,gFR_iArL],2); gFR_ArL = mean([gFR_ArL(:,1) gFR_Lb(:,2)],2);
    gFR_Ars = mean([gFR_tArs,gFR_iArs],2);
    
    
    good_FR = [gFR_Lb(:,[1 3]) gFR_Ab gFR_Ar gFR_ArL gFR_Ars];
    good_FR1 = good_FR(:,[1 3 5 6 7 2 4]);
    [within,dvn,xlabels] = make_within_table({'Cond'},[7]);
    dataT = make_between_table({good_FR1},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = [1 2 [4 5 6]-0.5 [8 9]-1];
    h(h==1) = 0;
    hf = get_figure(5,[8 7 1.5 1.5]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.5);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Lb','Ab','Ar','ArL','Ar*','Lb*','Ab*'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 50 100]); xtickangle(45);
%     changePosition(gca,[0.08 0.01 0.0 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{'Responsive Cells (%)',[-0.7 +15 0]});
    changePosition(gca,[0.06 0.01 0.03 -0.55]); put_axes_labels(gca,{[],[0 0 0]},{'Silent Cells (%)',[-0.25 15 0]});
    ha = gca; ptable = extras.pvalsTable;
    ha = display_p_table_img(ha,hbs,[0 0.25 0 0.35],ptable); %ytickangle(20)
    save_pdf(hf,mData.pdf_folder,sprintf('silent_cells_across_conditions_%d.pdf',ntrials),600);
    %%
    break;
end

%% compare percent silent cells
while 1
    props1 = get_props_Rs(o.Rs,50); si = si_seqG;
    good_FR = props1.N_Resp_Trials(:,si);
    per_silent =[]; per_active = [];
    for rr = 1:size(good_FR,1)
        for cc = 1:size(good_FR,2)
            tts = good_FR{rr,cc};
            per_silent(rr,cc) = 100*sum(tts == 0)/length(tts);
        end
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},[size(good_FR,2)]);
    dataT = make_between_table({per_silent},dvn);
    ra = RMA(dataT,within,{'lsd','hsd'});
    ra.ranova
%     ra.mauchly
%%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = xdataG;
    h(h==1) = 0;
    hf = get_figure(5,[8 7 2.2 1.5]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.5);
    
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rNaG;
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 50 100]); xtickangle(30);
    changePosition(gca,[0.06 0.01 0.03 -0.55]); put_axes_labels(gca,{[],[0 0 0]},{'Silent Cells (%)',[-1.25 -10 0]});
    ha = gca; ptable = extras.pvalsTable;
    display_p_table_img(ha,hbs,[0 0.23 0 0.37],ptable);%ytickangle(10);
    save_pdf(hf,mData.pdf_folder,sprintf('silent_cells_across_conditions.pdf'),600);
    %%
    break;
end


%% compare percent silent cells si_no_brake_all (both distance and time rasters included for trials and inter-trials)
while 1
    props1 = get_props_Rs(o.Rs,50); si = si_no_brake_allG;
    good_FR = props1.N_Resp_Trials(:,si);
    per_silent =[]; per_active = [];
    for rr = 1:size(good_FR,1)
        for cc = 1:size(good_FR,2)
            tts = good_FR{rr,cc};
            per_silent(rr,cc) = 100*sum(tts == 0)/length(tts);
        end
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},[size(good_FR,2)]);
    dataT = make_between_table({per_silent},dvn);
    ra = RMA(dataT,within,{'lsd','hsd'});
    ra.ranova
%     ra.mauchly
%%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
%     xdata = xdataG;
    h(h==1) = 0;
    hf = get_figure(5,[8 7 2.2 1.5]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.5);
    
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rasterNamesTxt(si);
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 50 100]); xtickangle(30);
    changePosition(gca,[0.06 0.01 0.03 -0.55]); put_axes_labels(gca,{[],[0 0 0]},{'Silent Cells (%)',[-1.25 -10 0]});
    ha = gca; ptable = extras.pvalsTable;
    display_p_table_img(ha,hbs,[0 0.23 0 0.37],ptable);%ytickangle(10);
    save_pdf(hf,mData.pdf_folder,sprintf('silent_cells_across_conditions_no_brake.pdf'),600);
    %%
    break;
end



%% compare percent unique cells
while 1
    ntrials = 50; 
    si = [Lb_T ArL_L_T Lbs_T];          props_light = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ab_T Abs_T];                  props_air_rest = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ar_t_D Ar_t_T];               props_Ar_trials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [ArL_t_D ArL_t_T];             props_ArL_trials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ars_t_D Ars_t_T];             props_Ars_trials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ar_i_D Ar_i_T];               props_Ar_itrials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [ArL_i_D ArL_i_T];             props_ArL_itrials = get_props_Rs(o.Rs(:,si),ntrials);
    si = [Ars_i_D Ars_i_T];             props_Ars_itrials = get_props_Rs(o.Rs(:,si),ntrials);
    
    gFR_Lb = props_light.good_FR; gFR_Ab = props_air_rest.good_FR;
    gFR_tAr = cell_list_op(props_Ar_trials.good_FR,[],'or');
    gFR_tArL = cell_list_op(props_ArL_trials.good_FR,[],'or');
    gFR_tArs = cell_list_op(props_Ars_trials.good_FR,[],'or');
    gFR_iAr = cell_list_op(props_Ar_itrials.good_FR,[],'or');
    gFR_iArL = cell_list_op(props_ArL_itrials.good_FR,[],'or');
    gFR_iArs = cell_list_op(props_Ars_itrials.good_FR,[],'or');
    
    gFR_Ar = cell_list_op(gFR_tAr,gFR_iAr,'or');
    gFR_ArL = cell_list_op(gFR_tArL,gFR_iArL,'or'); gFR_ArL = cell_list_op(gFR_ArL(:,1),gFR_Lb(:,2),'or');
    gFR_Ars = cell_list_op(gFR_tArs,gFR_iArs,'or');
    
    good_FR = [];
    good_FR = [gFR_Lb(:,[1 3]) gFR_Ab gFR_Ar(:,1) gFR_ArL(:,1) gFR_Ars(:,1)];
    per_unique =[];
    for rr = 1:size(good_FR,2)
        for cc = 1:size(good_FR,2)
            if rr == cc
                condMat(cc) = cc;
            else
                condMat(cc) = -cc;
            end
        end
        temp_unique = get_cell_list(good_FR,condMat,1);
        per_unique(:,rr) = temp_unique(:,1)*100;
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},[7]);
    dataT = make_between_table({per_unique},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = [1 2 [4 5 6]-0.5 [8 9]-1];
    h(h==1) = 0;
    hf = get_figure(5,[8 7 1.5 1.5]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.5);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Lb','Ab','Ar','ArL','Ar*','Lb*','Ab*'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 5 100]); xtickangle(45);
%     changePosition(gca,[0.08 0.01 0.0 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{'Responsive Cells (%)',[-0.7 +15 0]});
    changePosition(gca,[0.06 0.01 0.03 -0.55]); put_axes_labels(gca,{[],[0 0 0]},{'Unique Responsive Cells (%)',[-0.55 5 0]});
    ha = gca; ptable = extras.pvalsTable;
    ha = display_p_table_img(ha,hbs,[0 0.25 0 0.35],ptable); %ytickangle(20)
    save_pdf(hf,mData.pdf_folder,sprintf('unique_cells_across_conditions.pdf'),600);
    %%
    break;
end


%% compare percent unique cells in sets of conditions
while 1
     ntrials = 50; si = si_seq;
    props1 = get_props_Rs(o.Rs,ntrials);
    good_FR = props1.good_FR(:,si);
    per_unique =[];
    condMat = [-1 -2 3 -4 -5 -6 -7 -8 -9 -10 -11;...
                -1 -2 -3 -4 5 -6 -7 -8 -9 -10 -11;...
                -1 -2 -3 -4 -5 -6 7 -8 -9 -10 -11];
    per_unique = get_cell_list(good_FR,condMat,1)*100;
    [within,dvn,xlabels] = make_within_table({'Cond'},[length(si)]);
    dataT = make_between_table({per_unique},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    ra.mauchly
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    h(h==1) = 0;
    hf = get_figure(5,[8 7 1.7 1.5]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.5);
    
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rasterNamesTxt(si);
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    changePosition(gca,[0.08 0.01 0.0 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{'Unique Responsive Cells (%)',[-1 3 0]});
    ha = gca; ptable = extras.pvalsTable;
    display_p_table(ha,hbs,[0 -0.07 0 0.9],ptable);
    save_pdf(hf,mData.pdf_folder,sprintf('unique_cells_across_conditions.pdf'),600);
    %%
    break;
end

%% compare diference number of responsive trials across conditions
while 1
    props1 = get_props_Rs(o.Rs,50); si = si_seq(1:10);
    good_FR = props1.N_Resp_Trials(:,si);
    for rr = 1:size(good_FR,1)
        for cc = 2:size(good_FR,2)
            ttsp = good_FR{rr,cc-1}; tts = good_FR{rr,cc};
            diff_resp_trials{rr,cc-1} = tts-ttsp;
        end
    end
    var1 = exec_fun_on_cell_mat(diff_resp_trials,'mean');
    [within,dvn,xlabels] = make_within_table({'Cond'},[size(good_FR,2)-1]);
    dataT = make_between_table({var1},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    ra.mauchly
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.25 1 1]);
%     h(h==1) = 0;
    hf = get_figure(5,[8 7 1.75 1.5]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.17);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rasterNamesTxt(si);
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 100]); xtickangle(45)
    changePosition(gca,[0.035 0.01 0.05 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Silent Cells (%)',[0 0 0]});
    pos = get(gca,'Position');
    
    save_pdf(hf,mData.pdf_folder,sprintf('diff_responsive_trials.pdf'),600);
    %%
    break;
end


%% view all population vectors
while 1 
    props1 = get_props_Rs(o.Rs,5); 
    si = si_air_time_trials;
    si = si_seq;
    resp = props1.good_FR;%(:,si);
    all_conds = 1:size(resp,2);
    resp_o = get_cell_list(resp,[si]');
    resp_o1 = get_cell_list(resp,[-setdiff(all_conds,si)]);
    resp_ = cell_list_op(resp_o,resp_o1,'and');
    resp = props1.good_FR(:,si);
    resp_FR_or = cell_list_op(resp,[],'or'); resp_FR_and = cell_list_op(resp,[],'and');
    per_resp_FR_or = exec_fun_on_cell_mat(resp_FR_or,{'sum','length'});
    view_population_vector(o.Rs(:,si),o.mR(:,si),resp,100);
break;
end

%% Overlap Indices ImageSC Single animal
while 1
    ntrials = 50;
    si = [Ar_t_D,ArC_t_D,ArCB_t_D,ArB_t_D,Ar_i_T,ArC_i_T,ArCB_i_T,ArB_i_T];
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    good_FR = props1.good_FR;
    si = [Lb_T ArL_L_T Lbs_T Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
    props1 = get_props_Rs(o.Rs,ntrials);
    resp = [props1.good_FR(:,si)];% resp_speed];
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    mOI = OI_mat(:,:,4);
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = 0.6;%max([mOI(:);semOI(:)]);    
    minI = 0;%min([mOI(:);semOI(:)])
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = rasterNamesTxt(si); 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));   % imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.6 1.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 11.5],[0.5 11.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    text(-0.6,12.1,'Overlap Index (Representative Animal)','FontSize',5);
    set(gca,'Ydir','normal');ytickangle(20);
    box on;
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d.pdf',ntrials),600);
    break;
end

%% Overlap Indices ImageSC
while 1
    ntrials = 50;
    si = [Ar_t_D,ArC_t_D,ArCB_t_D,ArB_t_D,Ar_i_T,ArC_i_T,ArCB_i_T,ArB_i_T];
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    resp = props1.good_FR;

    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = rasterNamesTxt(si); 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.5 1.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    text(0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',ntrials),600);
    %%
    %
%     maxYss = [maxYs(9) maxYs(2) maxYs(7) maxYs(8) maxYs(7) maxYs(8) maxYs(7) maxYs(8) maxYs(9) maxYs(2) maxYs(11)];
    for sel_row = 1:sz
        sel_row
        for rr = 1:size(OI_mat,1)
            for cc = 1:size(OI_mat,2)
                if isnan(OI_mat(rr,cc,1))
                    continue;
                end
                if h_vals(rr,cc) == 1
        %             text(cc,rr,'*','Color','r','FontSize',12);
                end
            end
        end
    %     sel_row = 4;
        OIs = squeeze(OI_mat(sel_row,:,:))'; OIs(:,sel_row) = [];
        [within,dvn,xlabels1] = make_within_table({'Cond'},[size(OI_mat,1)-1]);
        dataT = make_between_table({OIs},dvn);
        ra = RMA(dataT,within);
        ra.ranova
        [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
        xdata = 1:size(OI_mat,1); 
        xdata(sel_row) = [];
%         xdata = xdataG;
        h(h==1) = 0;
        hf = get_figure(5,[8 7 1.95 1.5]);
        tcolors = mData.colors(setdiff(1:size(OI_mat,1),sel_row));
        [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
            'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.45);
        ylims = ylim;
        format_axes(gca);
        set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
        xticks = xdata; xticklabels = txl;
        set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
        changePosition(gca,[0.09 0.01 -0.03 -0.55]); put_axes_labels(gca,{[],[0 0 0]},{sprintf('Overlap Index of %s',txl{sel_row}),[-1.45 0.1 0]});
        ha = gca; ptable = extras.pvalsTable;
        display_p_table_img(ha,hbs,[0 0.23 0 0.37],ptable);ytickangle(0);
        save_pdf(hf,mData.pdf_folder,sprintf('OI_bar_%d.pdf',sel_row),600);
        maxYs(sel_row) = maxY;
    end
    %%
    break;
end

%% agglomerative hierarchical clustering
while 1
    mOI1 = mOI;
    mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1);
    tree = linkage(Di);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','right','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 1.25 2]);
    set(H,'linewidth',1);
    set(gca,'yticklabels',txl(TC));ytickangle(30);
    format_axes(gca);
    hx = xlabel('Eucledian Distance');%changePosition(hx,[-0.051 0 0]);
    changePosition(gca,[0 0.0 0.05 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);
    %%
    break;
end

%% Overlap Indices ImageSC for no brake conditions
while 1
    ntrials = 50;
    si = si_no_brake_allG;
    props1 = get_props_Rs(o.Rs,ntrials);
    resp = [props1.good_FR(:,si)];% resp_speed];
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)])
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = rasterNamesTxt(si); 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.6 1.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    text(0.5,13.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean_no_brake.pdf',ntrials),600);
    %%
    %
%     maxYss = [maxYs(9) maxYs(2) maxYs(7) maxYs(8) maxYs(7) maxYs(8) maxYs(7) maxYs(8) maxYs(9) maxYs(2) maxYs(11)];
    for sel_row = 1:11
        sel_row
        for rr = 1:size(OI_mat,1)
            for cc = 1:size(OI_mat,2)
                if isnan(OI_mat(rr,cc,1))
                    continue;
                end
                if h_vals(rr,cc) == 1
        %             text(cc,rr,'*','Color','r','FontSize',12);
                end
            end
        end
    %     sel_row = 4;
        OIs = squeeze(OI_mat(sel_row,:,:))'; OIs(:,sel_row) = [];
        [within,dvn,xlabels1] = make_within_table({'Cond'},[size(OI_mat,1)-1]);
        dataT = make_between_table({OIs},dvn);
        ra = RMA(dataT,within);
        ra.ranova
%         [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
%         xdata = 1:11; xdata(sel_row) = [];
% 
%     %     h(h==1) = 0;
%         hf = get_figure(5,[8 7 1.95 1.5]);
%         % s = generate_shades(length(bins)-1);
%         tcolors = mData.colors(setdiff(1:11,sel_row));
%         [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
%             'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
%             'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.01,'capsize',3);
%         ylims = ylim;
% %         maxY = maxYss(sel_row);
%         format_axes(gca); set(gca,'FontSize',16,'ytick',[0 0.2 0.4 0.6 0.8 1]);
%         set_axes_limits(gca,[0.35 11+.65],[ylims(1) maxY]); format_axes(gca);
%         xticks = xdata; xticklabels = txl;
%         set(gca,'xtick',1:11,'xticklabels',xticklabels); xtickangle(45)
%         changePosition(gca,[0.025 -0.01 0.05 0.05]); 
        [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
        xdata = xdataG;%1:size(OI_mat,1); 
        xdata(sel_row) = [];
%         xdata = xdataG;
        h(h==1) = 0;
        hf = get_figure(5,[8 7 1.95 1.5]);
        % s = generate_shades(length(bins)-1);
        tcolors = mData.colors(setdiff(1:size(OI_mat,1),sel_row));
        [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
            'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.45);

        ylims = ylim;
        format_axes(gca);
%         if ii < 11
        set_axes_limits(gca,[0.35 xdataG(end)+.65],[ylims(1) maxY]); format_axes(gca);
%         else
%         set_axes_limits(gca,[0.35 11.65],[ylims(1) maxY]); format_axes(gca);
%         end
        xticks = xdataG; xticklabels = txl;
        set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
        changePosition(gca,[0.09 0.01 -0.03 -0.55]); put_axes_labels(gca,{[],[0 0 0]},{sprintf('Overlap Index of %s',txl{sel_row}),[-1.45 0.1 0]});
        ha = gca; ptable = extras.pvalsTable;
        display_p_table_img(ha,hbs,[0 0.23 0 0.37],ptable);ytickangle(0);
        save_pdf(hf,mData.pdf_folder,sprintf('OI_bar_no_brake_%d.pdf',sel_row),600);
        maxYs(sel_row) = maxY;
    end
    %%
    break;
end
