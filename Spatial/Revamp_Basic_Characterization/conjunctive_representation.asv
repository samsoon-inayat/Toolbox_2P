function conjunctive_representation

%% Load Data
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei = evalin('base','ei'); 
    
    selContexts = [1 4 6 2 7 3 3 4 4 5 5 3 3 4 4 5 5 0 0 2 2 7 7];
    rasterNames = {'light22T','light22T','light22T','air55T','air55T','airD','airID','airD','airID','airD','airID','airT','airIT','airT','airIT','airT','airIT','motionOnsets','motionOffsets','airRT','airIRT','airRT','airIRT'};
    rasterNamesTxt = {'C1-T','C4L-T','C1''-T','C2-T','C2''-T','C3t-D','C3it-D','C4t-D','C4it-D','C3''t-D','C3''it-D','C3t-T','C3it-T','C4t-T','C4it-T','C3''t-T','C3''it-T','MOn-T','MOff-T','C2t-T','C2it-T','C2''t-T','C2''it-T'};
    o = get_data(ei,selContexts,rasterNames);
    
%     selContexts = [3 4 5];
%     rasterNames = {'air77T','air77T','air77T'};
%     opct = get_data(ei,selContexts,rasterNames);
    
    for ii = 1:length(selContexts)
        all_xl{ii} = sprintf('%c%c-%d',rasterNames{ii}(1),rasterNames{ii}(end),selContexts(ii));
    end
    
    [speedRs,resp_speed] = get_speed_response(ei);
    all_xl{ii+1} = 'sp';
    resp = [o.resp.vals resp_speed];
%     resp_o = cell_list_op(resp,[],'xor')
    si_light = [1 2 3];
    si_air_rest = [4 5];
    si_air_dist_trials = [6 8 10];
    si_air_dist_itrials = [7 9 11];
    si_air_time_trials = [12 14 16];
    si_air_time_itrials = [13 15 17];
    si_motion = [18 19];
    si_seq_m = [1 4 6 13 8 15 10 17 3 5 2 18 19];
    si_seq = [1 4 6 13 8 15 10 17 3 5 2];
    si_seqT = [1 4 12 13 14 15 16 17 3 5 2 18 19];
    si_seq_f = [1 20 21 6 13 8 15 10 17 3 22 23 2];
    si_no_brake = [6 13 8 15 10 17];
    si_no_brake_dist = [6 7 8 9 10 11];
    si_no_brake_time = [12 13 14 15 16 17];
    
    dzMI = prop_op(o.props.zMI(:,[si_air_dist_trials si_air_dist_itrials]),o.props.zMI(:,[si_air_time_trials si_air_time_itrials]),0.1);
    
    resp = [resp dzMI.resp_D_g_T(:,1) dzMI.resp_T_g_D(:,1)];
    resp_OR = cell_list_op(resp,[],'or');
    resp_AND = cell_list_op(resp,[],'and');
    all_xl = [all_xl {'D_g_T','T_g_D'}];
    
    break
end
n = 0;
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
    props1 = get_props_Rs(o.Rs,ntrials);
    si = si_seq;%si_no_brake_dist;
    good_FR = props1.good_FR(:,si);
    all_zMIs = props1.rs(:,si);
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
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
%     h(h==1) = 0;
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.561,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rasterNames(si);
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[-0.01 -0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. zMI',[0 0 0]});
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
    ntrials = 50; si = si_seq;
    props1 = get_props_Rs(o.Rs,ntrials);
    good_FR = props1.N_Resp_Trials(:,si);
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
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 20 40]); xtickangle(45);
    changePosition(gca,[0.08 0.01 0.0 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. Responsive Trials (%)',[-0.7 +25 0]});
    ha = gca; ptable = extras.pvalsTable;
    display_p_table(ha,hbs,[0 -0.07 0 0.9],ptable);
    
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
    ntrials = 50; si = si_seq;
    props1 = get_props_Rs(o.Rs,ntrials);
    good_FR = props1.good_FR(:,si);
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
    ra.mauchly
    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
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
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 50 100]); xtickangle(45);
    changePosition(gca,[0.08 0.01 0.0 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{'Responsive Cells (%)',[-0.7 +15 0]});
    ha = gca; ptable = extras.pvalsTable;
    ha = display_p_table(ha,hbs,[0 -0.07 0 0.9],ptable);
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

%% compare percent silent cells
while 1
    props1 = get_props_Rs(o.Rs,50); si = si_seq;
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
    ra.mauchly
    %%
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
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 50 100]); xtickangle(45);
    changePosition(gca,[0.08 0.01 0.0 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{'Silent Cells (%)',[-0.5 0 0]});
    ha = gca; ptable = extras.pvalsTable;
    display_p_table(ha,hbs,[0 -0.07 0 0.9],ptable);
    save_pdf(hf,mData.pdf_folder,sprintf('silent_cells_across_conditions.pdf'),600);
    %%
    break;
end


%% compare percent unique cells
while 1
     ntrials = 50; si = si_seq;
    props1 = get_props_Rs(o.Rs,ntrials);
    good_FR = props1.good_FR(:,si);
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
    [within,dvn,xlabels] = make_within_table({'Cond'},[length(si)]);
    dataT = make_between_table({per_unique},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    ra.mauchly
    %%
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
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 3 100]); xtickangle(45);
    changePosition(gca,[0.08 0.01 0.0 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{'Unique Responsive Cells (%)',[-1 3 0]});
    ha = gca; ptable = extras.pvalsTable;
    display_p_table(ha,hbs,[0 -0.07 0 0.9],ptable);
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
    props1 = get_props_Rs(o.Rs,50); si = si_seq;
    good_FR = props1.N_Resp_Trials(:,si);
    per_silent =[]; per_active = [];
    for rr = 1:size(good_FR,1)
        for cc = 1:size(good_FR,2)
            tts = good_FR{rr,cc};
            per_silent(rr,cc) = 100*sum(tts == 10)/length(tts);
        end
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},[size(good_FR,2)]);
    dataT = make_between_table({per_silent},dvn);
    ra = RMA(dataT,within,{'lsd','hsd'});
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
    
    save_pdf(hf,mData.pdf_folder,sprintf('silent_cells_across_conditions.pdf'),600);
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
    si = si_seq;
    props1 = get_props_Rs(o.Rs,ntrials);
    resp = [props1.good_FR(:,si)];% resp_speed];
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    mOI = OI_mat(:,:,4);
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = 0;%min([mOI(:);semOI(:)])
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = rasterNamesTxt(si); 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    imAlpha(isnan(mask))=0.25; imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.6 1.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 11.5],[0.5 11.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    text(0.5,-0.1,'Overlap Index (Representative Animal)','FontSize',5);
    box on;
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d.pdf',ntrials),600);
    break;
end

%% Overlap Indices ImageSC
while 1
    ntrials = 50;
    si = si_seq;
    props1 = get_props_Rs(o.Rs,ntrials);
    resp = [props1.good_FR(:,si)];% resp_speed];
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = 0;%min([mOI(:);semOI(:)])
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = rasterNamesTxt(si); 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    imAlpha(isnan(mask))=0.25; imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.6 1.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 11.5],[0.5 11.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',[],'Ydir','reverse'); xtickangle(45);
    text(0.5,-0.1,'Average Overlap Index (N = 5 animals)','FontSize',5);
    box on
    hc = putColorBar(gca,[0.08 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',ntrials),600);
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
        [within,dvn,xlabels1] = make_within_table({'Cond'},[11-1]);
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
        xdata = 1:11; xdata(sel_row) = [];
        h(h==1) = 0;
        hf = get_figure(5,[8 7 1.95 1.5]);
        % s = generate_shades(length(bins)-1);
        tcolors = mData.colors(setdiff(1:11,sel_row));
        [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
            'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.5);

        ylims = ylim;
        format_axes(gca);
        set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
        xticks = xdata; xticklabels = txl;
        set(gca,'xtick',1:11,'xticklabels',xticklabels); xtickangle(45);
        changePosition(gca,[0.03 0.01 0.045 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{sprintf('Overlap Index of %s',txl{sel_row}),[-0.1 0.15 0]});
        ha = gca; ptable = extras.pvalsTable;
        display_p_table(ha,hbs,[0 -0.07 0 0.9],ptable);
        save_pdf(hf,mData.pdf_folder,sprintf('OI_bar_%d.pdf',sel_row),600);
        maxYs(sel_row) = maxY;
    end
    %%
    break;
end
