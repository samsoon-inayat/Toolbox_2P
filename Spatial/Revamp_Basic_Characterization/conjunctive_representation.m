function conjunctive_representation

%% Load Data
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei = evalin('base','ei'); 
    
    selContexts = [1 4 6 2 7 3 3 4 4 5 5 3 3 4 4 5 5 0 0];
    rasterNames = {'light22T','light22T','light22T','air55T','air55T','airD','airID','airD','airID','airD','airID','airT','airIT','airT','airIT','airT','airIT','motionOnsets','motionOffsets'};
    o = get_data(ei,selContexts,rasterNames);
    
%     selContexts = [3 4 5];
%     rasterNames = {'air77T','air77T','air77T'};
%     opct = get_data(ei,selContexts,rasterNames);
    
    for ii = 1:length(selContexts)
        xlabels{ii} = sprintf('%c%c-%d',rasterNames{ii}(1),rasterNames{ii}(end),selContexts(ii));
    end
    
    [speedRs,resp_speed] = get_speed_response(ei);
    xlabels{ii+1} = 'sp';
    resp = [o.resp.vals resp_speed];
%     resp_o = cell_list_op(resp,[],'xor')
    si_light = [1 2 3];
    si_air_rest = [4 5];
    si_air_dist_trials = [6 8 10];
    si_air_dist_itrials = [7 9 11];
    si_air_time_trials = [12 14 16];
    si_air_time_itrials = [13 15 17];
    si_motion = [18 19];
    si_seq = [1 4 6 13 8 15 10 17 3 5 18 19];
    si_no_brake = [6 13 8 15 10 17];
    si_no_brake_dist = [6 7 8 9 10 11];
    si_no_brake_time = [12 13 14 15 16 17];
    
    dzMI = prop_op(o.props.zMI(:,[si_air_dist_trials si_air_dist_itrials]),o.props.zMI(:,[si_air_time_trials si_air_time_itrials]),0.1);
    
    resp = [resp dzMI.resp_D_g_T(:,1) dzMI.resp_T_g_D(:,1)];
    resp_OR = cell_list_op(resp,[],'or');
    resp_AND = cell_list_op(resp,[],'and');
    xlabels = [xlabels {'D_g_T','T_g_D'}];
    
    break
end
n = 0;
%% compare the MIs
while 1
    props1 = get_props_Rs(o.Rs,3);
    si = si_no_brake;
    good_FR = props1.good_FR(:,si);
    all_zMIs = props1.MI(:,si);
    zMIs = [];
    for rr = 1:size(all_zMIs,1)
        for cc = 1:size(all_zMIs,2)
            resp = good_FR{rr,cc};
            tzmis = all_zMIs{rr,cc};
            zMIs(rr,cc) = nanmean(tzmis(resp));
        end
    end
    [within,dvn,xlabels] = make_within_table({'Cond','DT'},[3,2]);
    dataT = make_between_table({zMIs},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    ra.mauchly
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'DT','bonferroni'},[1 1 1]);
%     h(h==1) = 0;
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
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

%% compare percent active cells
while 1
    all_trial_scores = o.props.trial_scores(:,si_seq);
    all_zMIs = o.props.zMI(:,si_seq);
    all_rs = o.props.rs(:,si_seq);
    per_silent =[]; per_active = []; cttszmi = [];
    for rr = 1:size(all_trial_scores,1)
        tts = [];
        zMIs = [];
        rs = [];
        for cc = 1:size(all_trial_scores,2)
            tts(:,cc) = all_trial_scores{rr,cc};
            zMIs(:,cc) = all_zMIs{rr,cc};
            rs(:,cc) = all_rs{rr,cc};
        end
%         tts1 = tts(:,1:9); tts2 = tts(:,2:10);
%         pctts = 100*(tts2 - tts1)./tts1;
%         dtts = diff(tts,[],2);
%         mdtts(rr,:) = mean(dtts);
        per_silent(rr,:) = 100*sum(tts == 0)/size(tts,1);
        per_active(rr,:) = 100*sum(tts > 0.5)/size(tts,1);
        zMIs(isnan(zMIs)) = 0;
        corr_tts_zMI = corr(tts,zMIs);
        mask = triu(ones(size(corr_tts_zMI)),0) & ~triu(ones(size(corr_tts_zMI)),1);
        cttszmi(rr,:) = corr_tts_zMI(mask);
        all_zMIs_sel{rr} = (tts > 0.5) .* zMIs;
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},[12]);
    dataT = make_between_table({per_active},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    ra.mauchly
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','bonferroni'},[1 1 1]);
    h(h==1) = 0;
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rasterNames(si_seq);
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[-0.01 -0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Active Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('change_in_trial_scores_across_conditions.pdf'),600);
    %%
    break;
end

%% compare percent silent cells
while 1
    all_trial_scores = o.props.trial_scores(:,si_seq);
    all_zMIs = o.props.zMI(:,si_seq);
    all_rs = o.props.rs(:,si_seq);
    per_silent =[]; per_active = []; cttszmi = [];
    for rr = 1:size(all_trial_scores,1)
        tts = [];
        zMIs = [];
        rs = [];
        for cc = 1:size(all_trial_scores,2)
            tts(:,cc) = all_trial_scores{rr,cc};
            zMIs(:,cc) = all_zMIs{rr,cc};
            rs(:,cc) = all_rs{rr,cc};
        end
%         tts1 = tts(:,1:9); tts2 = tts(:,2:10);
%         pctts = 100*(tts2 - tts1)./tts1;
%         dtts = diff(tts,[],2);
%         mdtts(rr,:) = mean(dtts);
        per_silent(rr,:) = 100*sum(tts == 0)/size(tts,1);
        per_active(rr,:) = 100*sum(tts > 0.5)/size(tts,1);
        zMIs(isnan(zMIs)) = 0;
        corr_tts_zMI = corr(tts,zMIs);
        mask = triu(ones(size(corr_tts_zMI)),0) & ~triu(ones(size(corr_tts_zMI)),1);
        cttszmi(rr,:) = corr_tts_zMI(mask);
        all_zMIs_sel{rr} = (tts > 0.5) .* zMIs;
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},[12]);
    dataT = make_between_table({per_silent},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    ra.mauchly
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','bonferroni'},[1 1 1]);
    h(h==1) = 0;
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rasterNames(si_seq);
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[-0.01 -0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Silent Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('silent_cells_across_conditions.pdf'),600);
    %%
    break;
end

%% compare percentage of cells active in atleast 1 condition
while 1
    all_trial_scores = o.props.trial_scores(:,si_seq);
    all_zMIs = o.props.zMI(:,si_seq);
    all_rs = o.props.rs(:,si_seq);
    per_silent =[]; per_active = []; cttszmi = [];
    for rr = 1:size(all_trial_scores,1)
        tts = [];
        zMIs = [];
        rs = [];
        for cc = 1:size(all_trial_scores,2)
            tts(:,cc) = all_trial_scores{rr,cc};
            zMIs(:,cc) = all_zMIs{rr,cc};
            rs(:,cc) = all_rs{rr,cc};
        end
%         tts1 = tts(:,1:9); tts2 = tts(:,2:10);
%         pctts = 100*(tts2 - tts1)./tts1;
%         dtts = diff(tts,[],2);
%         mdtts(rr,:) = mean(dtts);
        per_silent(rr,:) = 100*sum(tts == 0)/size(tts,1);
        per_active(rr,:) = 100*sum(tts > 0.5)/size(tts,1);
        zMIs(isnan(zMIs)) = 0;
        corr_tts_zMI = corr(tts,zMIs);
        mask = triu(ones(size(corr_tts_zMI)),0) & ~triu(ones(size(corr_tts_zMI)),1);
        cttszmi(rr,:) = corr_tts_zMI(mask);
        all_zMIs_sel{rr} = (tts > 0.5) .* zMIs;
        pocaaoc{rr} = sum(tts > 0.5,2);
        pocaaoc1(rr) = sum(sum(tts >= 0.5,2)>=1)/size(tts,1);
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},[12]);
    dataT = make_between_table({per_silent},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    ra.mauchly
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','bonferroni'},[1 1 1]);
    h(h==1) = 0;
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rasterNames(si_seq);
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[-0.01 -0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Silent Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('silent_cells_across_conditions.pdf'),600);
    %%
    break;
end
%% make distribution of cells responding in different environments

while 1
    distD = pocaaoc';
    [distDo,allVals] = getAveragesAndAllValues(distD);
    minBin = min(allVals);
    maxBin = max(allVals);
    incr = 1; %maxBin =
    hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    %    [ha,hb,hca] = plotDistributions(distD,'colors',tcolors,'maxY',maxBin,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    [ha,hb,hca] = plotAverageDistributions(distD,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'pdf_or_cdf','pdf','graphType','bar');
    set_axes_limits(gca,[minBin-0.75 maxBin+0.75],[]);
    break
end

%%
props1 = get_props_Rs(o.Rs,4);
si = si_air_time_trials;
si = si_seq;
all_conds = 1:size(resp,2);
resp_o = get_cell_list(resp,[si]');
resp_o1 = get_cell_list(resp,[-setdiff(all_conds,si)]);
resp_ = cell_list_op(resp_o,resp_o1,'and');
resp = props1.good_FR(:,si);
resp_FR_or = cell_list_op(resp,[],'or'); resp_FR_and = cell_list_op(resp,[],'and');
per_resp_FR_or = exec_fun_on_cell_mat(resp_FR_or,{'sum','length'});
view_population_vector(o.Rs(:,si),o.mR(:,si),resp,100);
%%
si = si_seq;
props1 = get_props_Rs(o.Rs,3);
resp = props1.good_FR(:,si);
[OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
maxOI = max([mOI(:);semOI(:)]); 
minOI = min([mOI(:);semOI(:)]);
maxOI = max(OI_mat,[],3);
txl = xlabels(si);
figure(100);clf;
% subplot 121
fOI = maxOI;
imagesc(fOI,[min(fOI(:)) max(fOI(:))]);
grid on;
set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);colorbar;
axis equal
for rr = 1:size(OI_mat,1)
    for cc = 1:size(OI_mat,2)
        if isnan(OI_mat(rr,cc,1))
            continue;
        end
        if h_vals(rr,cc) == 1
            text(cc,rr,'*','Color','r','FontSize',12);
        end
    end
end
OIs = squeeze(OI_mat(3,:,:))';
    [within,dvn,xlabels1] = make_within_table({'Cond'},[5]);
    dataT = make_between_table({OIs(:,4:8)},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    ra.mauchly
    %%
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','bonferroni'},[1 1 1]);
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
    xticks = xdata; xticklabels = txl(4:end);
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[-0.01 -0.02 -0.05 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Active Cells (%)',[0 0 0]});
