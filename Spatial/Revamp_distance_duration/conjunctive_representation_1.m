function conjunctive_representation

%% Load Data
tic
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei = evalin('base','ei'); 
    
    selContexts = [1 4 6 2 7 3 3 4 4 5 5 3 3 4 4 5 5 0 0 2 2 7 7 2 7 2 7 3 4 5 3 4 5 3 4 5 3 4 5];
    rasterNames = {'light22T','light22T','light22T','air55T','air55T','airD','airID','airD','airID','airD','airID','airT','airIT','airT','airIT','airT','airIT','motionOnsets','motionOffsets','airRT','airIRT','airRT','airIRT',...
        'airOnsets22T','airOnsets22T','airOffsets22T',...
                    'airOffsets22T','airOnsets22T','airOnsets22T','airOnsets22T',...
                    'airOffsets22T','airOffsets22T','airOffsets22T','beltD','beltD','beltD','beltT','beltT','beltT'};
    rasterNamesTxt = {'Lb-T','ArL-L-T','Lb*-T','Ab-T','Ab*-T','Ar-t-D','Ar-i-D','ArL-t-D','ArL-i-D','Ar*-t-D','Ar*-i-D','Ar-t-T','Ar-i-T','ArL-t-T','ArL-i-T','Ar*-t-T','Ar*-i-T','MOn-T','MOff-T','Ab-t-T','Ab-i-T','Ab*-t-T','Ab*-i-T',...
        '2-AOn','7-AOn','2-AOff','7-AOff',...
        '3-AOn','4-AOn','5-AOn','3-AOff','4-AOff','5-AOff','3-B-D','4-B-D','5-B-D','3-B-T','4-B-T','5-B-T'};
    xlabelsSeq = {'Lb','ArL-L','Lb*','Ab','Ab*','Ar-t','ArL-t','Ar*-t','Ar-i','ArL-i','Ar*-i'};

    o = get_data(ei,selContexts,rasterNames);
    for ii = 1:length(selContexts)
        all_xl{ii} = sprintf('%c%c-%d',rasterNames{ii}(1),rasterNames{ii}(end),selContexts(ii));
    end
    Lb_T = 1; ArL_L_T = 2; Lbs_T = 3; Ab_T = 4; Abs_T = 5; Ab_t_T = 20; Abs_t_T = 22;  Ab_i_T = 21; Abs_i_T = 23;
    Ar_t_D = 6; Ar_t_T = 12; ArL_t_D = 8; ArL_t_T = 14; Ars_t_D = 10; Ars_t_T = 16;
    Ar_i_D = 7; Ar_i_T = 13; ArL_i_D = 9; ArL_i_T = 15; Ars_i_D = 11; Ars_i_T = 17;
    MOn_T = 18; MOff_T = 19;
    Ab_On = 24; Abs_On = 25; Ab_Off = 26; Abs_Off = 27; 
    Ar_On = 28; ArL_On = 29; Ars_On = 30; Ar_Off = 31; ArL_Off = 32; Ars_Off = 33; Ar_B_D = 34; ArL_B_D = 35; Ars_B_D = 36;
    Ar_B_T = 37; ArL_B_T = 38; Ars_B_T = 39;
    
%     [speedRs,resp_speed,speed_percent,resp_speedAcc] = load_speed_response(ei);
%     all_xl{ii+1} = 'sp';
%     resp = [o.resp.vals resp_speed];
  
    M = [18 19];
   
%     dzMI = prop_op(o.props.zMI(:,[Ar_t_D Ar_i_D]),o.props.zMI(:,[Ar_t_T Ar_i_T]),0.1);
    break
end
toc
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
%%
an = 4; cn = 3;
% respC = cell_list_op(FD_conj{2},dzMI_FD.resp_complex,'and'); %all_responsive_cells{an,cn}
respC = FD_Dis_comp{2}; respC = FD_conj{2};
figure(1000);clf;subplot 141;imagesc(RsTt{an,cn}.speed); set(gca,'Ydir','normal'); subplot 142;imagesc(RsDt{an,cn}.speed);set(gca,'Ydir','normal'); subplot 143;imagesc(RsTi{an,cn}.speed); set(gca,'Ydir','normal'); subplot 144;imagesc(RsDi{an,cn}.speed);set(gca,'Ydir','normal');
plotRasters_dis_dur({RsTt{an,cn},RsDt{an,cn},RsTi{an,cn},RsDi{an,cn}},find(respC{an,cn}));
%% Show sample rasters dist
while 1
   Rs = RsDt; RsT = RsTt;
%    Rs = RsDi;
   an = 4; cn = 3;
    % plotRasters_dis_dur({Rs{an,cn},RsT{an,cn}},find(FT_conj{an,cn}))
    % find(resp_valsC{an}(:,cn));
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 5],...
        'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
        [-75 -475]);
    set(gcf,'color','w'); set(gcf,'Position',[10 4 4 1]);
    ff = sample_rasters(Rs{an,cn},[112,72,92,6,123],ff);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersD'),600);
    break;
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
        [-75 -475]);
    set(gcf,'color','w'); set(gcf,'Position',[10 4 3.25 1]);
    ff = sample_rasters(Rs{an,cn},[328 518 567 436],ff);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersDT'),600);
    break;
end
%% Show sample rasters dur
while 1
   Rs = RsTi;
   an = 3; cn = 3;
   props1 = get_props_Rs(Rs,50);
    % plotRasters_simplest(Rs{an,cn},find(props1.vals{an,cn}))
    % find(resp_valsC{an}(:,cn));
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 5],...
        'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
        [-75 -475]);
    set(gcf,'color','w'); set(gcf,'Position',[10 4 4 1]);
    ff = sample_rasters(Rs{an,cn},[64,117,311,136,53],ff);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersT'),600);
    break;
end

%% compare the zMIs
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    si = [Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    lnsi = length(si);
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    good_FR = [propsD.vals propsTt.vals propsDi.vals propsT.vals];
    good_FR = [propsD.vals propsD.vals propsDi.vals propsDi.vals];
    all_zMIs = props1.zMI;
    zMIs = [];
    for rr = 1:size(all_zMIs,1)
        for cc = 1:size(all_zMIs,2)
            resp = good_FR{rr,cc};
            tzmis = all_zMIs{rr,cc};
            zMIs(rr,cc) = nanmean(tzmis(resp));
        end
    end
    [within,dvn,xlabels] = make_within_table({'TI','DT','Cond'},[2,2,3]);
    dataT = make_between_table({zMIs},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova
    
    %%
     ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    si = [Ar_t_D ArL_t_D Ars_t_D Ar_i_D ArL_i_D Ars_i_D Ar_t_D ArL_t_D Ars_t_D Ar_i_D ArL_i_D Ars_i_D];
    lnsi = length(si);
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    good_FR = [propsD.vals propsD.vals propsDi.vals propsDi.vals];
    all_zMIs = props1.zMI;
    zMIs = [];
    for rr = 1:size(all_zMIs,1)
        for cc = 1:size(all_zMIs,2)
            resp = good_FR{rr,cc};
            tzmis = all_zMIs{rr,cc};
            zMIs(rr,cc) = nanmean(tzmis(resp));
        end
    end
    [within,dvn,xlabels] = make_within_table({'CT','TI','Cond'},[2,2,3]);
    dataT = make_between_table({zMIs},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova
    
     %%
     ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    si = [Ar_t_T ArL_t_T Ars_t_T Ar_i_T ArL_i_T Ars_i_T Ar_t_T ArL_t_T Ars_t_T Ar_i_T ArL_i_T Ars_i_T];
    lnsi = length(si);
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    good_FR = [propsTt.vals propsTt.vals propsT.vals propsT.vals];
    all_zMIs = props1.zMI;
    zMIs = [];
    for rr = 1:size(all_zMIs,1)
        for cc = 1:size(all_zMIs,2)
            resp = good_FR{rr,cc};
            tzmis = all_zMIs{rr,cc};
            zMIs(rr,cc) = nanmean(tzmis(resp));
        end
    end
    [within,dvn,xlabels] = make_within_table({'CT','TI','Cond'},[2,2,3]);
    dataT = make_between_table({zMIs},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova
    
    
%%    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT_by_TI','hsd'},[1 1 1]);
    p = ones(size(p)); p(1) = ra.MC.hsd.CT_by_TI.pValue(1); p(end) = ra.MC.hsd.CT_by_TI.pValue(end); h = p<0.05;
    xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 1.75 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'yspacing',0.5);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Trials','Intertrials'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 round(maxY/2) round(maxY)]); xtickangle(45);
    changePosition(gca,[0.01 -0.01 0.06 0.01]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. zMI',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('zMIs_good_FR_all_Conditions.pdf'),600);
%%    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'DT','hsd'},[1 1 1]);
    p = ones(size(p)); p(1) = ra.MC.hsd.DT_by_TI.pValue(1); p(end) = ra.MC.hsd.DT_by_TI.pValue(end); h = p<0.05;
    xdata = make_xdata([2],[1 1.5]);
    hf = get_figure(5,[8 7 1.75 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'yspacing',0.5);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dis','Dur'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 round(maxY/2) round(maxY)]); xtickangle(45);
    changePosition(gca,[0.01 -0.01 0.06 0.01]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. zMI',[0 0 0]});
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
    xticks = xdata; xticklabels = rasterNamesTxt(si);
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
    si = [Ar_t_D ArL_t_D Ars_t_D Ar_i_D ArL_i_D Ars_i_D Ar_t_T ArL_t_T Ars_t_T Ar_i_T ArL_i_T Ars_i_T];
    lnsi = length(si);
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    good_FR = props1.vals;
    good_FR_any = cell_list_op(good_FR,[],'or');
    good_FR_all = cell_list_op(good_FR,[],'and');
    per_active =[]; 
    for rr = 1:size(good_FR,1)
        for cc = 1:size(good_FR,2)
            tts = good_FR{rr,cc};
            per_active(rr,cc) = 100*sum(tts)/length(tts);
        end
        tts = good_FR_any{rr,1};        per_active_any(rr) = 100*sum(tts)/length(tts);
        tts = good_FR_all{rr,1};        per_active_all(rr) = 100*sum(tts)/length(tts);
    end
    [within,dvn,xlabels] = make_within_table({'DT','TI','Cond'},[2,2,3]);
    dataT = make_between_table({per_active},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova
    %%
    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
 %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_DT','hsd'},[1.5 1 1]);
    p = ones(size(p)); 
    ind = ismember(combs,[1 2],'rows');  p(ind) = ra.MC.hsd.TI_by_DT.pValue(1); 
    ind = ismember(combs,[3 4],'rows'); p(ind) = ra.MC.hsd.TI_by_DT.pValue(end); h = p<0.05;
    
    xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%)',round(mra),round(semra)),'FontSize',6);
%     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'T','I'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.01 0.01 -0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
    xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dis','Dur'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.01 0.01 -0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end


%% compare percent trials fidelity
while 1
     ntrials = 50; 
    si = [Ar_t_D ArL_t_D Ars_t_D Ar_i_D ArL_i_D Ars_i_D Ar_t_T ArL_t_T Ars_t_T Ar_i_T ArL_i_T Ars_i_T];
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    gFR = props1.vals; rf = props1.N_Resp_Trials;
    all_gFR = exec_fun_on_cell_mat(rf,'nanmean',gFR);

    [within,dvn,xlabels] = make_within_table({'DT','TI','Cond'},[2,2,3]);
    dataT = make_between_table({all_gFR},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova
    %%
    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
 %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_DT','hsd'},[1.5 1 1]);
    p = ones(size(p)); 
    ind = ismember(combs,[1 2],'rows');  p(ind) = ra.MC.hsd.TI_by_DT.pValue(1); 
    ind = ismember(combs,[3 4],'rows'); p(ind) = ra.MC.hsd.TI_by_DT.pValue(end); h = p<0.05;
    
    xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%)',round(mra),round(semra)),'FontSize',6);
%     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'T','I'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.01 0.01 -0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
    xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dis','Dur'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.01 0.01 -0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end


%% compare percent responsive cells
while 1
    ntrials = 50; 
    si = [Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
    props1 = get_props_Rs(o.Rs,ntrials);
    good_FR = props1.vals(:,si);
%     good_FR  = [respDT.inh respDT.exc propsD.good_FR propsT.good_FR]
%     good_FR = [respDT.exc respDT.inh];
    good_FR_any = cell_list_op(good_FR,[],'or');
    good_FR_all = cell_list_op(good_FR,[],'and');
    per_active =[]; 
    for rr = 1:size(good_FR,1)
        for cc = 1:size(good_FR,2)
            tts = good_FR{rr,cc};
            per_active(rr,cc) = 100*sum(tts)/length(tts);
        end
        tts = good_FR_any{rr,1};        per_active_any(rr) = 100*sum(tts)/length(tts);
        tts = good_FR_all{rr,1};        per_active_all(rr) = 100*sum(tts)/length(tts);
    end
    [within,dvn,xlabels] = make_within_table({'DT','Cond'},[2,3]);
    dataT = make_between_table({per_active},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'DT_by_Cond','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
    htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%)',round(mra),round(semra)),'FontSize',6);
%     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.01 0.01 -0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
    xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
%     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 50 100]); xtickangle(45)
    changePosition(gca,[0.01 0.01 -0.3 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end

%% compare percent responsive cells
while 1
    ntrials = 50; 
    si = [Ar_t_D ArL_t_D Ars_t_D];
    si = [Ar_i_T ArL_i_T Ars_i_T];
    props1 = get_props_Rs(o.Rs,ntrials);
    good_FR = props1.vals(:,si);
%     good_FR  = [respDT.inh respDT.exc propsD.good_FR propsT.good_FR]
%     good_FR = [respDT.exc respDT.inh];
    good_FR_any = cell_list_op(good_FR,[],'or');
    good_FR_all = cell_list_op(good_FR,[],'and');
    per_active =[]; 
    for rr = 1:size(good_FR,1)
        for cc = 1:size(good_FR,2)
            tts = good_FR{rr,cc};
            per_active(rr,cc) = 100*sum(tts)/length(tts);
        end
        tts = good_FR_any{rr,1};        per_active_any(rr) = 100*sum(tts)/length(tts);
        tts = good_FR_all{rr,1};        per_active_all(rr) = 100*sum(tts)/length(tts);
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},[3]);
    dataT = make_between_table({per_active},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 10;
    ylims = ylim;
    format_axes(gca);
    htxt = text(0.75,maxY-2,{'Any Condition',sprintf('(%d\x00B1%d%%)',round(mra),round(semra))},'FontSize',6);
%     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 25 50]); xtickangle(45)
    changePosition(gca,[0.04 0.01 -0.4 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end

%% compare percent of conjunctive and complementary cells
while 1
    ntrials = 50; 
    si = [Ar_t_D ArL_t_D Ars_t_D];
%     si = [Ar_i_T ArL_i_T Ars_i_T];
    props1 = get_props_Rs(o.Rs,ntrials);
    good_FR = props1.vals(:,si);
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(good_FR,0.5,0.05);
    mOI = mCI; semOI = semCI;
    
    clear respRV conjV comp1V comp2V
    for an = 1:5
        respRV(an,:) = diag(all_CI_mat(:,:,an));
        conjV(an,:) = diag(all_CI_mat(:,:,an),1);
        comp1V(an,:) = diag(uni(:,:,an),1);
        comp2V(an,:) = diag(uni(:,:,an),-1);
    end
    
    for an = 1:5
        conjV(an,3) = all_CI_mat(1,3,an);
        comp1V(an,3) = uni(1,3,an);
        comp2V(an,3) = uni(3,1,an);
    end
    
    %%
    [within,dvn,xlabels] = make_within_table({'Type','Cond'},[3,3]);
    dataT = make_between_table({conjV,comp1V,comp2V},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
    ra.ranova

    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Type_by_Cond','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3],[1 2]);
    hf = get_figure(5,[8 7 2.75 1]);
%     hf = get_figure(5,[8 7 2.5 1]);
    tcolors = repmat(mData.dcolors(1:3),1,3);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    maxY1 = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3-4','4-5','3-5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 5 10 15]); xtickangle(45)
    changePosition(gca,[-0.03 0.01 -0.0 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Type','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
    xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(4:6);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY1]); format_axes(gca);
    xticks = xdata; xticklabels = {'3-4','4-5','3-5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 5 10]); xtickangle(45);
    changePosition(gca,[0.04 0.01 -0.5 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end

%% compare percent of fidelity (Trials %)
while 1
    ntrials = 50; 
    si = [Ar_t_D ArL_t_D Ars_t_D];
%     si = [Ar_i_T ArL_i_T Ars_i_T];
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    gFR = props1.vals; rf = props1.N_Resp_Trials;
    all_gFR = exec_fun_on_cell_mat(rf,'nanmean',gFR);

    [within,dvn,xlabels] = make_within_table({'Cond'},[3]);
    dataT = make_between_table({all_gFR},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
    ra.ranova

    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
    xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 25 50 75]); xtickangle(45)
    changePosition(gca,[0.04 0.01 -0.4 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Trials (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end

%% compare zMI
while 1
    ntrials = 50; 
    si = [Ar_t_D ArL_t_D Ars_t_D];
%     si = [Ar_i_T ArL_i_T Ars_i_T];
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    gFR = props1.vals; rf = props1.zMI;
    all_gFR = exec_fun_on_cell_mat(rf,'nanmean',gFR);

    [within,dvn,xlabels] = make_within_table({'Cond'},[3]);
    dataT = make_between_table({all_gFR},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
    ra.ranova

    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
    xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY ;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 1 5]); xtickangle(45)
    changePosition(gca,[0.15 0.01 -0.4 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Z-Score'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end

%% compare PWs
while 1
    ntrials = 50; 
    si = [Ar_t_D ArL_t_D Ars_t_D];
%     si = [Ar_i_T ArL_i_T Ars_i_T];
    props1 = get_props_Rs(o.Rs(:,si),ntrials,1);
    gFR = props1.vals; rf = props1.PWs;
    all_gFR = exec_fun_on_cell_mat(rf,'nanmean',gFR);

    [within,dvn,xlabels] = make_within_table({'Cond'},[3]);
    dataT = make_between_table({all_gFR},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
    ra.ranova

    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY ;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
%     set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 5 10 20]); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.4 -0.05]); 
    put_axes_labels(gca,{[],[0 0 0]},{{'Width (cm)'},[0 0 0]});
%     put_axes_labels(gca,{[],[0 0 0]},{{'Width (sec)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end

%% Overlap Indices ImageSC
while 1
%     ntrials = 50; 
% %     si = [Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T];
% % %     si = [Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
%     si = [Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
% %     lnsi = length(si);
%     props1 = get_props_Rs(o.Rs(:,si),ntrials);
%     resp = props1.vals;
%     
%     si = [Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
% %     si = [Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
%     siG = si; RsG = o.Rs(:,si); propsG = get_props_Rs(RsG,ntrials); respG = propsG.vals;
%     
%     [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(respG,0.5,0.05);
    mOI = mCI; semOI = semCI;
%     mOI = mean(uni,3); semOI = std(uni,[],3)/sqrt(5);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:)]);%;semOI(:)]);    
    minI = min([mOI(:)]);%semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = {'3-T-Dis','4-T-Dis','5-T-Dis','3-T-Dur','4-T-Dur','5-T-Dur','3-T-Mix','4-T-Mix','5-T-Mix',...
        '3-I-Dis','4-I-Dis','5-I-Dis','3-I-Dur','4-I-Dur','5-I-Dur','3-I-Mix','4-I-Mix','5-I-Mix'};
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 3.5 3.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    for rr = 1:size(mCI,1)
        for cc = 1:size(mCI,1)
            if rr == cc
                text(rr-0.25,cc,sprintf('%.0f',mCI(rr,cc)),'FontSize',5);
            end
        end
    end
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
%     set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
%     text(-0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); 
    set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[-0.01 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.09 0.05]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',ntrials),600);
    %%
    break;
end
%% agglomerative hierarchical clustering
while 1
    %%
    mOI1 = mOI;
%     mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1,@naneucdist);
%     tree = linkage(mOI1,'average','euclidean');
    tree = linkage(Di,'average');
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 3.5 1.5]);
    set(H,'linewidth',1);
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel('Eucledian Distance');changePosition(hx,[0 -0.1 0]);
    changePosition(gca,[0.03 0.0 0.05 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);
    %%
    break;
end
