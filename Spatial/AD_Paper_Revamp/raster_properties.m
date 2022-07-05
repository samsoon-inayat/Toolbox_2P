function raster_properties

%% Three-Way RM-ANOVA Most Recent compare all variables pooled air on or off brake vs no-brake exc and inh (and controls) big ANOVA test
while 1
    prop_names = {'resp','N_Resp_Trials','zMI','zMINaN','HaFD','HiFD','cells_pooled'};
    event_type = {'1-D','2-D','3-D','4-D','1-T','2-T','3-T','4-T'};
    sic = {[C1_t_D];[C2_t_D];[C3_t_D];[C4_t_D];[C1_i_T];[C2_i_T];[C3_i_T];[C4_i_T]};
%     event_type = {'1-D','2-D','3-D','4-D'};
%     sic = {[C1_t_D];[C2_t_D];[C3_t_D];[C4_t_D]};
    pni = 3;
    [all_gFR_C,all_gV_C] = return_values_props(oC,sic,pni);
    [all_gFR_A,all_gV_A] = return_values_props(oA,sic,pni);
    %%
    %%
    varC = all_gV_C;
    varA = all_gV_A;
    [within,dvn,xlabels] = make_within_table({'TI','Cond'},[2,4]);
    dataT = make_between_table({varC;varA},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
    ra.ranova
    %%
    varC = all_gV_C;
    varA = all_gV_A;
    [within,dvn,xlabels] = make_within_table({'Cond'},[4]);
    dataT = make_between_table({varC;varA},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
    ra.ranova
    %%
    break;
end

%% compare the zMIs
while 1
    ntrials = 50; 
    si = [C1_t_D C2_t_D C3_t_D C4_t_D C1_i_T C2_i_T C3_i_T C4_i_T];
    o = oC; Rs = o.Rs(:,si);mR = o.mR(:,si); props1 = get_props_Rs(Rs,ntrials); good_FR_C = props1.zMI;
    o = oA; Rs = o.Rs(:,si);mR = o.mR(:,si); props1 = get_props_Rs(Rs,ntrials); good_FR_A = props1.zMI;
    
%     varC = find_percent(good_FR_C);
%     varA = find_percent(good_FR_A);
    
    varC = exec_fun_on_cell_mat(good_FR_C,'nanmean');
    varA = exec_fun_on_cell_mat(good_FR_A,'nanmean');
    
    [within,dvn,xlabels] = make_within_table({'TI','Cond'},[2,4]);
    dataT = make_between_table({varC;varA},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
    ra.ranova
    %%
    break;
end





%% compare the zMIs
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    si = [C1_t_D C2_t_D C3_t_D C4_t_D C1_i_T C2_i_T C3_i_T C4_i_T];
    lnsi = length(si);
    props1C = get_props_Rs(oC.Rs(:,si),ntrials); all_zMIsC = props1C.zMI; zMIsC = []; 
    props1A = get_props_Rs(oA.Rs(:,si),ntrials); all_zMIsA = props1A.zMI; zMIsA = []; 
    respC = props1C.good_FR_and_tuned;%cell_list_op(props1C.vals,props1C.good_FR,'and');
    respA = props1A.good_FR_and_tuned;%cell_list_op(props1A.vals,props1A.good_FR,'and');
    for rr = 1:size(all_zMIsA,1)
        for cc = 1:size(all_zMIsA,2)
            resp = respC{rr,cc}; tzmis = all_zMIsC{rr,cc}; zMIsC(rr,cc) = nanmean(tzmis(resp));
            resp = respA{rr,cc}; tzmis = all_zMIsA{rr,cc}; zMIsA(rr,cc) = nanmean(tzmis(resp));
        end
    end
    %%
    [within,dvn,xlabels] = make_within_table({'TI','Cond'},[2,4]);
    dataT = make_between_table({zMIsC;zMIsA},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
    ra.ranova
    
%%    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_TI_Cond','hsd'},[1 1 1]);
%     p = ones(size(p)); p(1) = ra.MC.hsd.CT_by_TI.pValue(1); p(end) = ra.MC.hsd.CT_by_TI.pValue(end); h = p<0.05;
    xdata = make_xdata([8 8],[1 1.5]);
    hf = get_figure(5,[8 7 1.75 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.colors,2,1);
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


