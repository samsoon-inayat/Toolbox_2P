function dist_dur


%% general for all properties including responsivity, response fidelity, zMI, Rs
while 1
    ntrials = 50; 
    si = [Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
%     si = [C1_i_T C2_i_T C3_i_T C4_i_T];
    Rs_C = o.Rs(:,si);mRs_C = o.mR(:,si);
    props_C = get_props_Rs(Rs_C,ntrials);
    pop_var_name = {'all','vals','valsT','Nvals','good_zMI','Ngood_zMI'};
    pop_var_name = {'vals'};
    sel_pop_C = cell_list_op(props_C,pop_var_name); 
    
    params = {'perc','N_Resp_Trials','zMI','rs','nan_zMI','nan_rs','HaFD','HiFD','PWs','centers','peak_locations'};
    varT = 1;%:length(params)
    for pii = varT
        if pii == 1
            mean_var_C = exec_fun_on_cell_mat(sel_pop_C,'percent'); 
        else
            eval(sprintf('var_C = get_vals(props_C.%s,sel_pop_C);',params{pii})); 
            if pii == 5 || pii == 6
                mean_var_C = exec_fun_on_cell_mat(sel_pop_C,'percent'); 
            else
                mean_var_C = exec_fun_on_cell_mat(var_C,'nanmean'); 
            end
        end
    end
    varC = mean_var_C;
    [within,dvn,xlabels,awithinD] = make_within_table({'TI','DT','Cond'},[2,2,3]);
    dataT = make_between_table({varC},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd'}});
    ra.ranova
break;
end
%% separate populations based on dist and time
while 1
%     sel_pop_C = cell_list_op(props_C,pop_var_name);

    dist_based_trials = sel_pop_C(:,awithinD(:,1) == 1 & awithinD(:,2) == 1);
    dist_based_itrials = sel_pop_C(:,awithinD(:,1) == 2 & awithinD(:,2) == 1);

    time_based_trials = sel_pop_C(:,awithinD(:,1) == 1 & awithinD(:,2) == 2);
    time_based_itrials = sel_pop_C(:,awithinD(:,1) == 2 & awithinD(:,2) == 2);

    dist_trials = cell_list_op(dist_based_trials,time_based_trials,'separate');
    time_trials = cell_list_op(time_based_trials,dist_based_trials,'separate');
    both_trials = cell_list_op(time_based_trials,dist_based_trials,'and');

    dist_itrials = cell_list_op(dist_based_itrials,time_based_itrials,'separate');
    time_itrials = cell_list_op(time_based_itrials,dist_based_itrials,'separate');
    both_itrials = cell_list_op(time_based_itrials,dist_based_itrials,'and');

    all_types_cat = [dist_trials both_trials time_trials dist_itrials both_itrials time_itrials];
    ii = 1; cc_types = {'D','B','T'}; rr_types = {'T','I'}; cn_types = {'1','2','3'};
    for rr = 1:2
        for cc = 1:3
            for cn = 1:3
                cell_types{1,ii} = sprintf('%s-%s-%s',cn_types{cn},cc_types{cc},rr_types{rr}); ii = ii + 1;
            end
        end
    end
    %%
    ntrials = 50; 
    si = [Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    Rs_C = o.Rs(:,si);mRs_C = o.mR(:,si);
    props_C = get_props_Rs(Rs_C,ntrials);
    
    params = {'perc','N_Resp_Trials','zMI','rs','nan_zMI','nan_rs','HaFD','HiFD','PWs','centers','peak_locations'};
    varT = 1;%:length(params)
    for pii = varT
        if pii == 1
            mean_var_C = exec_fun_on_cell_mat(sel_pop_C,'percent'); 
        else
            eval(sprintf('var_C = get_vals(props_C.%s,sel_pop_C);',params{pii})); 
            if pii == 5 || pii == 6
                mean_var_C = exec_fun_on_cell_mat(sel_pop_C,'percent'); 
            else
                mean_var_C = exec_fun_on_cell_mat(var_C,'nanmean'); 
            end
        end
    end
    varC = mean_var_C;
    %%
    varC = find_percent(all_types_cat);
    [within,dvn,xlabels,awithinD1] = make_within_table({'TI','DT','Cond'},[2,3,3]);
    dataT = make_between_table({varC},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd'}});
    ra.ranova
    
    %%
    alpha = 0.05/3;
    for cni = 1:3
        varC1 = varC(:,awithinD1(:,3) == cni);
        [within,dvn,xlabels] = make_within_table({'TI','DT'},[2,3]);
        dataT = make_between_table({varC1},dvn);
        ra_cond{cni} = RMA(dataT,within,{alpha,{'hsd'}});
        ra_cond{cni}.ranova
    end
    %%
    alpha = 0.05/2;
    for ti = 1:2
        varC1 = varC(:,awithinD1(:,1) == ti);
        [within,dvn,xlabels] = make_within_table({'DT','Cond'},[3,3]);
        dataT = make_between_table({varC1},dvn);
        ra_ti{ti} = RMA(dataT,within,{alpha,{'hsd'}});
        ra_ti{ti}.ranova
    end
    %%
    break;
end

%%
RsDt = o.Rs(:,[Ar_t_D ArL_t_D Ars_t_D]);  RsTt = o.Rs(:,[Ar_t_T ArL_t_T Ars_t_T]);
RsDi = o.Rs(:,[Ar_i_D ArL_i_D Ars_i_D]);  RsTi = o.Rs(:,[Ar_i_T ArL_i_T Ars_i_T]);
[dzMI_FD,dzMI_FT] = get_zMI_comp_dist_time(RsDt,RsTt,RsDi,RsTi);

%%
RsDtMC = o.RsMC(:,[Ar_t_D ArL_t_D Ars_t_D]);  RsTtMC = o.RsMC(:,[Ar_t_T ArL_t_T Ars_t_T]);
RsDiMC = o.RsMC(:,[Ar_i_D ArL_i_D Ars_i_D]);  RsTiMC = o.RsMC(:,[Ar_i_T ArL_i_T Ars_i_T]);

%%
while 1
    respfids = {[10 30],[40 100],[0 100],[50 100],[60 100]};
    raster_types = {'RsTt','RsDt','RsTi','RsDi'};
    % raster_types = {'RsTt','RsTi'};
    clear props
    for ii = 1:length(respfids)
        trialsR = respfids{ii};
        for jj = 1:length(raster_types)
            cmdTxt = sprintf('props{ii,jj} = get_props_Rs(%s,trialsR);',raster_types{jj});
            eval(cmdTxt);
        end
    end
    % find dis, dur, and mix cells

    clear FD_Dis_comp FD_Dur_comp FD_conj FT_Dis_comp FT_Dur_comp FT_conj
    for ii = 1:length(respfids)
        FD_Dur = cell_list_op(props{ii,1}.vals,props{ii,1}.good_FR,'and'); FD_Dis = cell_list_op(props{ii,2}.vals,props{ii,2}.good_FR,'and');
        FT_Dur = cell_list_op(props{ii,3}.vals,props{ii,3}.good_FR,'and'); FT_Dis = cell_list_op(props{ii,4}.vals,props{ii,4}.good_FR,'and');
%         FD_Dur = cell_list_op(props{ii,1}.vals,props{ii,1}.good_zMI,'and'); FD_Dis = cell_list_op(props{ii,2}.vals,props{ii,2}.good_zMI,'and');
%         FT_Dur = cell_list_op(props{ii,3}.vals,props{ii,3}.good_zMI,'and'); FT_Dis = cell_list_op(props{ii,4}.vals,props{ii,4}.good_zMI,'and');
%         FD_Dur = props{ii,1}.vals; FD_Dis = props{ii,2}.vals;
%         FT_Dur = props{ii,3}.vals; FT_Dis = props{ii,4}.vals;
        
        cellP1 = FD_Dis; cellP2 = FD_Dur;
        FD_Dis_comp{ii} = cell_list_op(cellP1,cell_list_op(cellP2,[],'not'),'and');
        FD_Dur_comp{ii} = cell_list_op(cell_list_op(cellP1,[],'not'),cellP2,'and');
        FD_conj{ii} = cell_list_op(cellP1,cellP2,'and');
        
        cellP1 = FT_Dis; cellP2 = FT_Dur;
        FT_Dis_comp{ii} = cell_list_op(cellP1,cell_list_op(cellP2,[],'not'),'and');
        FT_Dur_comp{ii} = cell_list_op(cell_list_op(cellP1,[],'not'),cellP2,'and');
        FT_conj{ii} = cell_list_op(cellP1,cellP2,'and');
%         an = 4; cn = 3; respC = FT_Dur_comp{2}; tempCL = respC{an,cn};
    end
break;
end

%%
while 1
    %%
    % for one RF
    cni = 1:3;
    rfi = 1;
    resp = [FD_Dur_comp{rfi}(:,cni) FD_Dis_comp{rfi}(:,cni) FD_conj{rfi}(:,cni) FT_Dur_comp{rfi}(:,cni) FT_Dis_comp{rfi}(:,cni) FT_conj{rfi}(:,cni)];
    rfi = 2;
    resp = [resp FD_Dur_comp{rfi}(:,cni) FD_Dis_comp{rfi}(:,cni) FD_conj{rfi}(:,cni) FT_Dur_comp{rfi}(:,cni) FT_Dis_comp{rfi}(:,cni) FT_conj{rfi}(:,cni)];

    per_resp = 100*exec_fun_on_cell_mat(resp,'sum')./exec_fun_on_cell_mat(resp,'length');

%     [within,dvn,xlabels] = make_within_table({'RF','TI','CT'},[2,2,3]);
%     dataT = make_between_table({per_resp},dvn);
%     ra = RMA(dataT,within,{'hsd'});
%     ra.ranova

    [within,dvn,xlabels] = make_within_table({'RF','TI','CT','Cond'},[2,2,3,3]);
    dataT = make_between_table({per_resp},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova

    any_cells = cell_list_op(resp,[],'or',1);
    per_active_any = 100*exec_fun_on_cell_mat(any_cells,'sum')./exec_fun_on_cell_mat(any_cells,'length');
    any_cells = cell_list_op(resp,[],'and',1);
    per_active_all = 100*exec_fun_on_cell_mat(any_cells,'sum')./exec_fun_on_cell_mat(any_cells,'length');

    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
    %%
    cni = 1:3;
    rfi = 1;
    resp = [FD_Dur_comp{rfi}(:,cni) FD_Dis_comp{rfi}(:,cni) FD_conj{rfi}(:,cni) FT_Dur_comp{rfi}(:,cni) FT_Dis_comp{rfi}(:,cni) FT_conj{rfi}(:,cni)];
    rfi = 2;
    resp = [resp FD_Dur_comp{rfi}(:,cni) FD_Dis_comp{rfi}(:,cni) FD_conj{rfi}(:,cni) FT_Dur_comp{rfi}(:,cni) FT_Dis_comp{rfi}(:,cni) FT_conj{rfi}(:,cni)];

    per_resp = 100*exec_fun_on_cell_mat(resp,'sum')./exec_fun_on_cell_mat(resp,'length');

    [within,dvn,xlabels] = make_within_table({'RF','TI','CT','Cond'},[2,2,3,3]);
    dataT = make_between_table({per_resp},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova

    any_cells = cell_list_op(resp,[],'or',1);
    per_active_any = 100*exec_fun_on_cell_mat(any_cells,'sum')./exec_fun_on_cell_mat(any_cells,'length');
    any_cells = cell_list_op(resp,[],'and',1);
    per_active_all = 100*exec_fun_on_cell_mat(any_cells,'sum')./exec_fun_on_cell_mat(any_cells,'length');

    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'RF','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([2],[1 1.5]);
    hf = get_figure(5,[8 7 3.5 1]);
    tcolors = repmat(mData.colors(1:9),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[-0.04 0.01 0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_all.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_CT','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3],[1 1.5]);
    hf = get_figure(5,[8 7 3.5 1]);
    tcolors = repmat(mData.colors(1:9),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dur','Dis','Mix'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[-0.04 0.01 0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_all.pdf'),600);
    %%
    break
end

%% check the difference in zMI for Dis and Dur for the different cell types DurC, DisC, and DDM
while 1
    %%
    mean_dzMI = [];
    FD_Prop = dzMI_FD.diff_T_D; FT_Prop = dzMI_FT.diff_T_D; 
%     FD_Prop = dzMI_FD.rs.diff_T_D; FT_Prop = dzMI_FT.rs.diff_T_D; 
%     FD_Prop = dzMI_FD.HaFD.diff_T_D; FT_Prop = dzMI_FT.HaFD.diff_T_D; 
%     FD_Prop = dzMI_FD.HiFD.diff_T_D; FT_Prop = dzMI_FT.HiFD.diff_T_D; 
    for rfi = 2
        TD = FD_Prop;
%         cell_resp = FD_Dur_comp{rfi};
%         mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp)];
%         cell_resp = FD_Dis_comp{rfi};
%         mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp)];
        cell_resp = FD_conj{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp)];

        TD = FT_Prop;
        cell_resp = FT_Dur_comp{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp)];
%         cell_resp = FT_Dis_comp{rfi};
%         mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp)];
%         cell_resp = FT_conj{rfi};
%         mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp)];
    end
    
    [within,dvn,xlabels] = make_within_table({'TI','Cond'},[2,3]);
    dataT = make_between_table({mean_dzMI},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova
    %%
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI','hsd'},[1.5 1 1]);
    xdata = make_xdata([2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(7:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.75,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1)-1 maxY+1]); format_axes(gca);
    xticks = xdata; xticklabels = {'T','I'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.4 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('dzMI_CT.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI','hsd'},[1.5 1 1]);
    xdata = make_xdata([2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(10:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.04,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.025);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'T','I'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.5 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]})
    %%
    break;
end


%% Overlap Indices ImageSC
while 1
    %%
%     si = [Ar_On ArL_On Ars_On]; props_ON = get_props_Rs(o.Rs(:,si)); resp_ON = cell_list_op(props_ON.vals,[],'or',1);
%     si = [Ar_Off ArL_Off Ars_Off]; props_OFF = get_props_Rs(o.Rs(:,si)); resp_OFF = cell_list_op(props_OFF.vals,[],'or',1);
    si = [Ab_On Abs_On]; props_ON = get_props_Rs(o.Rs(:,si)); resp_ON = cell_list_op(props_ON.vals,[],'or',1);
    si = [Ab_Off Abs_Off]; props_OFF = get_props_Rs(o.Rs(:,si)); resp_OFF = cell_list_op(props_OFF.vals,[],'or',1);
    rfi = 3;
    respAll = [resp_ON cell_list_op(FD_Dur_comp{rfi},[],'or',1) cell_list_op(FD_Dis_comp{rfi},[],'or',1) cell_list_op(FD_conj{rfi},[],'or',1) ...
        cell_list_op(FT_Dur_comp{rfi},[],'or',1) cell_list_op(FT_Dis_comp{rfi},[],'or',1) cell_list_op(FT_conj{rfi},[],'or',1) resp_OFF];
%     respAll = [resp_ON FD_conj{rfi} FT_Dur_comp{rfi} resp_OFF];
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(respAll,0.5,0.05);
    mOI = mCI; semOI = semCI;
%     mOI = mean(uni,3); semOI = std(uni,[],3)/sqrt(5);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:)]);%;semOI(:)]);    
    minI = min([mOI(:)]);%semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = {'A-On','T-Dur','T-Dis','T-Mix',...
        'I-Dur','I-Dis','I-Mix','A-Off'};
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
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',rfi),600);
    %%
    break;
end

%% Overlap Indices ImageSC
while 1
    %%
%     si = [Ar_On ArL_On Ars_On]; props_ON = get_props_Rs(o.Rs(:,si)); resp_ON = cell_list_op(props_ON.vals,[],'or',1);
%     si = [Ar_Off ArL_Off Ars_Off]; props_OFF = get_props_Rs(o.Rs(:,si)); resp_OFF = cell_list_op(props_OFF.vals,[],'or',1);
    si = [Ab_On Abs_On]; props_ON = get_props_Rs(o.Rs(:,si)); resp_ON = cell_list_op(props_ON.vals,[],'or',1);
    si = [Ab_Off Abs_Off]; props_OFF = get_props_Rs(o.Rs(:,si)); resp_OFF = cell_list_op(props_OFF.vals,[],'or',1);
    rfi = 3;
    respAll = [resp_ON FD_Dur_comp{rfi} FD_Dis_comp{rfi} FD_conj{rfi} FT_Dur_comp{rfi} FT_Dis_comp{rfi} FT_conj{rfi} resp_OFF];
%     respAll = [resp_ON FD_conj{rfi} FT_Dur_comp{rfi} resp_OFF];
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(respAll,0.5,0.05);
    mOI = mCI; semOI = semCI;
%     mOI = mean(uni,3); semOI = std(uni,[],3)/sqrt(5);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:)]);%;semOI(:)]);    
    minI = min([mOI(:)]);%semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = {'A-On','3-T-Dur','4-T-Dur','5-T-Dur','3-T-Dis','4-T-Dis','5-T-Dis','3-T-Mix','4-T-Mix','5-T-Mix',...
        '3-I-Dur','4-I-Dur','5-I-Dur','3-I-Dis','4-I-Dis','5-I-Dis','3-I-Mix','4-I-Mix','5-I-Mix','A-Off'};
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
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',rfi),600);
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

%% distribution of dzMI for trials on the belt OR intertrials
while 1
    %%
    choose_Trial = 1;
    trialsR = [60 100];
    mean_dzMI = [];
    FD_Prop = dzMI_FD.diff_T_D; FT_Prop = dzMI_FT.diff_T_D; 
    propsD = get_props_Rs(RsDt,trialsR); propsT = get_props_Rs(RsTi,trialsR);
    for ani = 1:5
        mean_dzMIC = [];
        for cni = 1:3
            if choose_Trial
                tplbin = propsD.peak_location_bin{ani,cni};
            else
                tplbin = propsT.peak_location_bin{ani,cni};
            end
            binNs = unique(tplbin); binNs = binNs(~isnan(binNs));
            if choose_Trial
                TD = FD_Prop(ani,cni);
            else
                TD = FT_Prop(ani,cni);
            end
            for bni = 1:length(binNs)
                cell_resp = {tplbin == binNs(bni)};
                mean_dzMIC = [mean_dzMIC exec_fun_on_cell_mat(TD,'nanmean',cell_resp)];
            end
        end
        mean_dzMI = [mean_dzMI;mean_dzMIC];
    end
    
    [within,dvn,xlabels] = make_within_table({'Cond','Bins'},[3,4]);
    dataT = make_between_table({mean_dzMI},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Bins','hsd'},[1.5 1 1]);
    xdata = make_xdata([4],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(10:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.04,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.025);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'B1','B2','B3','B4'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.5 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]})
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_Bins','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([4 4 4],[1 2]);
    hf = get_figure(5,[8 7 2.5 1]);
    tcolors = repmat(mData.dcolors(1:4),1,3);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.04,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.025);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'B1','B2','B3','B4'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]})
    %%
    break;
end


%% distribution of percentage of cells on the belt OR intertrials
while 1
    %%
    rfi = 2;
    choose_Trial = 0;
    trialsR = [60 100];
    mean_dzMI = [];
    FD_Prop = dzMI_FD.diff_T_D; FT_Prop = dzMI_FT.diff_T_D; 
    propsD = get_props_Rs(RsDt,trialsR); propsT = get_props_Rs(RsTi,trialsR);
    for ani = 1:5
        mean_dzMIC = [];
        for cni = 1:3
            if choose_Trial
                tplbin = propsD.peak_location_bin{ani,cni}; good_FR_T = FD_conj{rfi}(ani,cni);%propsD.good_FR(ani,cni);
            else
                tplbin = propsT.peak_location_bin{ani,cni}; good_FR_T = FT_Dur_comp{rfi}(ani,cni);%propsD.good_FR(ani,cni);
            end
            binNs = unique(tplbin); binNs = binNs(~isnan(binNs));
            if choose_Trial
                TD = FD_Prop(ani,cni);
            else
                TD = FT_Prop(ani,cni);
            end
            for bni = 1:length(binNs)
                cell_resp = cell_list_op({tplbin == binNs(bni)},good_FR_T,'and');
%                 cell_resp = {tplbin == binNs(bni)};
                mean_dzMIC = [mean_dzMIC find_percent(cell_resp)];
%                 mean_dzMIC = [mean_dzMIC exec_fun_on_cell_mat(TD,'nanmean',cell_resp)];
            end
        end
        mean_dzMI = [mean_dzMI;mean_dzMIC];
    end
    
    [within,dvn,xlabels] = make_within_table({'Cond','Bins'},[3,2]);
    dataT = make_between_table({mean_dzMI},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Bins','hsd'},[1.5 1 1]);
    xdata = make_xdata([2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(10:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.025);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'B1','B2','B3','B4'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.5 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]})
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_Bins','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([2 2 2],[1 2]);
    hf = get_figure(5,[8 7 2.5 1]);
    tcolors = repmat(mData.dcolors(1:4),1,3);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.04,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.025);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'B1','B2','B3','B4'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]})
    %%
    break;
end

%% Speed Figure
while 1
    Rs = [RsDt(:,1) RsTi(:,1) RsDt(:,2) RsTi(:,2) RsDt(:,3) RsTi(:,3)];
    for ii = 1:length(ei)
        b1 = ei{ii}.b;
        for jj = 1:10
            alds(ii,jj) = b1.dist(b1.stim_r(jj+10)) - b1.dist(b1.air_puff_r(jj+20));
        end
    end
    ald = round(mean(alds(:)));
     ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 6],'spaceRowsCols',[0 0],'rightUpShifts',[-0.05 0.13],'widthHeightAdjustment',[-3 -350]);
    set(gcf,'color','w'); set(gcf,'Position',[5 5 5 1]);
    [Y,E] = discretize(1:49,3);
    all_speeds = []; cTxt = {'3-Trials','3-Intertrials','4-Trials','4-Intertrials','5-Trials','5-Intertrials'}; 
    for cn = 1:6
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
        xs = Rs{1,cn}.xs; N = length(xs);
        mspeed = mean(mean_speed_over_trials(:,1:N)); semspeed = std(mean_speed_over_trials(:,1:N))/sqrt(5);
        plot(xs,mspeed);
        shadedErrorBar(xs,mspeed,semspeed);
        changePosition(gca,[0.1 0.15 -0.05 -0.15]);
        if cn == 1
            put_axes_labels(gca,{'',[0 0 0]},{{'Speed (cm/sec)'},[0 -1 0]});
        end
        if ismember(cn,[1 3 5])
            if cn == 3
                plot([ald ald],[0 30],'m:','linewidth',0.5);
            end
%             plot([50 50],[0 30],'k--','linewidth',0.25);
%             plot([100 100],[0 30],'k--','linewidth',0.25);
%             bTxt = {'dB1','dB2','dB3'}; 
%             xbTxt = [25 75 125]-7; ybTxt = 31;
%             for ii = 1:length(bTxt)
%                 text(xbTxt(ii),ybTxt,bTxt{ii},'FontSize',5);
%             end
            xlabel('Distance (cm)');
        else
%             plot([5 5],[0 30],'k--','linewidth',0.25);
%             plot([10 10],[0 30],'k--','linewidth',0.25);
%             bTxt = {'tB1','tB2','tB3'}; xbTxt = [2.5 7.5 12.5]-1; ybTxt = 31;
%             for ii = 1:length(bTxt)
%                 text(xbTxt(ii),ybTxt,bTxt{ii},'FontSize',5);
%             end
            xlabel('Time (sec)');
        end
        text(xbTxt(1),ybTxt+5,cTxt{cn},'FontSize',5);
        ylim([0 30]);
        box off;
        format_axes(gca);
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('speeds_345'),600);
%%
    all_speeds_Trials = all_speeds(:,[1:3 7:9 13:15]+0);
    [within,dvn,xlabels] = make_within_table({'Conds','Bins'},[3,3]);
    dataT = make_between_table({all_speeds_Trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Conds_by_Bins','hsd'},[1 0.5 1]);
    h(h==1) = 0;
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors(7:9); tcolors = repmat(tcolors,3,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 37]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[-0.01 0.02 0.05 -0.05])
    put_axes_labels(gca,{[],[0 0 0]},{'Avg. Speed (cm/sec)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('avg_speed_anova.pdf'),600);
    ra.ranova
    ra.mauchly
    %%
    all_speeds_Trials = all_speeds(:,[1:3 7:9 13:15]+0);
    [within,dvn,xlabels] = make_within_table({'Conds','Bins'},[3,3]);
    dataT = make_between_table({all_speeds_Trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Bins','hsd'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors(7:9);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 37]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'yticklabels',[]); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.35 -0.05])
%     put_axes_labels(gca,{[],[0 0 0]},{'Avg. Speed (cm/sec)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('avg_speed_anova_pooled.pdf'),600);
    %%
    all_speeds_Trials = all_speeds(:,[7:9]+0);
    [within,dvn,xlabels] = make_within_table({'Bins'},[3]);
    dataT = make_between_table({all_speeds_Trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Bins','hsd'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors(7:9);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 37]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'yticklabels',[]); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.35 -0.05])
%     put_axes_labels(gca,{[],[0 0 0]},{'Avg. Speed (cm/sec)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('avg_speed_anova_pooled.pdf'),600);
    %%
    all_speeds_Trials = all_speeds(:,[1:3 7:9 13:15]+3);
    [within,dvn,xlabels] = make_within_table({'Conds','Bins'},[3,3]);
    dataT = make_between_table({all_speeds_Trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Bins','hsd'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors([10 11 12]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 37]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'tB1','tB2','tB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.35 -0.05])
    put_axes_labels(gca,{[],[0 0 0]},{'Avg. Speed (cm/sec)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('avg_speed_anova_pooled_IT.pdf'),600);
    %%
    break;
end