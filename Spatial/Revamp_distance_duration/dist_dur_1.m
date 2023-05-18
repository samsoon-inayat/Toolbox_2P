function dist_dur_1

%% overall population comparisons
ntrials = 30; 
si = [Ar_t_T ArL_t_T Ars_t_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T Ar_i_D ArL_i_D Ars_i_D];
%     si = [C1_i_T C2_i_T C3_i_T C4_i_T];
Rs_C = o.Rs(:,si);mRs_C = o.mR(:,si);
props_C = get_props_Rs(Rs_C,ntrials);
sel_pop_C = [dur_cells_T dis_cells_T dur_cells_I dis_cells_I];
%%
params = {'perc','N_Resp_Trials','zMI','rs','nan_zMI','nan_rs','HaFD','HiFD','PWs','centers','peak_locations'};
varT = 4;%:length(params)
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
[within,dvn,xlabels,awithinD] = make_within_table({'TI','CT','Cond'},[2,2,3]);
dataT = make_between_table({varC},dvn);
ra = RMA(dataT,within,{0.05,{'bonferroni','hsd'}});
ra.ranova
print_for_manuscript(ra)
   %%
    resp = [dur_cells_T dis_cells_T dur_cells_I dis_cells_I];
    per_resp = find_percent(resp);

    [within,dvn,xlabels] = make_within_table({'TI','CT','Cond'},[2,2,3]);
    dataT = make_between_table({per_resp},dvn);
    ra = RMA(dataT,within,{'bonferroni'});
    ra.ranova
    print_for_manuscript(ra)

    any_cells = cell_list_op(resp,[],'or',1);
    per_active_any = 100*exec_fun_on_cell_mat(any_cells,'sum')./exec_fun_on_cell_mat(any_cells,'length');
    any_cells = cell_list_op(resp,[],'and',1);
    per_active_all = 100*exec_fun_on_cell_mat(any_cells,'sum')./exec_fun_on_cell_mat(any_cells,'length');

    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
    

    %% TI_CT for responsivity of Dur, Dis, and Ind cells
switch varT
  case 1
    MY = 50; mY = 0; ysp = 7; titletext = 'Responsiveness'; y_label = '% of cells';
  case 2
    MY = 100; mY = 0; ysp = 10; titletext = 'Response Fidelity'; y_label = '% of trials';
  case 3
    MY = 3; mY = 0; ysp = 0.4; titletext = 'Mutual Information'; y_label = 'z-score';
  case 4
    MY = 1.1; mY = 0; ysp = 0.12; titletext = 'Goodness-of-Fit'; y_label = 'r-squared';
end
ff = makeFigureRowsCols(107,[10 5 1.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.36],'widthHeightAdjustment',[10 -510]);
stp = 0.3*magfac; widths = ([0.95 0.55 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = repmat(mData.colors(4:5),1,6);

axes(ff.h_axes(1,1));

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_CT','bonferroni'},[1.5 1 1]);
xdata = make_xdata([2 2],[1 1.5]);
%     hf = get_figure(5,[8 7 1.5 1]);
% tcolors = repmat(mData.colors(1:3),1,2);
MmVar = max(mVar);
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.015);
maxY = maxY + 0;
ylims = ylim;
format_axes(gca);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);
xticks = xdata; xticklabels = {'Time','Dist'};
make_bars_hollow(hbs(3:end))
set(gca,'xtick',xticks,'xticklabels',xticklabels); %xtickangle(45)
%     changePosition(gca,[0.04 0.01 -0.1 0]); 
put_axes_labels(gca,{[],[0 0 0]},{y_label,[0 0 0]});
ht = set_axes_top_text_no_line(gcf,gca,titletext,[-0.05 -0.01 0.2 0]);set(ht,'FontWeight','Bold');
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Air','No-Air'});
format_axes(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% TI_CT_Cond
switch varT
  case 1
    MY = 30; mY = 0; ysp = 3; titletext = 'Responsiveness'; y_label = '% of cells';
  case 2
    MY = 100; mY = 0; ysp = 10; titletext = 'Response Fidelity'; y_label = '% of trials';
  end
  ff = makeFigureRowsCols(107,[10 5 2.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.36],'widthHeightAdjustment',[10 -510]);
  stp = 0.25*magfac; widths = ([2.2 0.55 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
  adjust_axes(ff,[mY MY],stp,widths,gap,{''});
  tcolors = repmat(mData.colors(4:5),1,6);
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_CT_Cond','bonferroni'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3 3],[1 1.5]);
%     hf = get_figure(5,[8 7 3.5 1]);
    tcolors = repmat(mData.colors(1:9),1,2);
        tcolors = repmat(mData.colors(7:9),5,1);

    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    make_bars_hollow(hbs(7:end))
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'C3','C4','C5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); %xtickangle(45)
%     changePosition(gca,[-0.06 0.01 0.1 -0.1]); put_axes_labels(gca,{[],[0 0 0]},{'z-score',[0 0 0]});
    put_axes_labels(gca,{[],[0 0 0]},{'z-score',[0 0 0]});
    ht = set_axes_top_text_no_line(gcf,gca,'Mutual Information',[0 0 0 0]);set(ht,'FontWeight','Bold');
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Time','Dist','Time','Dist'},{[0,0.03]});
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,6,{'Air','No-Air'},{[-0.1 00]});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graph.pdf'),600);
    
%% CT
ff = makeFigureRowsCols(107,[10 5 1.55 1],'RowsCols',[1 2],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.36],'widthHeightAdjustment',[10 -510]);
  MY = 4; mY = 0; ys = 0.52;
  stp = 0.27*magfac; widths = ([0.55 0.55 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
  adjust_axes(ff,[mY MY],stp,widths,gap,{''});
  tcolors = repmat(mData.colors(1:2),1,6);
  axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([2],[1 1.5]);
%     hf = get_figure(5,[8 7 1.25 1]);
    tcolors = repmat(mData.colors(4:6),2,1);
    MmVar = max(mVar);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ys,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);
    xticks = xdata; xticklabels = {'TE','DE','IE'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);% xtickangle(45)
put_axes_labels(gca,{[],[0 0 0]},{'z-score',[0 0 0]});
ht = set_axes_top_text_no_line(gcf,gca,'Mutual Information',[-0.05 -0.01 0.2 0]);set(ht,'FontWeight','Bold');
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Air','No-Air'});
format_axes(gca);
% save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

  axes(ff.h_axes(1,2));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([2],[1 1.5]);
%     hf = get_figure(5,[8 7 1.25 1]);
    tcolors = repmat(mData.dcolors(7:9),2,1);
    MmVar = max(mVar);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ys,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Air','No Air'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30)

% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Air','No-Air'});
format_axes(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% running RMANOVA on all neurons in one animal (animal wise)
an = 4;
pop_var_name = {'all'};
sel_pop_C = cell_list_op(props_C,pop_var_name); 
var_C = get_vals(props_C.zMI,sel_pop_C);

var_CM = [];
for cc = 1:size(var_C,2)
  var_CM = [var_CM var_C{an,cc}];
end

[within,dvn,xlabels,awithinD] = make_within_table({'TI','DT','Cond'},[2,2,3]);
dataT = make_between_table({var_CM},dvn);
ra = RMA(dataT,within,{0.05,{'bonferroni'}});
ra.ranova



%%
var_C1 = [dzMI_FD.diff_T_D dzMI_FT.diff_T_D];
sel_pop_dzMI_t = cell_list_op(sel_pop_C(:,1:3),sel_pop_C(:,4:6),'or');
sel_pop_dzMI_i = cell_list_op(sel_pop_C(:,7:9),sel_pop_C(:,10:12),'or');
var_C = get_vals(var_C1,[sel_pop_dzMI_t sel_pop_dzMI_i]);
mean_var_C1 = exec_fun_on_cell_mat(var_C,'nanmean');
[within,dvn,xlabels,awithinD] = make_within_table({'TI','Cond'},[2,3]);
dataT = make_between_table({mean_var_C1},dvn);
ra = RMA(dataT,within,{0.05,{'bonferroni'}});
ra.ranova

%%
sel_pop_C = [dzMI_FD.resp_T_g_D dzMI_FD.resp_D_g_T dzMI_FT.resp_T_g_D dzMI_FT.resp_D_g_T];

params = {'perc','N_Resp_Trials','zMI','rs','nan_zMI','nan_rs','HaFD','HiFD','PWs','centers','peak_locations'};
varT = 3;%:length(params)
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
ra = RMA(dataT,within,{0.05,{'bonferroni'}});
ra.ranova




%% responsivity
    % for one RF
    cni = 1:3;
    rfi = 4;
    resp = [FD_Dur_comp{rfi}(:,cni) FD_Dis_comp{rfi}(:,cni) FD_conj{rfi}(:,cni) FT_Dur_comp{rfi}(:,cni) FT_Dis_comp{rfi}(:,cni) FT_conj{rfi}(:,cni)];
%     rfi = 2;
%     resp = [resp FD_Dur_comp{rfi}(:,cni) FD_Dis_comp{rfi}(:,cni) FD_conj{rfi}(:,cni) FT_Dur_comp{rfi}(:,cni) FT_Dis_comp{rfi}(:,cni) FT_conj{rfi}(:,cni)];

    per_resp = 100*exec_fun_on_cell_mat(resp,'sum')./exec_fun_on_cell_mat(resp,'length');

    [withinD,dvn,xlabels] = make_within_table({'TI','CT','Cond'},[2,3,3]);
    dataT = make_between_table({per_resp},dvn);
    ra = RMA(dataT,withinD,{'bonferroni'});
    ra.ranova

    any_cells = cell_list_op(resp,[],'or',1);
    per_active_any = 100*exec_fun_on_cell_mat(any_cells,'sum')./exec_fun_on_cell_mat(any_cells,'length');
    any_cells = cell_list_op(resp,[],'and',1);
    per_active_all = 100*exec_fun_on_cell_mat(any_cells,'sum')./exec_fun_on_cell_mat(any_cells,'length');

    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
    %% subsequent ANOVAs to figure out clarity about what is happening
    alpha = 0.05/2;
    dataT_T = dataT(:,(awithinD(:,1) == 2));
    [within,dvn,xlabels] = make_within_table({'CT','Cond'},[2,3]);
    ra = RMA(dataT_T,within,{alpha,{'bonferroni'}});
    ra.ranova
    %%
    alpha = 0.05/2;
    dataT_T = dataT(:,(awithinD(:,2) == 2));
    [within,dvn,xlabels] = make_within_table({'TI','Cond'},[2,3]);
    ra = RMA(dataT_T,within,{alpha,{'bonferroni'}});
    ra.ranova
    %%
    alpha = 0.05/6;
    for tii = 1:2
      for cii = 1:3
        dataT_T_C1 = dataT(:,(withinD{:,1} == tii) & (withinD{:,3} == cii));
        [within,dvn,xlabels] = make_within_table({'CT'},[3]);
        ra = RMA(dataT_T_C1,within,{});
    %     ra = RMA(dataT_T_C1,within,{alpha,{'bonferroni'}});
        p_v(tii,cii) = ra.ranova{3,9};
      end
    end
    
    alpha = 0.05/3;
    dataT_C = dataT(:,(withinD{:,3} == 3));
    [within,dvn,xlabels] = make_within_table({'TI','CT'},[2,3]);
    ra = RMA(dataT_C,within,{alpha,{'bonferroni'}});
    ra.ranova
    
    %% not using
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


%% check the difference in zMI for Dis and Dur for the different cell types DurC, DisC, and DDM
while 1

    %%
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(7:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.75,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1)-1 maxY+1]); format_axes(gca);
    xticks = xdata; xticklabels = {'Ti','Di','TrDr'};
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


%% for only one condition e.g., 3 check the difference in zMI for Dis and Dur for the different cell types DurC, DisC, and DDM
while 1
    %%
    cni = 1;
    mean_dzMI = [];
    FD_Prop = dzMI_FD.diff_T_D; FT_Prop = dzMI_FT.diff_T_D; 
%     FD_Prop = dzMI_FD.rs.diff_T_D; FT_Prop = dzMI_FT.rs.diff_T_D; 
%     FD_Prop = dzMI_FD.HaFD.diff_T_D; FT_Prop = dzMI_FT.HaFD.diff_T_D; 
%     FD_Prop = dzMI_FD.HiFD.diff_T_D; FT_Prop = dzMI_FT.HiFD.diff_T_D; 
    for rfi = 1:2
        TD = FD_Prop(:,cni);
        cell_resp = FD_Dur_comp{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp(:,cni))];
        cell_resp = FD_Dis_comp{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp(:,cni))];
        cell_resp = FD_conj{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp(:,cni))];

        TD = FT_Prop(:,cni);
        cell_resp = FT_Dur_comp{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp(:,cni)) ];
        cell_resp = FT_Dis_comp{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp(:,cni))];
        cell_resp = FT_conj{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp(:,cni))];
    end
    
%     [within,dvn,xlabels] = make_within_table({'TI','CT'},[2,3]);
%     dataT = make_between_table({mean_dzMI},dvn);
%     ra = RMA(dataT,within,{'hsd'});
%     ra.ranova
    
    [within,dvn,xlabels] = make_within_table({'RF','TI','CT'},[2,2,3]);
    dataT = make_between_table({mean_dzMI},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'RF_TI_CT','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3 3],[1 1.5]);
    hf = get_figure(5,[8 7 3.5 1]);
    tcolors = mData.dcolors(7:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.015,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dur','Dis','Mix'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('dzMI_CT.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_CT','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3],[1 1.5]);
    hf = get_figure(5,[8 7 3.5 1]);
    tcolors = mData.dcolors(7:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.015,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dur','Dis','Mix'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('dzMI_CT.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1.5 1 1]);
    xdata = make_xdata([3],[1 1.5]);
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
    rfi = 2;
    respAll = [FD_Dur_comp{rfi} FD_Dis_comp{rfi} FD_conj{rfi} FT_Dur_comp{rfi} FT_Dis_comp{rfi} FT_conj{rfi}];
    respAll = [resp_ON FD_conj{rfi} FT_Dur_comp{rfi} resp_OFF];
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
    txl = {'A-On','3-T-Mix','4-T-Mix','5-T-Mix',...
        '3-I-Dur','4-I-Dur','5-I-Dur','A-Off'};
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
