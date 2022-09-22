function dist_dur_2

%% overall population comparisons
ntrials = 30; 
si = [Ar_t_D ArL_t_D Ars_t_D Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T Ar_i_T ArL_i_T Ars_i_T ];
si = [Ar_t_D ArL_t_D Ars_t_D Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T Ar_i_T ArL_i_T Ars_i_T ];
% si = [Ar_t_T ArL_t_T Ars_t_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T Ar_i_D ArL_i_D Ars_i_D];
Rs_C = o.Rs(:,si);mRs_C = o.mR(:,si);
props_C = get_props_Rs(Rs_C,ntrials);
sel_pop_C = [dur_cells_T dis_cells_T dur_cells_I dis_cells_I];
%%
params = {'perc','N_Resp_Trials','zMI','rs','nan_zMI','nan_rs','HaFD','HiFD','PWs','centers','peak_locations'};
varT = 11;%:length(params)
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
%%
varC = mean_var_C;
[within,dvn,xlabels,awithinD] = make_within_table({'TI','CT','Cond'},[2,2,3]);
dataT = make_between_table({varC},dvn);
ra = RMA(dataT,within,{0.05,{'bonferroni','hsd'}});
ra.ranova
%%
varC = mean_var_C(:,1:6);
[within,dvn,xlabels,awithinD] = make_within_table({'CT','Cond'},[2,3]);
dataT = make_between_table({varC},dvn);
ra = RMA(dataT,within,{0.05,{'bonferroni','hsd'}});
ra.ranova
%%
varC = mean_var_C(:,7:12);
[within,dvn,xlabels,awithinD] = make_within_table({'CT','Cond'},[2,3]);
dataT = make_between_table({varC},dvn);
ra = RMA(dataT,within,{0.05,{'bonferroni','hsd'}});
ra.ranova

%% CT zMI
ff = makeFigureRowsCols(107,[10 5 1.255 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.36],'widthHeightAdjustment',[10 -510]);
  MY = 120; mY = 0; ys = 15;
  stp = 0.27*magfac; widths = ([0.55 0.55 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
  adjust_axes(ff,[mY MY],stp,widths,gap,{''});
  tcolors = repmat(mData.colors(4:5),1,6);
  axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([2],[1 1.5]);
%     hf = get_figure(5,[8 7 1.25 1]);
    tcolors = repmat(mData.colors(4:9),2,1);
    MmVar = max(mVar);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ys,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Time','Dist'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);% xtickangle(45)
put_axes_labels(gca,{[],[0 0 0]},{'cm',[0 0 0]});
ht = set_axes_top_text_no_line(gcf,gca,'Field Center',[-0.05 -0.01 0.2 0]);set(ht,'FontWeight','Bold');
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'No-Air'});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Air'});
format_axes(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);