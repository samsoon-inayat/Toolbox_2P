function distance_duration_pops

%%
ntrials = [50];
RsDt = o.Rs(:,[Ar_t_D ArL_t_D Ars_t_D]);  RsTt = o.Rs(:,[Ar_t_T ArL_t_T Ars_t_T]);
RsDi = o.Rs(:,[Ar_i_D ArL_i_D Ars_i_D]);  RsTi = o.Rs(:,[Ar_i_T ArL_i_T Ars_i_T]);
[dzMI_FD,dzMI_FT] = get_zMI_comp_dist_time(RsDt,RsTt,RsDi,RsTi);

[allRs_FD,allmRs_FD] = get_trial_Rs(o,[Ar_t_D ArL_t_D Ars_t_D],1:10);
[allRs_FT,allmRs_FT] = get_trial_Rs(o,[Ar_i_T ArL_i_T Ars_i_T],1:10);
[allRs_Ab,allmRs_Ab] = get_trial_Rs(o,[Ab_T Abs_T],1:10);
respDTC = combine_distance_time_rasters(o.Rs(:,[Ar_t_D ArL_t_D Ars_t_D]),o.Rs(:,[Ar_i_T ArL_i_T Ars_i_T]),ntrials);


mRsDt = o.mR(:,[Ar_t_D ArL_t_D Ars_t_D]);  mRsTt = o.mR(:,[Ar_t_T ArL_t_T Ars_t_T]);
mRsDi = o.mR(:,[Ar_i_D ArL_i_D Ars_i_D]);  mRsTi = o.mR(:,[Ar_i_T ArL_i_T Ars_i_T]);

Rs_Ab = o.Rs(:,[Ab_T Abs_T]); mRs_Ab = o.mR(:,[Ab_T Abs_T]);

trialsR = [50];
scale = 1000;
props_Ab = get_props_Rs(Rs_Ab,trialsR); propsD = get_props_Rs(RsDt,trialsR,scale); propsT = get_props_Rs(RsTi,trialsR,scale);
propsDi = get_props_Rs(RsDi,trialsR,scale); propsTt = get_props_Rs(RsTt,trialsR,scale);
%%

all_responsive_cells = cell_list_op(propsD.vals,propsT.vals,'or');
all_responsive_cells = cell_list_op(all_responsive_cells,propsTt.vals,'or');
all_responsive_cells = cell_list_op(all_responsive_cells,propsDi.vals,'or');

trialsR = [10 30];
props_Ab13 = get_props_Rs(Rs_Ab,trialsR); propsD13 = get_props_Rs(RsDt,trialsR); propsT13 = get_props_Rs(RsTi,trialsR);

trialsR = [40 70];
props_Ab47 = get_props_Rs(Rs_Ab,trialsR); propsD47 = get_props_Rs(RsDt,trialsR); propsT47 = get_props_Rs(RsTi,trialsR);

trialsR = [80 100];
props_Ab810 = get_props_Rs(Rs_Ab,trialsR); propsD810 = get_props_Rs(RsDt,trialsR); propsT810 = get_props_Rs(RsTi,trialsR);

%%
while 1
FD_Dis = propsD.vals; FD_Dur = propsTt.vals;
FT_Dis = propsDi.vals; FT_Dur = propsT.vals;

for rr = 1:size(FD_Dis,1)
    for cc = 1:size(FD_Dis,2)
        cellP1 = FD_Dis(rr,cc); cellP2 = FD_Dur(rr,cc);
        FD_Dis_comp(rr,cc) = cell_list_op(cellP1,cell_list_op(cellP2,[],'not'),'and');
        FD_Dur_comp(rr,cc) = cell_list_op(cell_list_op(cellP1,[],'not'),cellP2,'and');
        FD_conj(rr,cc) = cell_list_op(cellP1,cellP2,'and');
        
        cellP1 = FT_Dis(rr,cc); cellP2 = FT_Dur(rr,cc);
        FT_Dis_comp(rr,cc) = cell_list_op(cellP1,cell_list_op(cellP2,[],'not'),'and');
        FT_Dur_comp(rr,cc) = cell_list_op(cell_list_op(cellP1,[],'not'),cellP2,'and');
        FT_conj(rr,cc) = cell_list_op(cellP1,cellP2,'and');
    end
end
%%
resp = [FD_Dis_comp FD_Dur_comp FD_conj FT_Dis_comp FT_Dur_comp FT_conj];

per_resp = 100*exec_fun_on_cell_mat(resp,'sum')./exec_fun_on_cell_mat(resp,'length');

[within,dvn,xlabels] = make_within_table({'TI','CT','Cond'},[2,3,3]);
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
break;
end

%% graph percentages
while 1
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1.5 1 1]);
    xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(7:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',4,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 30]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 5 10 15 20]); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.3 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
    xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(10:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',4,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 30]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 5 10 15 20]); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.3 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
%%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_CT','hsd'},[1.5 1 1]);
     h(h==1) = 0;
     p = ones(size(p)); 
    ind = ismember(combs,[2 3],'rows');  p(ind) = ra.MC.hsd.CT_by_TI.pValue(4); 
    ind = ismember(combs,[4 5],'rows'); p(ind) = ra.MC.hsd.CT_by_TI.pValue(7); 
    ind = ismember(combs,[4 6],'rows'); p(ind) = ra.MC.hsd.CT_by_TI.pValue(8); 
    ind = ismember(combs,[5 6],'rows'); p(ind) = ra.MC.hsd.CT_by_TI.pValue(10); 
    h = p<0.05;
    xdata = make_xdata([3 3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%)',round(mra),round(semra)),'FontSize',6);
%     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 30]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 5 10 15 20]); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);

    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_CT_Cond','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3 3 3 3],[1 1.5]);
    hf = get_figure(5,[8 7 3.5 1]);
    tcolors = repmat(mData.colors(1:9),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0.75,maxY-3,sprintf('Any Condition or Type (%d\x00B1%d%%)',round(mra),round(semra)),'FontSize',6);
%     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 30]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 5 10 15 20]); xtickangle(45)
    changePosition(gca,[-0.04 0.01 0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
break;
end
%% compare zMI
while 1
    ntrials = 50; 
    si = [Ar_t_D ArL_t_D Ars_t_D Ar_t_D ArL_t_D Ars_t_D Ar_t_D ArL_t_D Ars_t_D Ar_t_D ArL_t_D Ars_t_D Ar_t_D ArL_t_D Ars_t_D Ar_t_D ArL_t_D Ars_t_D];
    si = [Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ars_i_D Ar_i_D ArL_i_D Ar_i_T ArL_i_T Ars_i_T];
%     si = [Ar_i_T ArL_i_T Ars_i_T Ar_i_T ArL_i_T Ars_i_T Ar_i_T ArL_i_T Ars_i_T Ar_i_T ArL_i_T Ars_i_T Ar_i_T ArL_i_T Ars_i_T Ar_i_T ArL_i_T Ars_i_T];
%     si = [Ar_i_T ArL_i_T Ars_i_T];
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    gFR = [FD_Dis_comp FD_Dur_comp FT_Dis_comp FT_Dur_comp]; rf = props1.zMI;
    all_gFR = exec_fun_on_cell_mat(rf,'nanmean',gFR);

    [within,dvn,xlabels] = make_within_table({'TI','CT','Cond'},[2,2,3]);
    dataT = make_between_table({all_gFR},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_CT_Cond','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3 3],[1 1.5]);
    hf = get_figure(5,[8 7 3.5 1]);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY+1;
    ylims = ylim;
    format_axes(gca);
    htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%)',round(mra),round(semra)),'FontSize',6);
%     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 1 2 3]); xtickangle(45)
    changePosition(gca,[0.01 0.01 -0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT_by_TI','hsd'},[1.5 1 1]);
    p = ones(size(p)); 
    ind = ismember(combs,[1 2],'rows');  p(ind) = ra.MC.hsd.CT_by_TI.pValue(1); 
    ind = ismember(combs,[3 4],'rows'); p(ind) = ra.MC.hsd.CT_by_TI.pValue(end); 
    h = p<0.05;
    xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.015);
    maxY = maxY + 1;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'T','I'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.01 0.01 -0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Z-score',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
    xdata = make_xdata([2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY ;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dis','Dur','Conj'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 2.5 5]); xtickangle(45)
    changePosition(gca,[0.15 0.01 -0.4 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Z-Score'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end


%% compare zMI conj cells
while 1
    ntrials = 50; 
    si = [Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ars_i_D Ar_i_D ArL_i_D Ar_i_T ArL_i_T Ars_i_T];
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    gFR = [FD_conj FD_conj FT_conj FT_conj]; rf = props1.zMI;
    all_gFR = exec_fun_on_cell_mat(rf,'nanmean',gFR);

    [within,dvn,xlabels] = make_within_table({'TI','DT','Cond'},[2,2,3]);
    dataT = make_between_table({all_gFR},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_CT_Cond','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
    xdata = make_xdata([3 3],[1 1.5]);
    hf = get_figure(5,[8 7 3.5 1]);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY+1;
    ylims = ylim;
    format_axes(gca);
    htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%)',round(mra),round(semra)),'FontSize',6);
%     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 1 2 3]); xtickangle(45)
    changePosition(gca,[0.01 0.01 -0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_DT','hsd'},[1.5 1 1]);
    p = ones(size(p)); 
    ind = ismember(combs,[1 2],'rows');  p(ind) = ra.MC.hsd.TI_by_DT.pValue(1); 
    ind = ismember(combs,[3 4],'rows'); p(ind) = ra.MC.hsd.TI_by_DT.pValue(3); 
    ind = ismember(combs,[5 6],'rows'); p(ind) = ra.MC.hsd.TI_by_DT.pValue(end); 
    h = p<0.05;
    xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.015);
    maxY = maxY + 1;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'T','I'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.01 0.01 -0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
    xdata = make_xdata([2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY ;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dis','Dur','Conj'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 2.5 5]); xtickangle(45)
    changePosition(gca,[0.15 0.01 -0.4 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Z-Score'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end


%%
an = 5; cn = 3;
diffzMI = [dzMI_FD.diff_D_T{an,cn} dzMI_FT.diff_D_T{an,cn}];
bins = -11:0.5:11;
dists = (hist(diffzMI,bins));
dists = cumsum(dists./sum(dists));
figure(1000);clf;plot(bins,dists);
hold on;plot([0 0],ylim,'g');
% legend('FD','FT','','northwest');
xlim([bins(1) bins(end)]);

%% compare percent responsive cells
while 1
    cells_list = [dzMI_FD.resp_D_g_T dzMI_FD.resp_T_g_D dzMI_FD.resp_complex dzMI_FT.resp_D_g_T dzMI_FT.resp_T_g_D dzMI_FT.resp_complex];
    
    per_active = 100*exec_fun_on_cell_mat(cells_list,'sum')./exec_fun_on_cell_mat(cells_list,'length');
    
    [within,dvn,xlabels] = make_within_table({'TI','CT','Cond'},[2,3,3]);
    dataT = make_between_table({per_active},dvn);
    ra = RMA(dataT,within);
    ra.ranova

    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_CT','hsd'},[1.5 1 1]);
    ind = ismember(combs,[1 2],'rows');    p = ones(size(p)); p(ind) = ra.MC.hsd.TI_by_CT.pValue(1); 
    ind = ismember(combs,[3 4],'rows'); p(ind) = ra.MC.hsd.TI_by_CT.pValue(3); p(end) = ra.MC.hsd.TI_by_CT.pValue(end); h = p<0.05;
    xdata = make_xdata([2 2 2],[1 1.5]);
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
    changePosition(gca,[0.02 0.01 -0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
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



%% visualizing different cell groups for the same condition
while 1
    %%
    an = 3; cn = 1;
    
    tRs = [RsDt(an,cn) allRs_FD{cn}(an,:)];  tmRs = [mRsDt(an,cn) allmRs_FD{cn}(an,:)];
    tRsi = [RsTi(an,cn) allRs_FT{cn}(an,:)];  tmRsi = [mRsTi(an,cn) allmRs_FT{cn}(an,:)];

    resp = respDT.inh(an,cn);
%     resp = respDT.exc_and_good_FR(an,cn);
    resp = propsD13.good_FR(an,cn);
    ff = plot_pop_vec_trials(100,[1 1 14 3],tRs,tmRs,resp)
%     colormap parula
    
    resp = respDT.inh_and_good_FR(an,cn);
%     resp = respDT.exc_not_good_FR(an,cn);
    resp = propsD47.good_FR(an,cn);
    ff = plot_pop_vec_trials(101,[1 4 14 3],tRsi,tmRsi,resp)

    
    resp = propsD810.good_FR(an,cn);
    ff = plot_pop_vec_trials(102,[1 7 14 3],tRsi,tmRsi,resp)

    %%
    break;
end

%%
while 1
    %%
    cns = 3;
    an = 1; cn = cns(1);
    cnab = 2;
    resp = respDT.exc(an,cn);
    resp = props_Ab.good_FR(an,cnab);
    resp = propsD.good_FR(an,cn);
%     resp = dzMI_FD.resp_D_g_T_and_good_FR(an,cn);
    
    
    tRs = [RsDt(an,cn) allRs_FD{cn}(an,:)];  tmRs = [mRsDt(an,cn) allmRs_FD{cn}(an,:)];
    ff = plot_pop_vec_trials(100,[1 1 14 3],tRs,tmRs,resp)
    colormap jet
    
    tRs = [RsTi(an,cn) allRs_FT{cn}(an,:)];  tmRs = [mRsTi(an,cn) allmRs_FT{cn}(an,:)];
    ff = plot_pop_vec_trials(101,[1 4 14 3],tRs,tmRs,resp)
    colormap jet
    
    tRs = [Rs_Ab(an,cnab) allRs_Ab{cnab}(an,:)];  tmRs = [mRs_Ab(an,cnab) allmRs_Ab{cnab}(an,:)];
%     resp = respDT.exc(an,cn);
    ff = plot_pop_vec_trials(102,[1 7 14 3],tRs,tmRs,resp)
    colormap jet
    %%
    break;
end;

%% Overlap Indices ImageSC overall
while 1
%     cn = 1;
%     si = [Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
%     resp  = [respDT.inh respDT.exc propsD.good_FR propsT.good_FR]
%     resp  = [respDT.inh respDT.exc]
%     resp = [props_Ab.good_FR(:,cn),propsD.good_FR(:,cn),propsT.good_FR(:,cn)];
% %     resp = [props_Ab13.good_FR,propsD13.good_FR,propsT13.good_FR];
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(resp,0.5,0.05);
    mOI = mCI;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:)]);%semOI(:)]);    
    minI = min([mOI(:)]);%semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1;
%     ptxl = {'D','D','D','T','T','T','C','C','C','D','D','D','T','T','T','C','C','C'};
    ptxl = {'D','D','D','T','T','T','','','','','',''};
    ptxl = {'D','D','D','','',''};
%     ptxl = {'D','T','C','D','T','C'};
    for ii = 1:length(ptxl)
        txl{ii} = sprintf('%s-%s',ptxl{ii},rasterNamesTxt{si(ii)});
    end
    txl = rasterNamesTxt([Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T]);
    txl = rasterNamesTxt([Ab_T Ar_t_D Ar_i_T]);
%     txl = {'1-13','2-13','3-13','1-47','2-47','3-47','1-81','2-81','3-81'};
   
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 3 3.5 3.5]);
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
    axis equal
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
%     set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
%     text(-0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); 
    set(gca,'Ydir','normal');ytickangle(25);
    box on
    changePosition(gca,[0.0 0 -0.04 0]);
    hc = putColorBar(gca,[0.08 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.08 0.11 0.15]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean1_conjunctionP.pdf',ntrials),600);
    break;
end

%% Overlap Indices ImageSC
while 1

    mOI = mCI;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:)]);%semOI(:)]);    
    minI = min([mOI(:)]);%semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1;
%     ptxl = {'D','D','D','T','T','T','C','C','C','D','D','D','T','T','T','C','C','C'};
%     ptxl = {'D','T','C','D','T','C'};
%     for ii = 1:length(ptxl)
%         txl{ii} = sprintf('%s-%s',ptxl{ii},rasterNamesTxt{sii(ii)});
%     end

    
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 3 3.5 3.5]);
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
    axis equal
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
%     set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
%     text(-0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); 
    set(gca,'Ydir','normal');ytickangle(25);
    box on
    changePosition(gca,[0.0 0 -0.04 0]);
    hc = putColorBar(gca,[0.08 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.08 0.11 0.15]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean1_conjunctionP.pdf',ntrials),600);
    break;
end


%% Overlap Indices ImageSC
while 1
    sii = [Ar_t_D ArL_t_D Ars_t_D Ar_t_D ArL_t_D Ars_t_D Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T Ar_i_T ArL_i_T Ars_i_T Ar_i_T ArL_i_T Ars_i_T];
    resp = [FDgFR_D_g_T FDgFR_T_g_D FDgFR_Comp FTgFR_D_g_T FTgFR_T_g_D FTgFR_Comp];

%     [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(resp,0.5,0.05);
    mOI = mCI;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:)]);%semOI(:)]);    
    minI = min([mOI(:)]);%semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1;
    ptxl = {'D','D','D','T','T','T','C','C','C','D','D','D','T','T','T','C','C','C'};
    for ii = 1:length(ptxl)
        txl{ii} = sprintf('%s-%s',ptxl{ii},rasterNamesTxt{sii(ii)});
    end

    
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 3 3.5 3.5]);
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
    axis equal
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
%     set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
%     text(-0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); 
    set(gca,'Ydir','normal');ytickangle(25);
    box on
    changePosition(gca,[0.0 0 -0.04 0]);
    hc = putColorBar(gca,[0.08 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.08 0.11 0.15]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean1_conjunctionP.pdf',ntrials),600);
    break;
end


%% agglomerative hierarchical clustering
while 1
    mOI1 = mOI;
%     mOI1(isnan(mOI1)) = 1;
%     Di = pdist(mOI1);
    Di = pdist(mOI1,@naneucdist);
%     tree = linkage(mOI1,'average','euclidean');
    tree = linkage(Di,'average');
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 6.97 1]);
    set(H,'linewidth',1);
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel('Eucledian Distance');changePosition(hx,[0 -0.1 0]);
    changePosition(gca,[0.03 0.0 0.05 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster1.pdf'),600);
    %%
    break;
end
%% I want to explore the trial to trial changes when the animal is running a fixed distance
RsDC = combine_rasters_conditions(RsDt);
% plotRasters_simplest(RsDC{an,1},find(resp))

RsTC = combine_rasters_conditions(RsTi);
% plotRasters_simplest(RsTC{1,1},[])
%% combine the distance and time rasters horizontally and see individual cell rasters
RsDTC = combine_rasters_horizontally([RsDC,RsTC]);
% ccs = cell_list_op(respDT.inh(1,:),[],'or',1);
% plotRasters_simplest(RsDTC{an,1},find(ccs{an}));
%%
plotRasters_simplest(RsDTC{an,1},find(respA));
%% find mean over the 30 trials from conditions 3 4 5 and see population vector
mRDTC = calc_mean_rasters(RsDTC,[]); 
%%
an = 4;
Rs = RsDTC(an) ;mR = mRDTC(an);
ccs = cell_list_op(respDTC.resp(an,1),[],'or',1);
% ccs = cell_list_op(FT_conj(1,:),[],'or',1);
% ccs =  cell_list_op(ccs,{respA},'or');
% ccs = {respA};
ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[2 1],...
    'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.13],'widthHeightAdjustment',...
    [-50 -80]);    set(gcf,'color','w');  set(gcf,'Position',[10 3 3.5 5]);
[CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,ccs,0);
ff = show_population_vector_and_corr_combined(mData,ff,Rs,mRR,CRc,[],[]);
axes(ff.h_axes(2,1));
set(gca,'xtick',[1 49 51 100],'xticklabels',{'0','','','15'});
text(15,-4,'Distance (cm)','FontSize',6); text(65,-4,'Time (s)','FontSize',6);
text(44,-4,'150','FontSize',6); text(50,-4,'0','FontSize',6);
%%
an = 1;
rasters = RsDC{an}.sp_rasters1;
Dur = RsDC{an}.duration1;
rastersT = RsTC{an}.sp_rasters1;
DurT = RsTC{an}.duration1;
nbins = 4;
nShuffle = 500;
numCells = size(rasters,3);
parfor ii = 1:numCells
    rng(3,'twister');
    [OD(ii),~] = info_metrics_S_onlyMI(rasters(:,:,ii),[],nbins,Dur,nShuffle);
    rng(3,'twister');
    [OT(ii),~] = info_metrics_S_onlyMI(rastersT(:,:,ii),[],nbins,DurT,nShuffle);
end
%%
for ii = 1:numCells
    zMIsD(ii) = OD(ii).ShannonMI_Zsh;
    zMIsT(ii) = OT(ii).ShannonMI_Zsh;
end
%%
an = 1;
rasters = RsDC{an}.sp_rasters1;
rastersT = RsTC{an}.sp_rasters1;
rastersDT = RsDTC{an}.sp_rasters1;
numCells = size(rasters,3);
p = zeros(numCells,1);resp = logical(p);
pT = p; respT = logical(p); pDT = p; respDT = logical(p);
parfor ii = 1:numCells
    p(ii) = anova1(rasters(:,:,ii),1:size(rasters(:,:,ii),2),'off');    
    if p(ii) < 0.05
        resp(ii) = 1;
    end
    pT(ii) = anova1(rastersT(:,:,ii),1:size(rastersT(:,:,ii),2),'off');
    if pT(ii) < 0.05
        respT(ii) = 1;
    end
    pDT(ii) = anova1(rastersDT(:,:,ii),1:size(rastersDT(:,:,ii),2),'off');
    if pDT(ii) < 0.05
        respDT(ii) = 1;
    end
    
end
respA = resp | respT | respDT;
% plotRasters_simplest(RsDTC{an,1},find(respA));
%% fitting of 2D gaussian on the raster plots
statsetfitnlm = statset('fitnlm');
statsetfitnlm.MaxIter = 1000;
statsetfitnlm.TolFun = 1e-10;
% statsetfitnlm.Display = 'iter';
statsetfitnlm.TolX = statsetfitnlm.TolFun;
statsetfitnlm.UseParallel = 1;
% statsetfitnlm.RobustWgtFun = 'welsch';

%%
for ii = 1:numCells
    if ~respDT(ii)
        continue;
    end
%     raster = double(rasters(:,:,ii) > 0);
    raster = rastersDT(:,:,ii);
    rasterFilt = imgaussfilt(raster,[1 2]);
    [rasterF,mdl,coeff_rs] = do_gauss_fit2D(rasterFilt,statsetfitnlm,[1 1]);
    figure(10000);clf;
    subplot 131; imagesc(raster);set(gca,'YDir','normal')
    subplot 132; imagesc(rasterFilt);set(gca,'YDir','normal')
    subplot 133; imagesc(rasterF);set(gca,'YDir','normal')
    title(ii);
    pause(1);
%     pT(ii) = anova1(rastersT(:,:,ii),1:size(rastersT(:,:,ii),2),'off');    
%     pDT(ii) = anova1(rastersDT(:,:,ii),1:size(rastersDT(:,:,ii),2),'off');    
end
%%
respA = resp | respT | respDT;
plotRasters_simplest(RsDTC{an,1},find(respA));



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
    htxt = text(0.75,maxY-2,{'Any Condition',sprintf('(%d\x00B1%d%%)',round(mra),round(semra))},'FontSize',6);
%     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[0.04 0.01 0 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end

%% compare percent of conjunctive and complementary cells
while 1
    ntrials = 50; 
    si = [Ar_t_D ArL_t_D Ars_t_D];
    si = [Ar_i_T ArL_i_T Ars_i_T];
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

%%
while 1
FD_Dis = propsD.vals; FD_Dur = propsTt.vals;
FT_Dis = propsDi.vals; FT_Dur = propsT.vals;

for rr = 1:size(FD_Dis,1)
    for cc = 1:size(FD_Dis,2)
        cellP1 = FD_Dis(rr,cc); cellP2 = FD_Dur(rr,cc);
        FD_Dis_comp(rr,cc) = cell_list_op(cellP1,cell_list_op(cellP2,[],'not'),'and');
        FD_Dur_comp(rr,cc) = cell_list_op(cell_list_op(cellP1,[],'not'),cellP2,'and');
        FD_conj(rr,cc) = cell_list_op(cellP1,cellP2,'and');
        
        cellP1 = FT_Dis(rr,cc); cellP2 = FT_Dur(rr,cc);
        FT_Dis_comp(rr,cc) = cell_list_op(cellP1,cell_list_op(cellP2,[],'not'),'and');
        FT_Dur_comp(rr,cc) = cell_list_op(cell_list_op(cellP1,[],'not'),cellP2,'and');
        FT_conj(rr,cc) = cell_list_op(cellP1,cellP2,'and');
    end
end

resp = [FD_Dis_comp FD_Dur_comp FD_conj FT_Dis_comp FT_Dur_comp FT_conj];
resp = [FT_Dis_comp FT_Dur_comp FT_conj];

per_resp = 100*exec_fun_on_cell_mat(resp,'sum')./exec_fun_on_cell_mat(resp,'length');

[within,dvn,xlabels] = make_within_table({'CT','Cond'},[3,3]);
dataT = make_between_table({per_resp},dvn);
ra = RMA(dataT,within,{'hsd'});
ra.ranova

any_cells = cell_list_op(resp,[],'or',1);
per_active_any = 100*exec_fun_on_cell_mat(any_cells,'sum')./exec_fun_on_cell_mat(any_cells,'length');
any_cells = cell_list_op(resp,[],'and',1);
per_active_all = 100*exec_fun_on_cell_mat(any_cells,'sum')./exec_fun_on_cell_mat(any_cells,'length');

[mra,semra] = findMeanAndStandardError(per_active_any);
[mrall,semrall] = findMeanAndStandardError(per_active_all);
break;
end


%% compare percent of conjunctive and complementary cells
while 1
    resp = [FD_Dis_comp FD_Dur_comp FD_conj FT_Dis_comp FT_Dur_comp FT_conj];

    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(resp,0.5,0.05);
    mOI = mCI; semOI = semCI;
    
    clear respRV conjV comp1V comp2V
    for an = 1:5
        respRV(an,:) = diag(all_CI_mat(:,:,an));
        conjV(an,:) = diag(all_CI_mat(:,:,an),1);
        comp1V(an,:) = diag(uni(:,:,an),1);
        comp2V(an,:) = diag(uni(:,:,an),-1);
    end
    
%     for an = 1:5
%         conjV(an,3) = all_CI_mat(1,3,an);
%         comp1V(an,3) = uni(1,3,an);
%         comp2V(an,3) = uni(3,1,an);
%     end
    
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
