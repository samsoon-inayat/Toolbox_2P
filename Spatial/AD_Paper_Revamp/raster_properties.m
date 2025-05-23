function raster_properties

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;

prop_names = {'resp','N_Resp_Trials','zMI','zMINaN','HaFD','HiFD','cells_pooled'};
event_type = {'1-D','2-D','3-D','4-D','1-T','2-T','3-T','4-T'};
sic = {[C1_t_D];[C2_t_D];[C3_t_D];[C4_t_D];[C1_i_T];[C2_i_T];[C3_i_T];[C4_i_T]};
event_type = {'1-D','2-D','3-D','4-D'};
sic = {[C1_t_D];[C2_t_D];[C3_t_D];[C4_t_D]};
pni = 7;
[all_gFR_C,all_gV_C,good_zMI_C,good_zMI_MFR_C,good_zMI_MFR_Gauss_C,nan_zMI_C,all_C] = return_values_props(oC,sic,pni);
[all_gFR_A,all_gV_A,good_zMI_A,good_zMI_MFR_A,good_zMI_MFR_Gauss_A,nan_zMI_A,all_A] = return_values_props(oA,sic,pni);
%% Rs
while 1
ntrials = 50; 
sel_pop_C = good_zMI_C; sel_pop_A = good_zMI_A;
% sel_pop_C = good_zMI_MFR_C; sel_pop_A = good_zMI_MFR_A;
% sel_pop_C = good_zMI_MFR_Gauss_C; sel_pop_A = good_zMI_MFR_Gauss_A;
si = [C1_t_D C2_t_D C3_t_D C4_t_D];
props_C = get_props_Rs(oC.Rs(:,si),ntrials);
props_A = get_props_Rs(oA.Rs(:,si),ntrials);
var_C = get_vals(props_C.rs,sel_pop_C); var_A = get_vals(props_A.rs,sel_pop_A);
mean_var_C = exec_fun_on_cell_mat(var_C,'nanmean'); mean_var_A = exec_fun_on_cell_mat(var_A,'nanmean');

varC = mean_var_C;
varA = mean_var_A;
[within,dvn,xlabels] = make_within_table({'Cond'},[4]);
dataT = make_between_table({varC;varA},dvn);
ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
ra.ranova

ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -310]);
set(gcf,'color','w');    set(gcf,'Position',[10 3 1.9 1]);
MY = 8; ysp = 1; mY = 0; % responsive cells
stp = 0.3; widths = [1.7 1.3 1.3 1.3 1.3 0.5 0.5 0.5]-0.15; gap = 0.16;
adjust_axes(ff,[mY MY],stp,widths,gap,{'R-squared'});
tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_Cond','hsd'},[1.5 1 1]);
    h(h==1) = 0; xdata = make_xdata([4 4],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;

[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY 0.7]); format_axes_b(gca); xticks = xdata; 
    xticklabels = {'C1','C2','C3','C4'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(0);
    make_bars_hollow(hbs(5:end))
    ylabel('R-squared');
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Control','APP'});
    save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
break;
end
%% zMI
while 1
ntrials = 50; 
sel_pop_C = all_C; sel_pop_A = all_A;
% sel_pop_C = good_zMI_C; sel_pop_A = good_zMI_A;
% sel_pop_C = good_zMI_MFR_C; sel_pop_A = good_zMI_MFR_A;
% sel_pop_C = good_zMI_MFR_Gauss_C; sel_pop_A = good_zMI_MFR_Gauss_A;
si = [C1_t_D C2_t_D C3_t_D C4_t_D];
props_C = get_props_Rs(oC.Rs(:,si),ntrials);
props_A = get_props_Rs(oA.Rs(:,si),ntrials);
var_C = get_vals(props_C.zMI,sel_pop_C); var_A = get_vals(props_A.zMI,sel_pop_A);
mean_var_C = exec_fun_on_cell_mat(var_C,'nanmean'); mean_var_A = exec_fun_on_cell_mat(var_A,'nanmean');

varC = mean_var_C;
varA = mean_var_A;
[within,dvn,xlabels] = make_within_table({'Cond'},[4]);
dataT = make_between_table({varC;varA},dvn);
ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
ra.ranova

ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 2],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -310]);
set(gcf,'color','w');    set(gcf,'Position',[10 3 3 1]);
MY = 2; ysp = 0.2; mY = 0; % responsive cells
stp = 0.35; widths = [1.7 0.9 1.2 1.3 1.3 0.5 0.5 0.5]-0.15; gap = 0.05;
adjust_axes(ff,[mY MY],stp,widths,gap,{'R-squared'});
tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};

axes(ff.h_axes(1,1));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_Cond','hsd'},[1.5 1 1]);
    h(h==1) = 0; xdata = make_xdata([4 4],[1 1.5]);
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
    xticklabels = {'C1','C2','C3','C4'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(0);
    make_bars_hollow(hbs(5:end));
    put_axes_labels(gca,{'',[0 0 0]},{{'Mutual Information','(z-scored)'},[0 -0.1 0]});
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Control','APP'});
    
axes(ff.h_axes(1,2));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
     xdata = make_xdata([4],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
    xticklabels = {'C1','C2','C3','C4'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(0);
    make_bars_hollow(hbs(5:end));
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Pooled'});
    save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
break;
end
%% nan zMI
while 1
ntrials = 50; 
sel_pop_C = all_C; sel_pop_A = all_A;
% sel_pop_C = good_zMI_C; sel_pop_A = good_zMI_A;
% sel_pop_C = good_zMI_MFR_C; sel_pop_A = good_zMI_MFR_A;
% sel_pop_C = good_zMI_MFR_Gauss_C; sel_pop_A = good_zMI_MFR_Gauss_A;
si = [C1_t_D C2_t_D C3_t_D C4_t_D];
props_C = get_props_Rs(oC.Rs(:,si),ntrials);
props_A = get_props_Rs(oA.Rs(:,si),ntrials);
var_C = get_vals(props_C.nan_zMI,sel_pop_C); var_A = get_vals(props_A.nan_zMI,sel_pop_A);
mean_var_C = find_percent(var_C); mean_var_A = find_percent(var_A);

varC = mean_var_C;
varA = mean_var_A;
[within,dvn,xlabels] = make_within_table({'Cond'},[4]);
dataT = make_between_table({varC;varA},dvn);
ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
ra.ranova


axes(ff.h_axes(1,1));
ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -310]);
set(gcf,'color','w');    set(gcf,'Position',[10 3 1.9 1]);
MY = 60; ysp = 3; mY = 0; % responsive cells
stp = 0.35; widths = [1.7 1.3 1.3 1.3 1.3 0.5 0.5 0.5]-0.2; gap = 0.16;
adjust_axes(ff,[mY MY],stp,widths,gap,{'R-squared'});
tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_Cond','hsd'},[1.5 1 1]);
    xdata = make_xdata([4 4],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;

[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
    xticklabels = {'C1','C2','C3','C4'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(0);
    make_bars_hollow(hbs(5:end))
    ylabel({'Percent of','Silent Neurons'});
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Control','APP'});
    save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
break;
end
%% Spatially Tuned Cells
while 1
ntrials = 50; 
sel_pop_C = all_C; sel_pop_A = all_A;
sel_pop_C = good_zMI_C; sel_pop_A = good_zMI_A;
% sel_pop_C = good_zMI_MFR_C; sel_pop_A = good_zMI_MFR_A;
% sel_pop_C = good_zMI_MFR_Gauss_C; sel_pop_A = good_zMI_MFR_Gauss_A;
var_C = sel_pop_C; var_A = sel_pop_A;
mean_var_C = find_percent(var_C); mean_var_A = find_percent(var_A);

varC = mean_var_C;
varA = mean_var_A;
[within,dvn,xlabels] = make_within_table({'Cond'},[4]);
dataT = make_between_table({varC;varA},dvn);
ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
ra.ranova

ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 2],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -310]);
set(gcf,'color','w');    set(gcf,'Position',[10 3 3 1]);
MY = 25; ysp = 2; mY = 0; % responsive cells
stp = 0.37; widths = [1.6 0.9 1.1 1.3 1.3 0.5 0.5 0.5]-0.15; gap = 0.05;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};

axes(ff.h_axes(1,1));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_Cond','hsd'},[1.5 1 1]);
    h(h==1) = 0; xdata = make_xdata([4 4],[1 1.5]);
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
    xticklabels = {'C1','C2','C3','C4'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(0);
    make_bars_hollow(hbs(5:end));
    put_axes_labels(gca,{'',[0 0 0]},{{'Percent of Spatially','Tuned Cells'},[0 -2 0]});
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Control','APP'});
    
axes(ff.h_axes(1,2));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
     xdata = make_xdata([4],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
    xticklabels = {'C1','C2','C3','C4'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(0);
    make_bars_hollow(hbs(5:end));
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Pooled'});
    save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
break;
end

%% Spatially Tuned Cells Pooled
while 1
ntrials = 50; 
var_C = cell_list_op(sel_pop_C,[],'or',1); var_A = cell_list_op(sel_pop_A,[],'or',1);
mean_var_C = find_percent(var_C); mean_var_A = find_percent(var_A);


varC = mean_var_C; varA = mean_var_A;
[within,dvn,xlabels] = make_within_table({'Cond'},[2]);
dataT = make_between_table({[varC varA]},dvn);
ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
ra.ranova


ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[-350 -310]);
set(gcf,'color','w');    set(gcf,'Position',[10 3 1.1 1]);
MY = 50; ysp = 1; mY = 0; % responsive cells
stp = 0.45; widths = [1 1.3 1.3 1.3 1.3 0.5 0.5 0.5]-0.4; gap = 0.16;
adjust_axes(ff,[mY MY],stp,widths,gap,{'R-squared'});
tcolors = {'k','r'};
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
[h,p,tstat] = ttest2(varC,varA);
     xdata = make_xdata([2],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;

[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
    xticklabels = {'Control','APP','C3','C4'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(0);
    make_bars_hollow(hbs(1:end));
    put_axes_labels(gca,{'',[0 0 0]},{{'Spatially Tuned Cells','in any Condition (%)'},[0 -5 0]});
    format_axes_b(gca);
    save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
    
break;
end
