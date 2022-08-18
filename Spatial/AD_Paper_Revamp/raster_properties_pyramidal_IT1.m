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
%% general for all properties including responsivity, response fidelity, zMI, Rs
ntrials = 50;
% si = [C1_t_D C2_t_D C3_t_D C4_t_D];
si = [C1_i_T C2_i_T C3_i_T C4_i_T];
Rs_C = oC.Rs(:,si); Rs_A = oA.Rs(:,si); mRs_C = oC.mR(:,si); mRs_A = oA.mR(:,si);
props_C = get_props_Rs(oC.Rs(:,si),ntrials); props_A = get_props_Rs(oA.Rs(:,si),ntrials);
pop_var_name = {'all','vals','valsT','Nvals','good_zMI','Ngood_zMI'};
pop_var_name = {'good_zMI','good_Gauss','good_MFR'};
pop_var_name = {'all'};
pop_var_name = {'vals'};
sel_pop_C = cell_list_op(props_C,pop_var_name); sel_pop_A = cell_list_op(props_A,pop_var_name);
%%
cell_types = {'C1','C2','C3','C4'};
bins = [0 5 10 15];% bins = [0 150];
nbins = length(bins)-1;
%     sel_pop_C = bin_cells(props_C.peak_locations,bins,sel_pop_C); sel_pop_A = bin_cells(props_A.peak_locations,bins,sel_pop_A);
if nbins > 1
%     [sel_pop_C,perc_C] = bin_cells(props_C.peak_locations,bins,sel_pop_C); [sel_pop_A,perc_A] = bin_cells(props_A.peak_locations,bins,sel_pop_A);
    [sel_pop_C,perc_C] = bin_cells(props_C.centers,bins,sel_pop_C); [sel_pop_A,perc_A] = bin_cells(props_A.centers,bins,sel_pop_A);
    ind = 1;
    for ii = 1:4
        for jj = 1:nbins
            temp_tcolors{ind} = mData.colors{ii};
            ind = ind + 1;
        end
    end
end
%     


params = {'perc','N_Resp_Trials','zMI','rs','nan_zMI','nan_rs','HaFD','HiFD','PWs','centers','peak_locations','mean_FR','MFR'};
varT = 3;%:length(params)
for pii = varT
    if pii == 1
        mean_var_C = exec_fun_on_cell_mat(sel_pop_C,'percent'); mean_var_A = exec_fun_on_cell_mat(sel_pop_A,'percent');
    else
        eval(sprintf('var_C = get_vals(props_C.%s,sel_pop_C);',params{pii})); eval(sprintf('var_A = get_vals(props_A.%s,sel_pop_A);',params{pii}));
        if pii == 5 || pii == 6
            mean_var_C = exec_fun_on_cell_mat(sel_pop_C,'percent'); mean_var_A = exec_fun_on_cell_mat(sel_pop_A,'percent'); 
        else
            mean_var_C = exec_fun_on_cell_mat(var_C,'nanmean'); mean_var_A = exec_fun_on_cell_mat(var_A,'nanmean'); 
        end
    end
end
if nbins == 1
    varC = mean_var_C;
    varA = mean_var_A;
    [within,dvn,xlabels] = make_within_table({'Cond'},[4]);
    dataT = make_between_table({varC;varA},dvn);
    ra = RMA(dataT,within,{0.05,{'bonferroni'}});
    ra.ranova
else
    varC = mean_var_C;
    varA = mean_var_A;
    [within,dvn,xlabels] = make_within_table({'Cond','Bin'},[4,nbins]);
    dataT = make_between_table({varC;varA},dvn);
    ra = RMA(dataT,within,{0.05,{'bonferroni'}});
    ra.ranova
end
print_for_manuscript(ra)

%% one graph main effect
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[10 3 1.25 1.25],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
switch varT
    case 13
        MY = 3; ysp = 0.5; mY = 0; titletxt = ''; ylabeltxt = {'Firing Rate (A.U.)'};
end
stp = 0.36*magfac; widths = ([1 1.3 1.3 1.3 1.3 0.5 0.5 0.5]-0.51)*magfac; gap = 0.16*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = {'k','r'};%{colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group','hsd'},[1.5 1 1]);
    xdata = make_xdata([2],[1 1.5]);   
%         combs = [[1:2:12]' [2:2:12]']; p = ra.MC.bonferroni.Group_by_Cond{1:2:12,6}; h = p<0.05;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'Control','APP'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
make_bars_hollow(hbs(1:end));
put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
format_axes_b(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% one graph
magfac = mData.magfac; colors = mData.colors;
ff = makeFigureRowsCols(108,[10 5 1.5 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
switch varT
    case 1 % responsive cells 
        MY = 60; ysp = 3; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'Percent of Cells'};
    case 2
        MY = 70; ysp = 3; mY = 0; titletxt = 'Response Fidelity'; ylabeltxt = {'Percent of Trials'};
    case 4
        MY = 0.7; ysp = 3; mY = 0; titletxt = 'Goodness-of-Fit (1D Gauss)'; ylabeltxt = {'R-squared'};
    case 3
        MY = 7; ysp = 0.3; mY = 0; titletxt = ''; ylabeltxt = {'Mutual Information','(z-score)'};
    case 5
        MY = 60; ysp = 4; mY = 0; titletxt = ''; ylabeltxt = {'Percent of Silent','Neurons'};% for all cells (vals) MY = 70
    case 7
        MY = 1.7; ysp = 0.05; mY = 1; titletxt = 'Hausdorff Frac. Dim'; ylabeltxt = {'A.U.'};
    case 10
        MY = 70; ysp = 1; mY = 0; titletxt = ''; ylabeltxt = {'Place Field Center','Location (cm)'};
    case 11
        MY = 80; ysp = 1; mY = 0; titletxt = 'Peak Locations'; ylabeltxt = {'cm'};
end
stp = 0.3*magfac; widths = ([1.2 1.3 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.16*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_Cond','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([4 4],[1 1.5]);   
%         combs = [[1:2:12]' [2:2:12]']; p = ra.MC.bonferroni.Group_by_Cond{1:2:12,6}; h = p<0.05;
%     h(h==1) = 0;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'C1','C2','C3','C4'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(0);
make_bars_hollow(hbs(5:end));
put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Control','APP'});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.051 0 0]);
format_axes_b(gca);
set(ht,'FontWeight','Bold');
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% two graphs, all and cond

ff = makeFigureRowsCols(107,[10 5 2 1],'RowsCols',[1 2],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
switch varT
    case 1 % responsive cells 
        MY = 20; ysp = 2; mY = 0; titletxt = ''; ylabeltxt = {'Percent of Spatially','Tuned Cells'}; % for all cells (vals) MY = 80
    case 2
        MY = 90; ysp = 4; mY = 0; titletxt = 'Response Fidelity'; ylabeltxt = {'Percent of Trials'};% for all cells (vals) MY = 70
    case 3
        MY = 0.1; ysp = 0.01; mY = -0.15; titletxt = 'Mutual Information'; ylabeltxt = {'Mutual Information','Z-Score'};
        MY = 3.5; ysp = 0.3; mY = 0; titletxt = 'Mutual Information'; ylabeltxt = {'(z-score)'};
    case 9
        MY = 20; ysp = 1; mY = -0.15; titletxt = 'Tuning Width'; ylabeltxt = {'cm'};
end
stp = 0.225*magfac; widths = ([1.2 0.6 1.3 1.3 1.3 0.5 0.5 0.5]-0.05)*magfac; gap = 0.067*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
axes(ff.h_axes(1,1));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_Cond','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([4 4],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
%     h(h==1) = 0;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'C1','C2','C3','C4'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(0);
make_bars_hollow(hbs(5:end));
[~,hyl] = put_axes_labels(gca,{'',[]},{ylabeltxt,[]}); set(hyl,'FontWeight','bold');
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Control','APP'});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.051 0 0]);
set(ht,'FontWeight','Bold');
format_axes_b(gca);

axes(ff.h_axes(1,2));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([4],[1 1.5]);   
%     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
%     h(h==1) = 0;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'C1','C2','C3','C4'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(0);
%     make_bars_transparent(hbs,0.5);
%     hatch(hbs,0,'w','-',2,0.25); %hatch(obj,angle,color,style,step,width)
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Pooled'});
format_axes_b(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% pooled bin effect only
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[5 5 1.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.39],'widthHeightAdjustment',[10 -520]);
switch varT
    case 1 % responsive cells 
        MY = 12; ysp = 2; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'% of Cells'};
    case 2
        MY = 80; ysp = 3; mY = 0; titletxt = 'Response Fidelity'; ylabeltxt = {'% of Trials'};
    case 4
        MY = 1; ysp = 0.07; mY = 0; titletxt = 'Goodness-of-Fit'; ylabeltxt = {'R-squared'};
    case 3
        MY = 4.5; ysp = 0.5; mY = 0; titletxt = 'Mutual Information'; ylabeltxt = {'Z-Score'};% for all cells (vals) MY = 70
end
stp = 0.3;magfac; widths = ([0.4 0.5 1.3 1.3 0.5 0.5 0.5]+0.08)*magfac; gap = 0.1*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});


[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Bin','bonferroni'},[1.5 1 1]);
xdata = make_xdata([3],[1 1.5]);   
%     mVar = ra.est_marginal_means.Mean; semVar = ra.est_marginal_means.Formula_StdErr;
tcolors = mData.dcolors;
axes(ff.h_axes(1,1))
[hbs,maxY]  = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.01);
make_bars_hollow(hbs(13:end));
set_axes_limits(gca,[0.25 xdata(end)+.75],[mY MY]); format_axes_b(gca); xticks = xdata; 
set(gca,'xtick',xticks,'xticklabels',{'B1','B2','B3'});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[-0.15 -0.051 0.6 0]); set(ht,'FontWeight','Bold');
put_axes_labels(gca,{[],[0 0 0]},{ylabeltxt,[0 0 0]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Pooled'},{[0 0]});
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graph'),600);
%% pooled bin and cond effect only
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[5 5 1.4 1],'RowsCols',[1 2],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.39],'widthHeightAdjustment',[10 -520]);
switch varT
    case 1 % responsive cells 
        MY = 5; ysp = 1; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'% of Cells'};
    case 2
        MY = 80; ysp = 3; mY = 0; titletxt = 'Response Fidelity'; ylabeltxt = {'% of Trials'};
    case 4
        MY = 1; ysp = 0.07; mY = 0; titletxt = 'Goodness-of-Fit'; ylabeltxt = {'R-squared'};
    case 3
        MY = 8; ysp = 1; mY = 0; titletxt = 'Mutual Information'; ylabeltxt = {'Z-Score'};% for all cells (vals) MY = 70
end
stp = 0.25*magfac; widths = ([0.4 0.5 1.3 1.3 0.5 0.5 0.5]+0.08)*magfac; gap = 0.051*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});


[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Bin','bonferroni'},[1.5 1 1]);
xdata = make_xdata([3],[1 1.5]);   
%     mVar = ra.est_marginal_means.Mean; semVar = ra.est_marginal_means.Formula_StdErr;
tcolors = mData.dcolors;
axes(ff.h_axes(1,1))
[hbs,maxY]  = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.01);
make_bars_hollow(hbs(13:end));
set_axes_limits(gca,[0.25 xdata(end)+.75],[mY MY]); format_axes_b(gca); xticks = xdata; 
set(gca,'xtick',xticks,'xticklabels',{'B1','B2','B3'});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.051 0.6 0]); set(ht,'FontWeight','Bold');
put_axes_labels(gca,{[],[0 0 0]},{ylabeltxt,[0 0 0]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Pooled'},{[0 0]});

% set(ff.h_axes(1,3),'Visible','Off');
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','bonferroni'},[1.5 1 1]);
xdata = make_xdata([4],[1 1.5]);   
%     mVar = ra.est_marginal_means.Mean; semVar = ra.est_marginal_means.Formula_StdErr;
tcolors = colors;
axes(ff.h_axes(1,2))
[hbs,maxY]  = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.01);
make_bars_hollow(hbs(13:end));
set_axes_limits(gca,[0.25 xdata(end)+.75],[mY MY]); format_axes_b(gca); xticks = xdata; 
set(gca,'xtick',xticks,'xticklabels',{'C1','C2','C3','C4'});
put_axes_labels(gca,{[],[0 0 0]},{[],[0 0 0]});
%     changePosition(gca,[-0.03 0.03 0.11 -0.01]);
%     set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'C1','C2','C3','C4'});
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'C1','C2','C3','C4'},{[-0.02 0.01]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Pooled'},{[0 0]});


save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graph'),600);

%% pooled Cond_by_bin effect only
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[5 5 2.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.39],'widthHeightAdjustment',[10 -520]);
switch varT
    case 1 % responsive cells 
        MY = 5; ysp = 1; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'% of Cells'};
    case 2
        MY = 90; ysp = 10; mY = 0; titletxt = 'Response Fidelity'; ylabeltxt = {'% of Trials'};
    case 4
        MY = 1; ysp = 0.07; mY = 0; titletxt = 'Goodness-of-Fit'; ylabeltxt = {'R-squared'};
    case 3
        MY = 8; ysp = 1; mY = 0; titletxt = 'Mutual Information'; ylabeltxt = {'Z-Score'};% for all cells (vals) MY = 70
end
stp = 0.25*magfac; widths = ([1.75 0.5 1.3 1.3 0.5 0.5 0.5]+0.08)*magfac; gap = 0.051*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});


[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_Bin','bonferroni'},[1.5 1 1]);
xdata = make_xdata([3 3 3 3],[1 1.5]);   
%     mVar = ra.est_marginal_means.Mean; semVar = ra.est_marginal_means.Formula_StdErr;
tcolors = temp_tcolors;
axes(ff.h_axes(1,1))
[hbs,maxY]  = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.01);
make_bars_hollow(hbs(13:end));
set_axes_limits(gca,[0.25 xdata(end)+.75],[mY MY]); format_axes_b(gca); xticks = xdata; 
set(gca,'xtick',xticks,'xticklabels',{'B1','B2','B3'});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'C1','C2','C3','C4'},{[-0.01 0.02]});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.051 0 0]); set(ht,'FontWeight','Bold');
put_axes_labels(gca,{[],[0 0 0]},{ylabeltxt,[0 0 0]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,12,{'Pooled'},{[-0.12 0]});
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graph'),600);

%% cond effect only
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[5 5 1.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.39],'widthHeightAdjustment',[10 -520]);
switch varT
    case 1 % responsive cells 
        MY = 12; ysp = 2; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'% of Cells'};
    case 2
        MY = 80; ysp = 3; mY = 0; titletxt = 'Response Fidelity'; ylabeltxt = {'% of Trials'};
    case 4
        MY = 1; ysp = 0.07; mY = 0; titletxt = 'Goodness-of-Fit'; ylabeltxt = {'R-squared'};
    case 3
        MY = 8; ysp = 1; mY = 0; titletxt = 'Mutual Information'; ylabeltxt = {'Z-Score'};% for all cells (vals) MY = 70
    case 9
        MY = 25; ysp = 3; mY = 0; titletxt = 'Field Width'; ylabeltxt = {'cm'};% for all cells (vals) MY = 70
end
stp = 0.3*magfac; widths = ([0.5 0.5 1.3 1.3 0.5 0.5 0.5]+0.08)*magfac; gap = 0.051*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});

% set(ff.h_axes(1,3),'Visible','Off');
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','bonferroni'},[1.5 1 1]);
xdata = make_xdata([4],[1 1.5]);   
%     mVar = ra.est_marginal_means.Mean; semVar = ra.est_marginal_means.Formula_StdErr;
tcolors = colors;
axes(ff.h_axes(1,1))
[hbs,maxY]  = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.01);
make_bars_hollow(hbs(13:end));
set_axes_limits(gca,[0.25 xdata(end)+.75],[mY MY]); format_axes_b(gca); xticks = xdata; 
set(gca,'xtick',xticks,'xticklabels',{'C1','C2','C3','C4'});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.051 0 0]); set(ht,'FontWeight','Bold');
put_axes_labels(gca,{[],[0 0 0]},{ylabeltxt,[0 0 0]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Pooled'},{[0 0]});


save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graph'),600);

%% average distributions w.r.t centers for the two groups all bars and main effect of bins and Cond
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[5 5 6.9 1],'RowsCols',[1 3],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.39],'widthHeightAdjustment',[10 -520]);
switch varT
    case 1 % responsive cells 
        MY = 12; ysp = 2; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'% of Cells'};
    case 2
        MY = 90; ysp = 3; mY = 0; titletxt = 'Response Fidelity'; ylabeltxt = {'Percent of Trials'};
    case 4
        MY = 0.7; ysp = 3; mY = 0; titletxt = ''; ylabeltxt = {'R-squared'};
    case 3
        MY = 10; ysp = 0.5; mY = 0; titletxt = 'Mutual Information'; ylabeltxt = {'Z-Score'};% for all cells (vals) MY = 70
end
stp = 0.3;magfac; widths = ([3.3 0.4 0.5 1.3 1.3 0.5 0.5 0.5]+0.08)*magfac; gap = 0.051*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_Cond_Bin','bonferroni'},[1.5 1 1]);
xdata = make_xdata([12 12],[1 1.5]);   
%     mVar = ra.est_marginal_means.Mean; semVar = ra.est_marginal_means.Formula_StdErr;
tcolors = [temp_tcolors temp_tcolors];
%     combs = ra.mcs.combs; p = ra.mcs.p; h = p<0.00005;
%     xdata = [1:(length(mVar)/2) ((length(mVar)/2)+1+(1:(length(mVar)/2)))];
colors = mData.colors;
%     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 6.9 2],'color','w');
axes(ff.h_axes(1,1))
[hbs,maxY]  = plotBarsWithSigLines(mVar,semVar,[],[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.3);
make_bars_hollow(hbs(13:end));
set_axes_limits(gca,[0.25 xdata(end)+.75],[mY MY]); format_axes_b(gca); xticks = xdata; 
set(gca,'xtick',xticks,'xticklabels',{'B1','B2','B3'});
put_axes_labels(gca,{[],[0 0 0]},{ylabeltxt,[0 0 0]});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.051 0 0]); set(ht,'FontWeight','Bold');
%     changePosition(gca,[-0.03 0.03 0.11 -0.01]);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'C1','C2','C3','C4','C1','C2','C3','C4'},{[-0.01 0.02]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,12,{'Control','APP'},{[-0.12 0]});

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Bin','bonferroni'},[1.5 1 1]);
xdata = make_xdata([3],[1 1.5]);   
%     mVar = ra.est_marginal_means.Mean; semVar = ra.est_marginal_means.Formula_StdErr;
tcolors = mData.dcolors;
axes(ff.h_axes(1,2))
[hbs,maxY]  = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.01);
make_bars_hollow(hbs(13:end));
set_axes_limits(gca,[0.25 xdata(end)+.75],[mY MY]); format_axes_b(gca); xticks = xdata; 
set(gca,'xtick',xticks,'xticklabels',{'B1','B2','B3'});
put_axes_labels(gca,{[],[0 0 0]},{[],[0 0 0]});
%     changePosition(gca,[-0.03 0.03 0.11 -0.01]);
%     set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'C1','C2','C3','C4'});
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'C1','C2','C3','C4'},{[-0.02 0.01]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Pooled'},{[0 0]});

% set(ff.h_axes(1,3),'Visible','Off');
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','bonferroni'},[1.5 1 1]);
xdata = make_xdata([4],[1 1.5]);   
%     mVar = ra.est_marginal_means.Mean; semVar = ra.est_marginal_means.Formula_StdErr;
tcolors = colors;
axes(ff.h_axes(1,3))
[hbs,maxY]  = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.01);
make_bars_hollow(hbs(13:end));
set_axes_limits(gca,[0.25 xdata(end)+.75],[mY MY]); format_axes_b(gca); xticks = xdata; 
set(gca,'xtick',xticks,'xticklabels',{'C1','C2','C3','C4'});
put_axes_labels(gca,{[],[0 0 0]},{[],[0 0 0]});
%     changePosition(gca,[-0.03 0.03 0.11 -0.01]);
%     set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'C1','C2','C3','C4'});
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'C1','C2','C3','C4'},{[-0.02 0.01]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Pooled'},{[0 0]});


save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graph'),600);

%% average distributions w.r.t centers for the two groups all bars and main effect of bins
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[5 5 6.9 1],'RowsCols',[1 3],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.39],'widthHeightAdjustment',[10 -520]);
switch varT
    case 5 % responsive cells 
        MY = 5; ysp = 7; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'Percent of Cells'};
    case 9
        MY = 15; ysp = 1; mY = 0; titletxt = 'Field Width'; ylabeltxt = {'s'};
    case 4
        MY = 1; ysp = 3; mY = 0; titletxt = 'Goodness-of-Fit'; ylabeltxt = {'R-squared'};
    case 1
        MY = 1.5; ysp = 0.5; mY = 0; titletxt = 'Mutual Information'; ylabeltxt = {'Z-Score'};% for all cells (vals) MY = 70
end
stp = 0.3;magfac; widths = ([3.3 0.4 0.5 1.3 1.3 0.5 0.5 0.5]+0.08)*magfac; gap = 0.051*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_Cond_Bin','bonferroni'},[1.5 1 1]);
xdata = make_xdata([12 12],[1 1.5]);   
%     mVar = ra.est_marginal_means.Mean; semVar = ra.est_marginal_means.Formula_StdErr;
tcolors = [temp_tcolors temp_tcolors];
%     combs = ra.mcs.combs; p = ra.mcs.p; h = p<0.00005;
%     xdata = [1:(length(mVar)/2) ((length(mVar)/2)+1+(1:(length(mVar)/2)))];
colors = mData.colors;
%     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 6.9 2],'color','w');
axes(ff.h_axes(1,1))
[hbs,maxY]  = plotBarsWithSigLines(mVar,semVar,[],[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.3);
make_bars_hollow(hbs(13:end));
set_axes_limits(gca,[0.25 xdata(end)+.75],[mY MY]); format_axes_b(gca); xticks = xdata; 
set(gca,'xtick',xticks,'xticklabels',{'B1','B2','B3'});
put_axes_labels(gca,{[],[0 0 0]},{ylabeltxt,[0 0 0]});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.051 0 0]); set(ht,'FontWeight','Bold');
%     changePosition(gca,[-0.03 0.03 0.11 -0.01]);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'C1','C2','C3','C4','C1','C2','C3','C4'},{[-0.01 0.02]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,12,{'Control','APP'},{[-0.12 0]});

set(ff.h_axes(1,2),'Visible','Off');

set(ff.h_axes(1,3),'Visible','Off');
% [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
% xdata = make_xdata([4],[1 1.5]);   
% %     mVar = ra.est_marginal_means.Mean; semVar = ra.est_marginal_means.Formula_StdErr;
% tcolors = colors;
% axes(ff.h_axes(1,3))
% [hbs,maxY]  = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
%     'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
%     'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.01);
% make_bars_hollow(hbs(13:end));
% set_axes_limits(gca,[0.25 xdata(end)+.75],[mY MY]); format_axes_b(gca); xticks = xdata; 
% set(gca,'xtick',xticks,'xticklabels',{'C1','C2','C3','C4'});
% put_axes_labels(gca,{[],[0 0 0]},{[],[0 0 0]});
% %     changePosition(gca,[-0.03 0.03 0.11 -0.01]);
% %     set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'C1','C2','C3','C4'});
% % set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'C1','C2','C3','C4'},{[-0.02 0.01]});
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Pooled'},{[0 0]});
% 

save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graph'),600);


%% remapping across conditions

tic
ntrials = 50;
% si = [C1_t_D C2_t_D C3_t_D C4_t_D];
si = [C1_i_T C2_i_T C3_i_T C4_i_T];
Rs_C = oC.Rs(:,si); Rs_A = oA.Rs(:,si); mRs_C = oC.mR(:,si); mRs_A = oA.mR(:,si);
props_C = get_props_Rs(oC.Rs(:,si),ntrials); props_A = get_props_Rs(oA.Rs(:,si),ntrials);
pop_var_name = {'vals'};
sel_pop_C = cell_list_op(props_C,pop_var_name); sel_pop_A = cell_list_op(props_A,pop_var_name);
r_sel_pop_C = cell_list_op(sel_pop_C,[],'or',1);  r_sel_pop_A = cell_list_op(sel_pop_A,[],'or',1);
remap_C = find_population_vector_corr_remap(Rs_C,mRs_C,r_sel_pop_C);
remap_A = find_population_vector_corr_remap(Rs_A,mRs_A,r_sel_pop_A);
toc
%%
%%%%% Correlations across conditions
magfac = mData.magfac
selC = remap_C; selA = remap_A;
typeCorr = {'Spatial Correlation',{'Population Vector','Correlation'},'\Delta FR Score'};
FF = {'SP','PV','RR'};
ysp = [0.05 0.05 0.1];
ci = 2;
[within,dvn,xlabels] = make_within_table({'Cond'},3);
switch ci
    case 1
        var_C = arrayfun(@(x) mean(x{1}),selC.adj_SP_corr_diag);var_A = arrayfun(@(x) mean(x{1}),selA.adj_SP_corr_diag);
        MY = 0.172; ysp = 0.03; mY = 0; titletxt = 'Temporal Correlation'; ylabeltxt = {'A.U.'};
    case 2
        var_C = arrayfun(@(x) mean(x{1}),selC.adj_PV_corr_diag);var_A = arrayfun(@(x) mean(x{1}),selA.adj_PV_corr_diag);
        MY = 0.25; ysp = 0.02; mY = 0; titletxt = 'PV Correlation'; ylabeltxt = {'A.U.'};
    case 3
        var_C = arrayfun(@(x) nanmean(x{1}),selC.adj_RR_SP);var_A = arrayfun(@(x) nanmean(x{1}),selA.adj_RR_SP);
        MY = 0.97; ysp = 0.01; mY = 0.85; titletxt = 'Rate Remap'; ylabeltxt = {'A.U.'};
end
dataT = make_between_table({var_C;var_A},dvn);
ra = RMA(dataT,within);
ra.ranova
% just cond
ff = makeFigureRowsCols(108,[10 5 1.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.3],'widthHeightAdjustment',[10 -450]);
stp = 0.3*magfac; widths = ([0.5 1.3 1.3 1.3 0.5 0.5 0.5]+0.)*magfac; gap = 0.1*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};


[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
xdata = make_xdata([3],[1 1.5]);   
%         [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
colors = mData.dcolors;
axes(ff.h_axes(1,1));
tcolors = mData.dcolors(1:3); tcolors = repmat(tcolors,1,2);
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.01);
set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[mY MY],'FontSize',6,'FontWeight','Normal','TickDir','out');
xticks = xdata(1:end)+0; xticklabels = {'C12','C23','C34','C12','C23','C34'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
make_bars_hollow(hbs(4:end));
xtickangle(30);
put_axes_labels(gca,{[],[0 0 0]},{ylabeltxt,[0 0 0]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Pooled'});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.01 0 0]); set(ht,'FontWeight','Bold');
format_axes_b(gca);

save_pdf(ff.hf,mData.pdf_folder,sprintf('%s_correlation',FF{ci}),600);
%%
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_Cond','hsd'},[1.5 1 1]);
xdata = make_xdata([3 3],[1 1.5]);   
%         [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
colors = mData.dcolors;
axes(ff.h_axes(1,1));
tcolors = mData.dcolors(1:3); tcolors = repmat(tcolors,1,2);
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp(ci),'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
xticks = xdata(1:end)+0; xticklabels = {'C12','C23','C34','C12','C23','C34'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
make_bars_hollow(hbs(4:end));
for ii = 4:6
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
end
xtickangle(30);
put_axes_labels(gca,{[],[0 0 0]},{typeCorr{ci},[0 0 0]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Control','APP'});
format_axes_b(gca);
save_pdf(ff.hf,mData.pdf_folder,sprintf('%s_correlation',FF{ci}),600);
%%
%%%% Correlations across conditions
ff = makeFigureRowsCols(108,[10 7 2.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.3],'widthHeightAdjustment',[10 -450]);
stp = 0.3*magfac; widths = ([1.2 0.5 1.3 1.3 1.3 0.5 0.5 0.5]+0.)*magfac; gap = 0.1*magfac;
MY = 0.22; ysp = 0.01; mY = 0; titletxt = 'PV Correlation'; ylabeltxt = {'A.U.'};
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};


[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_Cond','hsd'},[1.5 1 1]);
xdata = make_xdata([3 3],[1 1.5]);   
%         [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
colors = mData.dcolors;
axes(ff.h_axes(1,1));
tcolors = mData.dcolors(1:3); tcolors = repmat(tcolors,1,2);
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0);
set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[mY MY],'FontSize',6,'FontWeight','Normal','TickDir','out');
xticks = xdata(1:end)+0; xticklabels = {'C12','C23','C34','C12','C23','C34'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
make_bars_hollow(hbs(4:end));
for ii = 4:6
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
end
xtickangle(30);
put_axes_labels(gca,{[],[0 0 0]},{ylabeltxt,[0 0 0]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Control','APP'});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 0 0 0]); set(ht,'FontWeight','Bold');
format_axes_b(gca);
save_pdf(ff.hf,mData.pdf_folder,sprintf('%s_correlation',FF{ci}),600);

%%
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
xdata = make_xdata([3],[1 1.5]);   
%         [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
colors = mData.dcolors;
axes(ff.h_axes(1,2));
tcolors = mData.dcolors(1:3); tcolors = repmat(tcolors,1,2);
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[mY MY],'FontSize',6,'FontWeight','Normal','TickDir','out');
xticks = xdata(1:end)+0; xticklabels = {'C12','C23','C34','C12','C23','C34'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
make_bars_hollow(hbs(4:end));
xtickangle(30);
put_axes_labels(gca,{[],[0 0 0]},{[],[0 0 0]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Pooled'});
format_axes_b(gca);

save_pdf(ff.hf,mData.pdf_folder,sprintf('%s_correlation',FF{ci}),600);

