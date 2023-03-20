function raster_properties

mData = evalin('base','mData'); tcolors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;

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
si = [C1_t_D C2_t_D C3_t_D C4_t_D];
%     si = [C1_i_T C2_i_T C3_i_T C4_i_T];
Rs_C = oC.Rs(:,si); Rs_A = oA.Rs(:,si); mRs_C = oC.mR(:,si); mRs_A = oA.mR(:,si);
props_C = get_props_Rs(oC.Rs(:,si),ntrials); props_A = get_props_Rs(oA.Rs(:,si),ntrials);
% pop_var_name = {'all','vals','valsT','Nvals','good_zMI','Ngood_zMI'};
pop_var_name = {'good_zMI','good_Gauss','good_MFR'};
pop_var_name = {'all'};
pop_var_name = {'vals'};
pop_var_name = {'vals'};
sel_pop_C = cell_list_op(props_C,pop_var_name); sel_pop_A = cell_list_op(props_A,pop_var_name);
% pop_var_name = {'Nvals'};
% sel_pop_CNR = cell_list_op(props_C,pop_var_name); sel_pop_ANR = cell_list_op(props_A,pop_var_name);

params = {'perc','N_Resp_Trials','zMI','rs','nan_zMI','nan_rs','HaFD','HiFD','PWs','centers','peak_locations','mean_FR','MFR'};
varT = 2;%:length(params)
[~,~,pop_C] = plotDistributions(sel_pop_C);  [~,~,pop_A] = plotDistributions(sel_pop_A);
eval(sprintf('var_CT = props_C.%s;',params{varT}));  eval(sprintf('var_AT = props_A.%s;',params{varT}));
[~,~,var_C] = plotDistributions(var_CT);  [~,~,var_A] = plotDistributions(var_AT);
var_Cv = get_vals(var_C,pop_C); var_Av = get_vals(var_A,pop_A);

for ci = 1:4
    [hall(ci),pa(ci),ks2stata(ci)] = kstest2(var_Cv{ci},var_Av{ci});
end
[hall',pa',ks2stata']

var_Cvp = get_vals(var_CT,sel_pop_C); var_Avp = get_vals(var_AT,sel_pop_A);
mean_var_C = exec_fun_on_cell_mat(var_Cvp,'nanmean'); mean_var_A = exec_fun_on_cell_mat(var_Avp,'nanmean'); 

for ci = 1:4
    [ht(ci),pt(ci),~] = ttest2(mean_var_C(:,ci),mean_var_A(:,ci));
end
[ht',pt']


%% 
    xlabels = {NaN,'Resp Fidelity','zMI'};
    incrs = [NaN 10];
    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[6 3 6.9 1],'RowsCols',[1 4],'spaceRowsCols',[0.03 -0.01],'rightUpShifts',[0.1 0.15],'widthHeightAdjustment',[10 -285]);
    MY = 8; ysp = 1; mY = 0; 
    stp = 0.25*magfac; widths = ([0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]+0.21)*magfac; gap = 0.69*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{'Cell #'});
    var_CvD = get_vals(var_CT,sel_pop_C); var_AvD = get_vals(var_AT,sel_pop_A); 
    for ci = 1:4
        distD = [var_CvD(:,ci) var_AvD(:,ci)];
        [~,~,var_Ctt] = plotDistributions(var_CvD(:,ci)); [~,~,var_Att] = plotDistributions(var_AvD(:,ci));
        [h,p,ks2stat] = kstest2(var_Ctt{1},var_Att{1});
        tcolors = {'k','r'};
        [distDo,allVals,allValsG] = plotDistributions(distD);
        minBin = min(allVals);
        maxBin = max(allVals);
        incr = incrs(varT); %maxBin =
%         hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.85 1],'color','w');
        axes(ff.h_axes(1,ci)); hold on;
    %    [ha,hb,hca] = plotDistributions(distD,'colors',tcolors,'maxY',maxBin,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
        [ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
        set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
        changePosition(gca,[0.129 0.15 -0.09 -0.13]);
        ylim([0 100]);
        if ci == 1
            put_axes_labels(gca,{xlabels{varT},[0 0 0]},{{'Neurons (%)'},[0 0 0]});
        else
            put_axes_labels(gca,{xlabels{varT},[0 0 0]},{{''},[0 0 0]});
        end
        if ci == 4
        pos = get(ff.h_axes(1,ci),'Position'); ha = axes; set(ha,'units','inches');set(ha,'Position',pos);
        changePosition(ha,[pos(1)+0.3 0 -0.5 0]);
        box off;
        end
%         [h,p,ks2stat] = kstest2(allValsG{1},allValsG{2})
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('Distribution_firing_rate'),600);
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
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'C-TG','A-TG'});
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.051 0 0]);
format_axes_b(gca);
set(ht,'FontWeight','Bold');
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% two graphs, all and cond

ff = makeFigureRowsCols(107,[10 5 2 1],'RowsCols',[1 2],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
switch varT
    case 1 % responsive cells 
        MY = 65; ysp = 2; mY = 0; titletxt = ''; ylabeltxt = {'Percent of Spatially','Tuned Cells'}; % for all cells (vals) MY = 80
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
    h(h==1) = 0;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'C1','C2','C3','C4'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(0);
make_bars_hollow(hbs(5:end));
[~,hyl] = put_axes_labels(gca,{'',[]},{ylabeltxt,[]}); set(hyl,'FontWeight','bold');
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'C-TG','A-TG'});
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
        MY = 10; ysp = 0.5; mY = 0; titletxt = 'Mutual Information'; ylabeltxt = {'Z-Score'};% for all cells (vals) MY = 70
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
        MY = 12; ysp = 2; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'% of Cells'};
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
        MY = 50; ysp = 1; mY = 0; titletxt = 'Field Width'; ylabeltxt = {'cm'};
    case 4
        MY = 0.7; ysp = 3; mY = 0; titletxt = ''; ylabeltxt = {'R-squared'};
    case 1
        MY = 1.5; ysp = 0.5; mY = 0; titletxt = 'Mutual Information'; ylabeltxt = {'Z-Score'};% for all cells (vals) MY = 70
end
stp = 0.3;magfac; widths = ([4.1 0.5 1.3 1.3 1.3 0.5 0.5 0.5]+0.18)*magfac; gap = 0.1*magfac;
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
set(gca,'xtick',xticks,'xticklabels',xticklabels);
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
set(gca,'xtick',xticks,'xticklabels',xticklabels);
put_axes_labels(gca,{[],[0 0 0]},{[],[0 0 0]});
%     changePosition(gca,[-0.03 0.03 0.11 -0.01]);
%     set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'C1','C2','C3','C4'});
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'C1','C2','C3','C4'},{[-0.02 0.01]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Pooled'},{[0 0]});

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

save_pdf(ff.hf,mData.pdf_folder,sprintf('%s_distributions_over_belt_%d',all_variables{vn},number_of_bins),600);


%% remapping across conditions

tic
ntrials = 50;
si = [C1_t_D C2_t_D C3_t_D C4_t_D];
%     si = [C1_i_T C2_i_T C3_i_T C4_i_T];
Rs_C = oC.Rs(:,si); Rs_A = oA.Rs(:,si); mRs_C = oC.mR(:,si); mRs_A = oA.mR(:,si); mRs_C1 = oC.mR1(:,si); mRs_A1 = oA.mR1(:,si);
props_C = get_props_Rs(oC.Rs(:,si),ntrials); props_A = get_props_Rs(oA.Rs(:,si),ntrials);
pop_var_name = {'vals','good_zMI'};
pop_var_name = {'vals'};
sel_pop_C = cell_list_op(props_C,pop_var_name); sel_pop_A = cell_list_op(props_A,pop_var_name);
r_sel_pop_C = cell_list_op(sel_pop_C,[],'or',1);  r_sel_pop_A = cell_list_op(sel_pop_A,[],'or',1);
remap_C = find_population_vector_corr_remap(Rs_C,mRs_C,r_sel_pop_C);
remap_A = find_population_vector_corr_remap(Rs_A,mRs_A,r_sel_pop_A);
% remap_C1 = find_population_vector_corr_remap(Rs_C,mRs_C1,r_sel_pop_C);
% remap_A1 = find_population_vector_corr_remap(Rs_A,mRs_A1,r_sel_pop_A);
toc
%%
%%%%% Correlations across conditions
selC = remap_C; selA = remap_A;
typeCorr = {'Spatial Correlation',{'Population Vector','Correlation'},'\Delta FR Score'};
FF = {'SP','PV','RR'};
ci = 3;
[within,dvn,xlabels] = make_within_table({'Cond'},3);
switch ci
    case 1
        var_C = arrayfun(@(x) mean(x{1}),selC.adj_SP_corr_diag);var_A = arrayfun(@(x) mean(x{1}),selA.adj_SP_corr_diag);
        MY = 0.35; ysp = 0.04; mY = 0; titletxt = 'Spatial Corr.'; ylabeltxt = {'A.U.'};
    case 2
        var_C = arrayfun(@(x) mean(x{1}),selC.adj_PV_corr_diag);var_A = arrayfun(@(x) mean(x{1}),selA.adj_PV_corr_diag);
        MY = 0.25; ysp = 0.05; mY = 0; titletxt = 'PV Corr.'; ylabeltxt = {'A.U.'};
    case 3
        var_C = arrayfun(@(x) nanmean(x{1}),selC.adj_RR_SP);var_A = arrayfun(@(x) nanmean(x{1}),selA.adj_RR_SP);
        MY = 0.51; ysp = 0.01; mY = 0; titletxt = 'Rate Remap'; ylabeltxt = {'A.U.'};
end
dataT = make_between_table({var_C;var_A},dvn);
ra = RMA(dataT,within);
ra.ranova

% %%%%% Correlations across conditions
% selC = remap_C1; selA = remap_A1;
% typeCorr = {'Spatial Correlation',{'Population Vector','Correlation'},'\Delta FR Score'};
% FF = {'SP','PV','RR'};
% ci = 3;
% [within,dvn,xlabels] = make_within_table({'Cond'},3);
% switch ci
%     case 1
%         var_C = arrayfun(@(x) mean(x{1}),selC.adj_SP_corr_diag);var_A = arrayfun(@(x) mean(x{1}),selA.adj_SP_corr_diag);
%         MY = 0.35; ysp = 0.04; mY = 0; titletxt = 'Spatial Corr.'; ylabeltxt = {'A.U.'};
%     case 2
%         var_C = arrayfun(@(x) mean(x{1}),selC.adj_PV_corr_diag);var_A = arrayfun(@(x) mean(x{1}),selA.adj_PV_corr_diag);
%         MY = 0.25; ysp = 0.05; mY = 0; titletxt = 'PV Corr.'; ylabeltxt = {'A.U.'};
%     case 3
%         var_C = arrayfun(@(x) nanmean(x{1}),selC.adj_RR_SP);var_A = arrayfun(@(x) nanmean(x{1}),selA.adj_RR_SP);
%         MY = 0.51; ysp = 0.01; mY = 0; titletxt = 'Rate Remap'; ylabeltxt = {'A.U.'};
% end
% dataT = make_between_table({var_C;var_A},dvn);
% ra1 = RMA(dataT,within);
% ra1.ranova
%%
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_Cond','hsd'},[1.5 1 1]);
xdata = make_xdata([3 3],[1 1.5]);   
%         [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
colors = mData.dcolors;
axes(ff.h_axes(1,1));
tcolors = mData.dcolors(1:3); tcolors = repmat(tcolors,1,2);
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp(1),'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
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
%%%% just cond
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[10 3 2.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.3],'widthHeightAdjustment',[10 -450]);
stp = 0.3*magfac; widths = ([0.5 1.3 1.3 1.3 0.5 0.5 0.5]+0.)*magfac; gap = 0.1*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};

switch ci
    case 1
        MY = 0.3; ysp = 0.04; mY = 0; titletxt = 'Spatial Corr.'; ylabeltxt = {'A.U.'};
    case 2
        MY = 0.25; ysp = 0.05; mY = 0; titletxt = 'PV Corr.'; ylabeltxt = {'A.U.'};
    case 3
        MY = 0.95; ysp = 0.001; mY = 0.8; titletxt = 'Rate Remap (ns)'; ylabeltxt = {'A.U.'};
end

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
xdata = make_xdata([3],[1 1.5]);   
%         [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
colors = mData.dcolors;
axes(ff.h_axes(1,1));
tcolors = mData.dcolors(1:3); tcolors = repmat(tcolors,1,2);
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
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
%%%% Correlations across conditions
ff = makeFigureRowsCols(108,[10 3 2.25 1.25],'RowsCols',[1 2],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
stp = 0.3*magfac; widths = ([1.2 0.5 1.3 1.3 1.3 0.5 0.5 0.5]+0.)*magfac; gap = 0.1*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};


[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_Cond','hsd'},[1.5 1 1]);
xdata = make_xdata([3 3],[1 1.5]);   
%         [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
colors = mData.dcolors;
axes(ff.h_axes(1,1));
tcolors = mData.dcolors(1:3); tcolors = repmat(tcolors,1,2);
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,[],[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp(ci),'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
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

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
xdata = make_xdata([3],[1 1.5]);   
%         [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
colors = mData.dcolors;
axes(ff.h_axes(1,2));
tcolors = mData.dcolors(1:3); tcolors = repmat(tcolors,1,2);
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp(ci),'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
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

