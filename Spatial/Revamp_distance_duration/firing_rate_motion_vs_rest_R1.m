function firing_rate_motion_vs_rest

Rs_C = o.Rs;% get_rasters_data(ei_C,selContexts,rasterNames);
% typeP = {'all','vals'
typeP  = 'all';
thr = -1;
%%
ntrials = 50;
si = [Ab_t_T Ab_i_T Ar_t_D Ar_i_T ArL_t_D ArL_i_T Ars_t_D Ars_i_T Abs_t_T Abs_i_T];
Rs_C = o.Rs(:,si); mRs_C = o.mR(:,si); 
props_C = get_props_Rs(o.Rs(:,si),ntrials);
% pop_var_name = {'all','vals','valsT','Nvals','good_zMI','Ngood_zMI'};
pop_var_name = {'good_zMI','good_Gauss','good_MFR'};
pop_var_name = {'all'};
% pop_var_name = {'vals'};
% pop_var_name = {'vals','good_zMI'};
sel_pop_C = cell_list_op(props_C,pop_var_name);

disp('Done')
%%
out_C = get_spike_rate(ei,thr,sel_pop_C);


   
% rest vs motion FR

    tcolors = {'k','r'};

    data_C = [out_C.m_sp_animal_level_rest' out_C.m_sp_animal_level_motion'];
    [within,dvn,xlabels] = make_within_table({'St'},[2]);
    dataT = make_between_table({data_C},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd'}});
%     ra.ranova
    print_for_manuscript(ra);

   magfac = mData.magfac;
ff = makeFigureRowsCols(108,[10 3 1.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.31 0.25],'widthHeightAdjustment',[10 -350]);
MY = 0.15; ysp = 0.03; mY = 0; titletxt = ''; ylabeltxt = {'Avg. FR','(A.U.)'};
stp = 0.45*magfac; widths = ([0.5 0.4 0.4 1.3 1.3 0.5 0.5 0.5])*magfac; gap = 0.051*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
colors = mData.colors;
tcolors = {'b','m'};%repmat({colors{1};colors{2};colors{3};colors{4}},4);
axes(ff.h_axes(1,1));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'St','hsd'},[1.5 1 1]);
xdata = make_xdata([2],[1 1.5]); 
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca); xticks = xdata; 
xticklabels = {'Rest','Motion'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30);
make_bars_hollow(hbs(9:end));
put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Air-On','Air-Off','Air-On','Air-Off'},{[-0.01 0.03]});
% ht = set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,8,{'C-TG','A-TG'},{[-0.14 -0]}); 
format_axes(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

    


%% spike rate
out_C = get_spike_rate_ext_dist_dur(ei,thr,sel_pop_C,selContexts(si));
% out_C_Ca = get_spike_rate_ext_Ca(ei_C,thr,sel_pop_C);
disp('Done');
%%
sel_c  = [9 10];
varC = exec_fun_on_cell_mat(out_C,'mean');
[within,dvn,xlabels] = make_within_table({'Ph'},[2]);
dataT = make_between_table({varC(:,sel_c)},dvn);
ra = RMA(dataT,within,{0.05,{'bonferroni'}});
ra.ranova
print_for_manuscript(ra)

%%
sel_c  = [1 2 3 4];
sel_c  = [3 4 5 6 7 8];
sel_c  = [7 8 9 10];
% sel_c  = [1 2 9 10];
varC = exec_fun_on_cell_mat(out_C,'mean');
[within,dvn,xlabels] = make_within_table({'Ph','Cond'},[2,length(sel_c)/2]);
dataT = make_between_table({varC(:,sel_c)},dvn);
ra = RMA(dataT,within,{0.05,{'bonferroni'}});
% ra.ranova
print_for_manuscript(ra)

%%
varC = exec_fun_on_cell_mat(out_C,'mean');
[within,dvn,xlabels] = make_within_table({'Ph','Cond'},[2,5]);
dataT = make_between_table({varC(:,:)},dvn);
ra = RMA(dataT,within,{0.05,{'bonferroni'}});
ra.ranova
print_for_manuscript(ra)


%%
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[10 3 3.6 1],'RowsCols',[1 3],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.31 0.40],'widthHeightAdjustment',[10 -450]);
MY = 0.2; ysp = 0.02; mY = 0; titletxt = ''; ylabeltxt = {'Average Firing','Rate (A.U.)'};
stp = 0.45*magfac; widths = ([2.25 0.4 0.4 1.3 1.3 0.5 0.5 0.5])*magfac; gap = 0.051*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
colors = mData.colors;
tcolors = repmat({colors{1};colors{2}},5);
axes(ff.h_axes(1,1));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Ph_by_Cond','bonferroni'},[1.5 1 1]);
xdata = make_xdata([2 2 2 2 2],[1 1.5]); 
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'C2','C2','C3','C3','C4','C4','C5','C5','C7','C7'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(0);
make_bars_hollow(hbs(2:2:end));
put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Air-On','Air-Off','Air-On','Air-Off'},{[-0.01 0.03]});
% ht = set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,8,{'C-TG','A-TG'},{[-0.14 -0]}); 
for ii = 1:length(ht) 
  set(ht(ii),'FontWeight','Bold');
end
format_axes_b(gca);

tcolors = {'k','r'};%{colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
axes(ff.h_axes(1,2));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([2],[1 2]);   
%         combs = [[1:2:12]' [2:2:12]']; p = ra.MC.bonferroni.Group_by_Cond{1:2:12,6}; h = p<0.05;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.45,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'C-TG','A-TG'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30);
make_bars_hollow(hbs(1:end));
ht = set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[-0.05 -0]}); 
format_axes_b(gca);

tcolors = {mData.colors{7};mData.colors{6}};
axes(ff.h_axes(1,3));

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Ph','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([2],[1 2]);   
%         combs = [[1:2:12]' [2:2:12]']; p = ra.MC.bonferroni.Group_by_Cond{1:2:12,6}; h = p<0.05;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.45,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'Air-On','Air-Off'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30);
format_axes_b(gca);
% set(gca,'ytick',[]);
ht = set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[-0.04 -0]}); 
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%%
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[10 3 1.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.31 0.40],'widthHeightAdjustment',[10 -450]);
MY = 0.15; ysp = 0.025; mY = 0; titletxt = ''; ylabeltxt = {''};

stp = 0.2*magfac; widths = ([1 1.3 1.3 1.3 1.3 0.5 0.5 0.5]-0.61)*magfac; gap = 0.16*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = {'k','r'};%{colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([2],[1 2]);   
%         combs = [[1:2:12]' [2:2:12]']; p = ra.MC.bonferroni.Group_by_Cond{1:2:12,6}; h = p<0.05;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'C-TG','A-TG'};set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[]); xtickangle(30);
make_bars_hollow(hbs(1:end));
put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
format_axes_b(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
%% one graph main effect
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[10 3 1.25 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.21 0.40],'widthHeightAdjustment',[10 -450]);

        MY = 15; ysp = 3; mY = 0; titletxt = ''; ylabeltxt = {''};

stp = 0.2*magfac; widths = ([1 1.3 1.3 1.3 1.3 0.5 0.5 0.5]-0.61)*magfac; gap = 0.16*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
tcolors = {mData.colors{7};mData.colors{6}};

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Ph','bonferroni'},[1.5 1 1]);
    xdata = make_xdata([2],[1 2]);   
%         combs = [[1:2:12]' [2:2:12]']; p = ra.MC.bonferroni.Group_by_Cond{1:2:12,6}; h = p<0.05;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'Air-On','Air-Off'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30);
% make_bars_hollow(hbs(1:end));
put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
format_axes_b(gca);
set(gca,'ytick',[]);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
%% transients graphs trials
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[4 5 6.9 1.45],'RowsCols',[1 8],'spaceRowsCols',[0.03 -0.01],'rightUpShifts',[0.1 0.3],'widthHeightAdjustment',[10 -585]);
MY = 10; ysp = 1; mY = 0; 
stp = 0.4*magfac; widths = ([0.5 0.5 0.5 0.5 0.25 0.25 0.25 0.25]+0.23)*magfac; gap = 0.175*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{'Cell #'});
shift_axes(ff,5:8,0.35,gap);
var_CvD = out_CT(:,1:4); var_AvD = out_AT(:,1:4);
var_CvD = out_CT(:,5:8); var_AvD = out_AT(:,5:8);
for ci = 1:4
    distD = [var_CvD(:,ci) var_AvD(:,ci)];
%     [~,~,var_Ctt] = plotDistributions(var_CvD(:,ci)); [~,~,var_Att] = plotDistributions(var_AvD(:,ci));
%     [h,p,ks2stat] = kstest2(var_Ctt{1},var_Att{1});
    tcolors = {'k','r'};
    [distDo,allVals,allValsG] = plotDistributions(distD);
    minBin = min(allVals);
    maxBin = max(allVals);
    incr = 0.1;
%         hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.85 1],'color','w');
    axes(ff.h_axes(1,ci)); hold on;
    [ha,hb,hca] = plotDistributions(allValsG,'colors',tcolors,'maxY',maxBin,'min',minBin,'incr',incr,'max',maxBin,'do_mean','No');
%     [ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'do_mean','Yes');
    set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
%     changePosition(gca,[0.129 0.15 -0.09 -0.13]);
    ylim([0 100]); xlim([minBin maxBin]);
    if ci == 1
        put_axes_labels(ha,{{'Avg. # of Trans.','per min'},[0 0 0]},{{'Cells (%)'},[0 0 0]});
    else
%         put_axes_labels(ha,{'Firing Rate (A.U.)',[0 0 0]},{{'Neurons (%)'},[0 0 0]});
        put_axes_labels(ha,{{'Avg. # of Trans.','per min'},[0 0 0]},{{''},[0 0 0]});
    end
    format_axes_b(ha);
    [h,p,ks2stat] = kstest2(allValsG{1},allValsG{2});
    ht = set_axes_top_text_no_line(gcf,ha,'KS-Test',[0.0 -0.01 0 0]);set(ht,'FontSize',7);
    titletxt = sprintf('%s',getNumberOfAsterisks(p));
    ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.061 -0.01 0 0]);set(ht,'FontSize',9);
    titletxt = sprintf('C%d',ci);
    ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.031 0.1 0 0]);set(ht,'FontSize',8);
    titletxt = sprintf('n = %d,',length(allValsG{1}));
    ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.015 -0.45 0 0]);set(ht,'FontSize',7,'Color','k');
    titletxt = sprintf('%d',length(allValsG{2}));
    ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.067 -0.45 0 0]);set(ht,'FontSize',7,'Color','r');
end
for ci = 1:4
    distD = [var_CvD(:,ci) var_AvD(:,ci)];
    tcolors = {'k','r'};
    [distDo,allVals,allValsG] = plotDistributions(distD);
    axes(ff.h_axes(1,ci+4)); hold on;
    ysp = 3;
    mVar(1,1) = mean(allValsG{1}); mVar(2,1) = mean(allValsG{2});
    semVar(1,1) = std(allValsG{1})/sqrt(length(allValsG{1})); semVar(2,1) = std(allValsG{2})/sqrt(length(allValsG{2}));
    mY = 0; MY = 15;%max([mVar+semVar]);
    tcolors = {'k','r'};%{colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
    xdata = make_xdata([2],[1 2]);   combs = [1 2];
    [h,p,tstat] = ttest2(allValsG{1},allValsG{2});
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',9,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); xticks = xdata; 
    if ci == 1
        hy = ylabel({'Avg. # of Trans.','per min'}); changePosition(hy,[-0.02 0 0]);
        set(gca,'ytick',[mY MY],'yticklabels',{'0','15'});
    end
    xticklabels = {'C-TG','A-TG'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    make_bars_hollow(hbs(1:end));
%     put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
    format_axes_b(gca);
    titletxt = sprintf('C%d',ci);
    ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0.021 0.1 0 0]);set(ht,'FontSize',8);
    ht = set_axes_top_text_no_line(gcf,gca,'t-Test',[0.0 -0.01 0 0]);set(ht,'FontSize',7);
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('Distribution_firing_rate'),600);

