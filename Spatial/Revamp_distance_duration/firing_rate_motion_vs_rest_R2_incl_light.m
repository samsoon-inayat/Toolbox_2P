function firing_rate_motion_vs_rest

Rs_C = o.Rs;% get_rasters_data(ei_C,selContexts,rasterNames);
% typeP = {'all','vals'
typeP  = 'all';
thr = -1;
%%
ntrials = 50;
si = [Lb_T Ab_t_T Ab_i_T Ar_t_D Ar_i_T ArL_t_D ArL_i_T Ars_t_D Ars_i_T Lbs_T Abs_t_T Abs_i_T];
Rs_C = o.Rs(:,si); mRs_C = o.mR(:,si); 
props_C = get_props_Rs(o.Rs(:,si),ntrials);
% pop_var_name = {'all','vals','valsT','Nvals','good_zMI','Ngood_zMI'};
pop_var_name = {'good_zMI','good_Gauss','good_MFR'};
pop_var_name = {'all'};
% pop_var_name = {'vals'};
% pop_var_name = {'vals','good_zMI'};
sel_pop_C = cell_list_op(props_C,pop_var_name);

disp('Done')
%% For configurations C2, C3, C4, C5, and C7 ... FR average over trials and cells and animals from raster plots
MFR = []; mFR = [];
for sii = 1:length(si)
    avgFR = [];
    for ani = 1:5
        tR = Rs_C{ani,sii};
        firing_rate = tR.sp_rasters1;
        cell_mean = mean(firing_rate, [1, 2], 'omitnan');  % Mean for each cell (dimension 3)
        cell_std = std(firing_rate, 0, [1, 2], 'omitnan');  % Std for each cell

        % Use repmat to expand the mean and std to match the dimensions of the original matrix
        mean_expanded = repmat(cell_mean, [size(firing_rate, 1), size(firing_rate, 2), 1]);
        std_expanded = repmat(cell_std, [size(firing_rate, 1), size(firing_rate, 2), 1]);
        
        % Z-score the matrix (subtract mean and divide by std)
        z_firing_rate = (firing_rate - mean_expanded) ./ std_expanded;
        
        % Handle cells where standard deviation is 0 (to avoid division by zero)
        z_firing_rate(std_expanded == 0) = 0;  % Set Z-score to 0 if the standard deviation is zero

        avgFR(ani,:) = (nanmean(squeeze(nanmean(z_firing_rate,1)),2))';
    end
    all_avgFR{sii} = avgFR;
    mall_avgFR{sii} = mean(avgFR);
    MFR = [MFR max(mall_avgFR{sii})];  mFR = [mFR min(mall_avgFR{sii})];
    xsFR{sii} = tR.xs;
end

MMFR = max(MFR);
mmFR = min(mFR);

xlabelsT = {'s','s','s','cm','s','cm','s','cm','s','s','s','s'};
aoo = {'Light','Air-On','Air-Off','Air-On','Air-Off','Air-On','Air-Off','Air-On','Air-Off','Light','Air-On','Air-Off'};
% colors = {'b','c','b','c','b','c','b','c','b','c'};
colors = [];
for ii = 1:7
    colors = [colors;[mData.colors(ii);mData.colors(ii)]]
end

magfac = mData.magfac;
ff = makeFigureRowsCols(108,[3 3 5.8 1.5],'RowsCols',[1 12],'spaceRowsCols',[0.01 -0.02],...
    'rightUpShifts',[0.31 0.4],'widthHeightAdjustment',[10 -650]);
MY = MMFR; ysp = 0.025; mY = mmFR; titletxt = ''; ylabeltxt = {'Avg. FR (A.U.)'};

stp = 0.3*magfac; widths = ([0.3*ones(1,12)])*magfac; gap = 0.16*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,ylabeltxt);

for gii = 1:12
    axes(ff.h_axes(1,gii)); ha = gca;
    tavgFR = all_avgFR{gii};
    m_avgFR = mean(all_avgFR{gii});  sem_avgFR = std(all_avgFR{gii})/sqrt(5);
    shadedErrorBar(1:length(m_avgFR),m_avgFR,sem_avgFR,{'color',colors{gii},'linewidth',0.25},0.7);
%     shadedErrorBar(xsFR{gii},m_avgFR,sem_avgFR,{'color',colors{gii},'linewidth',0.25},0.7);
    xlim([1 length(m_avgFR)]);ylim([mmFR MMFR]);
    if gii > 1
        set(gca,'YTick',[]);
    else
        ylabel(ylabeltxt);
    end
    xlabel(xlabelsT{gii});
    [within,dvn,xlabels] = make_within_table({'B'},[size(tavgFR,2)]);
    dataT = make_between_table({tavgFR},dvn);
    ra = RMA(dataT,within,{0.05,{''}});
%     ra.ranova;
    print_for_manuscript(ra)
    p_vals_FR(gii) = ra.ranova{3,ra.selected_pval_col};
    titletxt = sprintf('%s',aoo{gii});
    ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0 0.09 0 0]);set(ht,'FontSize',7);
    titletxt = sprintf('%s',getNumberOfAsterisks(p_vals_FR(gii)));
    ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.015 -0.05 0 0]);set(ht,'FontSize',7);
    format_axes(gca)
end

% set_sub_graph_text(ff,2,{'C1','C2','C3','C4','C5','C6','C7'},[0 0 0 0 0],[0.0 0.0 0 0 0]);
% set_sub_graph_text(ff,1,titletxts,[0 0.55 0 0],[0.02 0.07 0 0]);

for gii = 1:12
    axes(ff.h_axes(1,gii)); ha = gca;
    m_avgFR = mean(all_avgFR{gii});
    set(ha,'xtick',[1 length(m_avgFR)]);
    if gii == 1 || gii == 10
        set(ha,'xticklabel',{'0','2'});
    end
    if gii == 2 || gii == 11
        set(ha,'xticklabel',{'0','5'});
    end
    if ismember(gii,[3 12])
        set(ha,'xticklabel',{'0','10'});
    end
    if ismember(gii,[3 5 7])
        set(ha,'xticklabel',{'0','150'});
    end
    if ismember(gii,[4 6 8])
        set(ha,'xticklabel',{'0','15'});
    end
    format_axes(gca)
end

save_pdf(ff.hf,mData.pdf_folder,'allFR.pdf',600);

%% run to get firing rates rest vs motion % added on Sep 13 to check the difference between locomotion and rest periods firing rate during air-off phase
out_C_off = get_spike_rate_air_off(ei,thr,sel_pop_C);
out_C_on = get_spike_rate_air_on(ei,thr,sel_pop_C);
% %% cell level stats individually for separate animals
% ani = 5;
% data_C = [out_C_on.m_sp_animal_motion{ani} out_C_on.m_sp_animal_rest{ani} out_C_off.m_sp_animal_motion{ani} out_C_off.m_sp_animal_rest{ani}];
% [within,dvn,xlabels] = make_within_table({'Ph','St'},[2,2]);
% dataT = make_between_table({data_C},dvn);
% ra = RMA(dataT,within,{0.05,{'hsd'}});
% %     ra.ranova
% print_for_manuscript(ra);

%% animal level stats
data_C = [out_C_on.m_sp_animal_level_motion' out_C_on.m_sp_animal_level_rest' out_C_off.m_sp_animal_level_motion' out_C_off.m_sp_animal_level_rest'];
[within,dvn,xlabels] = make_within_table({'Ph','St'},[2,2]);
dataT = make_between_table({data_C},dvn);
ra = RMA(dataT,within,{0.05,{'hsd'}});
%     ra.ranova
print_for_manuscript(ra);
%%
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[10 3 2.45 1],'RowsCols',[1 2],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.31 0.40],'widthHeightAdjustment',[10 -450]);
MY = 0.25; ysp = 0.05; mY = 0; titletxt = ''; ylabeltxt = {'Avg. Peak','FR (A.U.)'};
stp = 0.4*magfac; widths = ([1 0.4 0.4 1.3 1.3 0.5 0.5 0.5])*magfac; gap = 0.0251*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
colors = mData.colors;
tcolors = [colors(1);colors(1);colors(2);colors(2);colors(3);colors(3);colors(4);colors(4);colors(5);colors(5)];
axes(ff.h_axes(1,1));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Ph','hsd'},[1.5 1 1]);
% xdata = make_xdata([2 2 2 2 2],[1 1.5]); 
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',9,'barWidth',0.5,'sigLinesStartYFactor',0.07,'sigAsteriskyshift',0.015);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'Motion','Rest'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30);
make_bars_hollow(hbs(2:2:end));
put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'AOn','AOff','C4','C5','C7'},{[-0.03 0.01]});
% ht = set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,8,{'C-TG','A-TG'},{[-0.14 -0]}); 
format_axes(gca);

tcolors = [mData.dcolors(1);mData.dcolors(1)];
axes(ff.h_axes(1,2));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Ph','hsd'},[1.5 1 1]);
    xdata = make_xdata([2],[1 2]);   
%         combs = [[1:2:12]' [2:2:12]']; p = ra.MC.bonferroni.Group_by_Cond{1:2:12,6}; h = p<0.05;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.45,'sigLinesStartYFactor',0.02,...
    'sigAsteriskyshift',0.01);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'AOn','AOff'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30);
make_bars_hollow(hbs(2));
ht = set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[-0.05 -0]}); 
format_axes(gca);

save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);


%% run to get firing rates rest vs motion
out_C = get_spike_rate(ei,thr,sel_pop_C);
%% plot distributions FR rest vs motion
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[10 3 1.35 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.31 0.29],...
    'widthHeightAdjustment',[10 -450]);
MY = 0.15; ysp = 0.03; mY = 0; titletxt = ''; ylabeltxt = {'Avg. FR','(Z-Score)'};
stp = 0.3*magfac; widths = ([0.8 0.4 0.4 1.3 1.3 0.5 0.5 0.5])*magfac; gap = 0.051*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
hold on;
distD = [out_C.m_sp_animal_rest' out_C.m_sp_animal_motion'];
plotDistributions(distD,'ks')
tcolors = {'b','m'};
[distDo,allVals,allValsG] = plotDistributions(distD);
minBin = min(allVals);
maxBin = max(allVals);
incr = 0.01;
% [ha,hb,hca] = plotDistributions(allValsG,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'do_mean','No');
[ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'do_mean','Yes');
set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
%     changePosition(gca,[0.129 0.15 -0.09 -0.13]);
ylim([0 100]); xlim([minBin maxBin]); xlim([minBin 0.25]);
put_axes_labels(ha,{{'Avg. Firing Rate (Z-score)'},[0 0 0]},{{'Cells (%)'},[0 0 0]});
format_axes(ha);
[ks.h,ks.p,ks.ks2stat] = kstest2(allValsG{1},allValsG{2});
ks.DF1 = length(allValsG{1}); ks.DF2 = length(allValsG{2});
print_for_manuscript(ks,'KS2')
ht = set_axes_top_text_no_line(gcf,ha,'KS-Test',[0.0 -0.01 0.1 0]);set(ht,'FontSize',7);
titletxt = sprintf('%s',getNumberOfAsterisks(ks.p));
ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.3 -0.01 0.1 0]);set(ht,'FontSize',9);
% titletxt = sprintf('n = %d,',length(allValsG{1}));
% ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.015 -0.45 0 0]);set(ht,'FontSize',7,'Color','k');
% titletxt = sprintf('%d',length(allValsG{2}));
% ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.067 -0.45 0 0]);set(ht,'FontSize',7,'Color','r');
save_pdf(ff.hf,mData.pdf_folder,'firing_rate.pdf',600);
   
%% rest vs motion FR average

    tcolors = {'k','r'};
    data_C = [out_C.m_sp_animal_level_rest' out_C.m_sp_animal_level_motion'];
    [within,dvn,xlabels] = make_within_table({'St'},[2]);
    dataT = make_between_table({data_C},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd'}});
%     ra.ranova
    print_for_manuscript(ra);

   magfac = mData.magfac;
ff = makeFigureRowsCols(108,[10 3 1 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.31 0.25],'widthHeightAdjustment',[10 -350]);
MY = 0.1; ysp = 0.03; mY = -0.051; titletxt = ''; ylabeltxt = {'Avg. FR (Z-score)'};
stp = 0.45*magfac; widths = ([0.4 0.4 0.4 1.3 1.3 0.5 0.5 0.5])*magfac; gap = 0.051*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
colors = mData.colors;
tcolors = {'b','m'};%repmat({colors{1};colors{2};colors{3};colors{4}},4);
axes(ff.h_axes(1,1));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'St','hsd'},[1.5 1 1]);
xdata = make_xdata([2],[1 1.5]); 
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.025,...
'sigAsteriskyshift',0.01);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca); xticks = xdata; 
xticklabels = {'Rest','Motion'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30);
make_bars_hollow(hbs(9:end));
put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Air-On','Air-Off','Air-On','Air-Off'},{[-0.01 0.03]});
% ht = set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,8,{'C-TG','A-TG'},{[-0.14 -0]}); 
format_axes(gca);
save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%% spike rate
[out_C,out_C_motion,out_C_rest] = get_spike_rate_ext_dist_dur(ei,thr,sel_pop_C,selContexts(si));
out_CM = get_spike_rate_ext_dist_dur_MAX(ei,thr,sel_pop_C,selContexts(si));
% out_C_Ca = get_spike_rate_ext_Ca(ei_C,thr,sel_pop_C);
disp('Done');
%%
varC_motion = exec_fun_on_cell_mat(out_C_motion(:,3:8),'nanmean');
varC_rest = exec_fun_on_cell_mat(out_C_rest(:,3:8),'nanmean');
[within,dvn,xlabels] = make_within_table({'St','Cond','Ph'},[2,3,2]);
dataT = make_between_table({[varC_motion varC_rest]},dvn);
ra = RMA(dataT,within,{0.05,{'hsd'}});
ra.ranova
print_for_manuscript(ra)
%%
sel_c  = [1:10];
varC = exec_fun_on_cell_mat(out_C,'mean');
[within,dvn,xlabels] = make_within_table({'Cond','Ph'},[5,2]);
dataT = make_between_table({varC(:,sel_c)},dvn);
ra = RMA(dataT,within,{0.05,{'hsd'}});
ra.ranova
print_for_manuscript(ra)

%% for NB conditions
sel_c  = [1 2 9 10];
avdata = [varC(:,[7 8]) varC(:,[9 10])];
avdata = [varC(:,[3 4 5 6 7 8])];
% avdata = [varC(:,[1 2 9 10])];
varC = exec_fun_on_cell_mat(out_C,'mean');
[within,dvn,xlabels] = make_within_table({'Cond','Ph'},[3,2]);
dataT = make_between_table({avdata},dvn);
ra = RMA(dataT,within,{0.05,{'hsd'}});
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
[within,dvn,xlabels] = make_within_table({'Cond','Ph'},[5,2]);
dataT = make_between_table({varC(:,:)},dvn);
ra = RMA(dataT,within,{0.05,{'hsd'}});
ra.ranova
print_for_manuscript(ra)


%%
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[10 3 4.6 1],'RowsCols',[1 3],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.31 0.40],'widthHeightAdjustment',[10 -450]);
MY = 0.05; ysp = 0.02; mY = 0; titletxt = ''; ylabeltxt = {'Avg. FR (A.U.)'};
stp = 0.3*magfac; widths = ([2 1 0.4 0.4 1.3 1.3 0.5 0.5 0.5])*magfac; gap = 0.051*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
colors = mData.colors;
tcolors = [colors(1);colors(1);colors(2);colors(2);colors(3);colors(3);colors(4);colors(4);colors(5);colors(5)];
axes(ff.h_axes(1,1));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Ph_by_Cond','hsd'},[1.5 1 1]);
xdata = make_xdata([2 2 2 2 2],[1 1.5]); 
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,[],[h p],'colors',tcolors,'sigColor','k',...
'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.0015,...
'sigAsteriskyshift',0);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'AOn','AOff'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30);
make_bars_hollow(hbs(2:2:end));
put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'C2','C3','C4','C5','C7'},{[-0.03 0.01]});
% ht = set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,8,{'C-TG','A-TG'},{[-0.14 -0]}); 
format_axes(gca);

tcolors = [mData.dcolors(1);mData.dcolors(1)];
axes(ff.h_axes(1,2));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Ph','hsd'},[1.5 1 1]);
    xdata = make_xdata([2],[1 2]);   
%         combs = [[1:2:12]' [2:2:12]']; p = ra.MC.bonferroni.Group_by_Cond{1:2:12,6}; h = p<0.05;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.45,'sigLinesStartYFactor',0.02,...
    'sigAsteriskyshift',0.01);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'AOn','AOff'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30);
make_bars_hollow(hbs(2));
ht = set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[-0.05 -0]}); 
format_axes(gca);


save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%%
varC = exec_fun_on_cell_mat(out_CM,'mean');
% varCR = [varC(:,[1:2:10]) varC(:,[2:2:10])];
[within,dvn,xlabels] = make_within_table({'Cond','Ph'},[5,2]);
dataT = make_between_table({varC(:,:)},dvn);
ra = RMA(dataT,within,{0.05,{'hsd'}});
ra.ranova
print_for_manuscript(ra)
%%
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[10 3 2.45 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.31 0.40],'widthHeightAdjustment',[10 -450]);
MY = 35; ysp = 0.02; mY = 0; titletxt = ''; ylabeltxt = {'Avg. Peak','FR (A.U.)'};
stp = 0.4*magfac; widths = ([2 0.4 0.4 1.3 1.3 0.5 0.5 0.5])*magfac; gap = 0.051*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
colors = mData.colors;
tcolors = [colors(1);colors(1);colors(2);colors(2);colors(3);colors(3);colors(4);colors(4);colors(5);colors(5)];
axes(ff.h_axes(1,1));
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Ph_by_Cond','hsd'},[1.5 1 1]);
xdata = make_xdata([2 2 2 2 2],[1 1.5]); 
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,[],[h p],'colors',tcolors,'sigColor','k',...
'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',7,'barWidth',0.5,'sigLinesStartYFactor',0.05);
set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
xticklabels = {'AOn','AOff'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30);
make_bars_hollow(hbs(2:2:end));
put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'C2','C3','C4','C5','C7'},{[-0.03 0.01]});
% ht = set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,8,{'C-TG','A-TG'},{[-0.14 -0]}); 
format_axes(gca);


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

