function plot_Vars_wrt_centers_AD
%% Old Code ... still works
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_C = evalin('base','ei10_C'); 
ei_A = evalin('base','ei10_A'); 


selContexts = [1 2 3 4];
rasterNames = {'airD','airD','airD','airD'};

RsC = get_rasters_data(ei_C,selContexts,rasterNames);
mRsC = calc_mean_rasters(RsC,1:10);
RsC = find_responsive_rasters(RsC,1:10);
% view_population_vector(Rs,mRs,300);
[resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC] = get_responsive_fraction(RsC);

RsA = get_rasters_data(ei_A,selContexts,rasterNames);
mRsA = calc_mean_rasters(RsA,1:10);
RsA = find_responsive_rasters(RsA,1:10);
% view_population_vector(Rs,mRs,400);
[resp_fractionA,resp_valsA,OIA,mean_OIA,resp_ORA,resp_OR_fractionA,resp_ANDA,resp_AND_fractionA] = get_responsive_fraction(RsA);

%% Prepare data
RsC = Rs_C; RsA = Rs_A;
ntrials = 50;
props_C = get_props_Rs(RsC,ntrials); props_A = get_props_Rs(RsA,ntrials);
pop_var_name = {'all','vals','valsT','Nvals','good_zMI','Ngood_zMI'};
pop_var_name = {'good_zMI','good_Gauss','good_MFR'};
pop_var_name = {'good_zMI'};
sel_pop_CA = cell_list_op(props_C,pop_var_name); sel_pop_AA = cell_list_op(props_A,pop_var_name);
colors = mData.colors;
n = 0;
%% Run the Statistical test

all_variables = {'all_zMIs','all_fFR','all_fwidths','all_frs','P'};
% ylabels = {{'Mutual Information','(z-score)'},'Firing Rate (AU)','Field Widths (cm)','R-Squared',{'Spatially Tuned', 'Cells (%)'}};

vn = 3;

number_of_bins = 3;
[all_valsC,all_vals_NC] = get_values(RsC,number_of_bins,all_variables{vn},sel_pop_CA);
[all_valsA,all_vals_NA] = get_values(RsA,number_of_bins,all_variables{vn},sel_pop_AA);

if 1
    all_valsC = all_vals_NC;
    all_valsA = all_vals_NA;
    vn = 5;
end
if number_of_bins > 1
    ind = 1;
    w1 = [];
    w2 = [];
    for ii = 1:(size(all_valsC,2)/number_of_bins)
        for jj = 1:number_of_bins
            varNames{ind} = sprintf('C%dB%d',ii,jj);
            xticklabels{ind} = sprintf('B%d',jj);
            temp_tcolors{ind} = mData.colors{ii};
            w1 = [w1 ii];
            w2 = [w2 jj];
            ind = ind + 1;
        end
    end
    
    [within,dvn,xlabels] = make_within_table({'Cond','Bin'},[4,number_of_bins]);
    dataT = make_between_table({all_valsC;all_valsA},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
    ra.ranova
else
    temp_tcolors = repmat(colors(1:4),2,1);
    varNames = {'C1','C2','C3','C4'};
    xticklabels = varNames;
    dataT = array2table([[1;1;1;1;1;2;2;2;2;2] [all_valsC;all_valsA]]);
    dataT.Properties.VariableNames = {'Group',varNames{:}};
    within = array2table([1 2 3 4]');
    within.Properties.VariableNames = {'Cond'};
    within.Cond = categorical(within.Cond);
    dataT.Group = categorical(dataT.Group);
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    ra.ranova
end
% writetable(dataT,fullfile(mData.pdf_folder,sprintf('%s_values.xlsx',all_variables{vn})));
n = 0;
%% average distributions w.r.t centers for the two groups

magfac = mData.magfac;
ff = makeFigureRowsCols(108,[5 5 6.9 1.25],'RowsCols',[1 2],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.39],'widthHeightAdjustment',[10 -520]);
switch vn
    case 5 % responsive cells 
        MY = 150; ysp = 10; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'Percent of Cells'};
    case 2
        MY = 70; ysp = 3; mY = 0; titletxt = 'Response Fidelity'; ylabeltxt = {'Percent of Trials'};
    case 4
        MY = 0.7; ysp = 3; mY = 0; titletxt = ''; ylabeltxt = {'R-squared'};
    case 1
        MY = 10; ysp = 0.5; mY = 0; titletxt = 'Mutual Information'; ylabeltxt = {'Z-Score'};% for all cells (vals) MY = 70
end
stp = 0.3;magfac; widths = ([4.1 1.9 1.3 1.3 1.3 0.5 0.5 0.5]+0.18)*magfac; gap = 0.1*magfac;
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
ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.051 0 0]);
%     changePosition(gca,[-0.03 0.03 0.11 -0.01]);
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'C1','C2','C3','C4','C1','C2','C3','C4'},{[-0.01 0.02]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,12,{'Control','APP'},{[-0.12 0]});

[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_Bin','hsd'},[1.5 1 1]);
xdata = make_xdata([3 3 3 3],[1 1.5]);   
%     mVar = ra.est_marginal_means.Mean; semVar = ra.est_marginal_means.Formula_StdErr;
tcolors = [temp_tcolors temp_tcolors];
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
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'C1','C2','C3','C4'},{[-0.02 0.01]});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,12,{'Pooled'},{[-0.14 0]});


save_pdf(ff.hf,mData.pdf_folder,sprintf('%s_distributions_over_belt_%d',all_variables{vn},number_of_bins),600);
%% average distributions w.r.t centers for the two groups all bars and main effect of bins
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[5 5 6.9 1],'RowsCols',[1 3],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.39],'widthHeightAdjustment',[10 -520]);
switch vn
    case 5 % responsive cells 
        MY = 5; ysp = 7; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'Percent of Cells'};
    case 2
        MY = 70; ysp = 3; mY = 0; titletxt = 'Response Fidelity'; ylabeltxt = {'Percent of Trials'};
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
%% average distributions w.r.t centers for the two groups all bars and main effect of bins and Cond
magfac = mData.magfac;
ff = makeFigureRowsCols(108,[5 5 6.9 1],'RowsCols',[1 3],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.39],'widthHeightAdjustment',[10 -520]);
switch varT
    case 5 % responsive cells 
        MY = 80; ysp = 7; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'Percent of Cells'};
    case 2
        MY = 70; ysp = 3; mY = 0; titletxt = 'Response Fidelity'; ylabeltxt = {'Percent of Trials'};
    case 4
        MY = 0.7; ysp = 3; mY = 0; titletxt = ''; ylabeltxt = {'R-squared'};
    case 1
        MY = 10; ysp = 0.5; mY = 0; titletxt = 'Mutual Information'; ylabeltxt = {'Z-Score'};% for all cells (vals) MY = 70
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
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,3,{'Pooled'},{[0 0]});


save_pdf(ff.hf,mData.pdf_folder,sprintf('%s_distributions_over_belt_%d',all_variables{vn},number_of_bins),600);


