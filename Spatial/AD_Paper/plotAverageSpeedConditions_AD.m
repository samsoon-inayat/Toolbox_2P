function plotAverageSpeedConditions(b,markers1,markers2,fn)
%%
ei_C = evalin('base','ei10_C');
ei_A = evalin('base','ei10_A');

selContexts1 = [1 2 3 4 1 2 3 4 ];
rasterNames1 = {'airT','airT','airT','airT','airIT','airIT','airIT','airIT'};
[~,ei_C] = get_rasters_data(ei_C,selContexts1,rasterNames1); 
[~,ei_A] = get_rasters_data(ei_A,selContexts1,rasterNames1); 

%%
mData = evalin('base','mData');
colors = mData.colors;

n = 0;
%%
out = find_speeds(ei_C)
moas = out.meanSpeedTrialsAnimalsT;
moasi = out.meanSpeedTrialsAnimalsIT;

out = find_speeds(ei_A)
moas_A = out.meanSpeedTrialsAnimalsT;
moasi_A = out.meanSpeedTrialsAnimalsIT;

%% R1 (during revision)
var_C = [moas(:,1) moasi(:,1) moas(:,2) moasi(:,2) moas(:,3) moasi(:,3) moas(:,4) moasi(:,4)]; 
var_A = [moas_A(:,1) moasi_A(:,1) moas_A(:,2) moasi_A(:,2) moas_A(:,3) moasi_A(:,3) moas_A(:,4) moasi_A(:,4)];
[within,dvn,xlabels] = make_within_table({'Conds','Trials'},[4,2]);
dataT = make_between_table({var_C;var_A},dvn);
ra1 = RMA(dataT,within,{0.05,{'bonferroni','hsd','lsd'}});
ra1.ranova
print_for_manuscript(ra1)

%%
    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[10 3 3.4 1],'RowsCols',[1 2],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.35],'widthHeightAdjustment',[10 -410]);
    MY = 40; ysp = 7; mY = 0; titletxt = ''; ylabeltxt = {'Average','Speed (cm/s)'};
    stp = 0.36*magfac; widths = ([1.9 1 1.3 1.3 1.3 0.5 0.5 0.5])*magfac; gap = 0.1*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{''});

   [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra1,{'Group_by_Conds_Trials','hsd'},[1.5 1 1]);
    xdata = make_xdata([2 2 2 2 2 2 2 2],[1 1.5 1]);   
    axes(ff.h_axes(1,1));    hold on;
    tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
    tcolors = repmat(tcolors,2,1)';
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,[],[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);

    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 MY],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = xdata(1:2:end)+0.5; xticklabels = {'C1','C2','C3','C4','C1','C2','C3','C4'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(0);
    for ii = 2:2:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,8,{'C-TG','A-TG'});
%     rectangle(gca,'Position',[0.75 30.5 1 3.5],'edgecolor','k','facecolor','k');
%     text(1.85,30.5,'Trials','FontSize',5);
%     rectangle(gca,'Position',[6 30.5 1 3.5],'edgecolor','k');
%     text(7.2,30.5,'Inter-Trials','FontSize',5);
% For pooled
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra1,{'Conds_by_Trials','hsd'},[1.5 1 1]);
    rowsToOmit = [11 13 20 22 23 28];
    h(rowsToOmit) = 0;
    xdata = make_xdata([2 2 2 2],[1 1.5]);   
    axes(ff.h_axes(1,2));
    hold on;
    colors = mData.colors;
    tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3}};
    tcolors = repmat(tcolors,2,1)';
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);

    set(gca,'xlim',[0.25 max(xdata)+0.75],'ylim',[0 MY],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = [1.5 4.5 7.5]; xticklabels = {'C1','C2','C3','C4'}; 
    set(gca,'xtick',xdata(1:2:end)+0.5,'xticklabels',xticklabels);
    xtickangle(0);
    for ii = 2:2:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,8,{'Pooled'});
    save_pdf(ff.hf,mData.pdf_folder,'Figure_1_behavior_anova_speed.pdf',600);


%%
% writetable(between,'Training_Data.xls');
if 1
mVar = ra.est_marginal_means.Mean;
semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < -0.05;

xdata = [1 2 4 5 7 8 10 11 [13 14 16 17 19 20 22 23]+2]; maxY = 22;
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 2.5 1],'color','w');
hold on;
tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4};colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
hbs = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'maxY',maxY,'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.1,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',-0.1);
for ii = 2:2:length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
end
% plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 40],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
% xticks = [1.5 4.5 7.5 10.5]; xticklabels = {'C1','C2','C3','C4'};
xticks = xdata(1:2:end)+0.5; xticklabels = {'C1','C2','C3','C4','C1','C2','C3','C4'};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
% changePosition(gca,[0.03 0.03 0.02 -0.11]);
changePosition(gca,[0.05 0.02 0 -0.011])
put_axes_labels(gca,{[],[0 0 0]},{{'Average','Speed (cm/s)'},[0 0 0]});
% rectangle(gca,'Position',[0.75 21 1 2],'edgecolor','k','facecolor','k');
% text(1.85,22,'Trials','FontSize',5);
% rectangle(gca,'Position',[6 21 1 2],'edgecolor','k');
% text(7.2,22,'Inter-Trials','FontSize',5);
% text(3.5,18,'Control','FontSize',7);
% text(18.5,18,'APP','FontSize',7);
% applyhatch_plusC(gcf
save_pdf(hf,mData.pdf_folder,'AverageSpeedConditions',600);
return;
end
%%
if 1
    mVar = ra.est_marginal_means_wf.Mean;
    semVar = ra.est_marginal_means_wf.Formula_StdErr;
    combs = ra.mcs_lsd_wf.combs; p = ra.mcs_lsd_wf.p; h = p < 0.05;
    combs = ra.mcs_wf.combs; p = ra.mcs_wf.p; h = p < 0.05;
%     combs = ra.mcs_lsd_wf.combs; p = ra.mcs_lsd_wf.p; h = p < 0.05;
    rowsToOmit = [13 20 22 23 28];
    rowsToOmit = [11 13 20 22 23 28];
    [[1:28]' combs h]
    h(rowsToOmit) = 0;
    xdata = [1 2 4 5 7 8 10 11]; 

    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
    tcolors = repmat(tcolors,2,1)';
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',7,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.7,'sigLinesStartYFactor',0.01);

    set(gca,'xlim',[0.25 max(xdata)+0.75],'ylim',[0 40],'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
    xticks = [1.5 4.5 7.5 10.5]; xticklabels = {'C1','C2','C3','C4'}; 
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30);
    for ii = 2:2:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    changePosition(gca,[0.13 0.02 -0.1 -0.011])
%     put_axes_labels(gca,{[],[0 0 0]},{{'Average','Speed (cm/s)'},[0.1 0 0]});

    save_pdf(hf,mData.pdf_folder,'Figure_1_behavior_anova_speed__testing_pooled.pdf',600);
    get(gca, 'FontName')
return;
end


function out = find_speeds(ei)
for an = 1:length(ei)
    for cc = 1:4
        thisspeed = ei{an}.plane{1}.contexts(cc).rasters.airT.speed;
        meanSpeedTrials(:,cc) = nanmean(thisspeed,2);
        thisspeed = ei{an}.plane{1}.contexts(cc).rasters.airIT.speed;
        meanSpeedTrialsI(:,cc) = nanmean(thisspeed,2);
    end
    
    meanSpeedTrialsAnimalsT(an,:) = mean(meanSpeedTrials);
    meanSpeedTrialsAnimalsIT(an,:) = mean(meanSpeedTrialsI);
end
out.meanSpeedTrialsAnimalsT = meanSpeedTrialsAnimalsT;
out.meanSpeedTrialsAnimalsIT = meanSpeedTrialsAnimalsIT;