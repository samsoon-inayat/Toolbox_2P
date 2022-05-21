function air_place_cells_IT_figure

while 1
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei = evalin('base','ei'); 

selContexts = [3 4 5];
rasterNames = {'airIT','airIT','airIT'};
% Rs = get_rasters_data(ei_11_15,selContexts,rasterNames);
RsT = get_rasters_data(ei,selContexts,rasterNames);
% Rs = filterRasters(Rs);
RsT = find_responsive_rasters(RsT,1:9);
% [resp_fractionCS,resp_valsCS,OIC,mean_OICS,resp_ORCS,resp_OR_fractionCS,resp_ANDCS,resp_AND_fractionCS] = get_responsive_fraction(Rs);
resp_FRT = get_responsive_fraction_FR(RsT);
respAT = get_responsive_cells(RsT);
mRT = calc_mean_rasters(RsT,1:9);


selContexts = [3 4 5];
rasterNames = {'airID','airID','airID'};
% Rs = get_rasters_data(ei_11_15,selContexts,rasterNames);
Rs = get_rasters_data(ei,selContexts,rasterNames);
% Rs = filterRasters(Rs);
Rs = find_responsive_rasters(Rs,1:9);
% [resp_fractionCS,resp_valsCS,OIC,mean_OICS,resp_ORCS,resp_OR_fractionCS,resp_ANDCS,resp_AND_fractionCS] = get_responsive_fraction(Rs);
resp_FR = get_responsive_fraction_FR(Rs);
respA = get_responsive_cells(Rs); % [resp_fractionCSA,resp_valsCSA,OICA,mean_OICSA,resp_ORCSA,resp_OR_fractionCSA,resp_ANDCSA,resp_AND_fractionCSA] = get_responsive_fraction(Rs(:,1:3));
% [resp_fractionCSB,resp_valsCSB,OICB,mean_OICSB,resp_ORCSB,resp_OR_fractionCSB,resp_ANDCSB,resp_AND_fractionCSB] = get_responsive_fraction(Rs(:,4:6));
mR = calc_mean_rasters(Rs,1:9);

break
end
n = 0;

%% binarize cells by finding which ones have higher mutual information for distance or time
while 1
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            R = Rs{rr,cc}; zMI_D = R.info_metrics.ShannonMI_Zsh'; zMIsC{rr,cc} = zMI_D;
            R = RsT{rr,cc}; zMI_T = R.info_metrics.ShannonMI_Zsh'; zMIsCT{rr,cc} = zMI_T;
            MI_D_g_T{rr,cc} = zMI_D > zMI_T; MI_T_g_D{rr,cc} = zMI_T > zMI_D;
            D_g_T(rr,cc) = 100*sum(MI_D_g_T{rr,cc} == 1)/length(zMI_D); T_g_D(rr,cc) = 100*sum(MI_T_g_D{rr,cc} == 1)/length(zMI_D);
            diff_D_T{rr,cc} = zMI_D - zMI_T;
            resp_diff_D_g_T{rr,cc} = diff_D_T{rr,cc} > 1.65; resp_diff_T_g_D{rr,cc} = diff_D_T{rr,cc} < -1.65;
        end
    end
    
    
    CN = 1;
    tcolors = {'c'};
    distD(:,1) = diff_D_T(:,CN);
    [distDo,allVals] = getAveragesAndAllValues(distD);
    minBin = min(allVals);
    maxBin = max(allVals);
    incr = 1; %maxBin =
    hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    %    [ha,hb,hca] = plotDistributions(distD,'colors',tcolors,'maxY',maxBin,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    [ha,hb,hca] = plotAverageDistributions(distD,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'pdf_or_cdf','pdf');
    break;
end

%% Find the peak of mean response location of cells on the belt which have higher distance MI versus which have higher time MI

while 1
     for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            R = Rs{rr,cc}; 
            vals = R.peak_location(resp_diff_D_g_T{rr,cc}); mvalsD(rr,cc) = mean(vals);
            vals = R.peak_location(resp_diff_T_g_D{rr,cc}); mvalsT(rr,cc) = mean(vals);
            mzMIsC(rr,cc) = nanmean(zMIsC{rr,cc}); mzMIsCT(rr,cc) = nanmean(zMIsCT{rr,cc});
        end
     end
     [within,dvn,xlabels] = make_within_table({'DT','Cond'},[2,3]);
     dataT = make_between_table({mvalsD,mvalsT},dvn);
    dataT = make_between_table({mzMIsC,mzMIsCT},dvn);
    ra = RMA(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'DT','bonferroni'},[1 1 1]);
%     s = generate_shades(3); tcolors = s.g;
    tcolors = colors;%mData.colors(1:3);
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
%     make_bars_hollow(hbs(4:6))
    format_axes(gca);
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[0 maxY])
    xticks = [xdata(1:end)]; xticklabels = {'D-A','D-B','T-A','T-B'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.17 0.02 -0.1 -0.011])
%     set_title(gca,'AFoR',[1.3,25],5); set_title(gca,'BFoR',[5.3,25],5,'r');
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
%     put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('peak_location_IT_D_vs_T'),600);
    break;
end
%% Speed Figure
if 1
    for ii = 1:length(ei)
        b1 = ei{ii}.b;
        for jj = 1:10
            alds(ii,jj) = b1.dist(b1.stim_r(jj+10)) - b1.dist(b1.air_puff_r(jj+20));
        end
    end
    ald = round(mean(alds(:)));
     ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 6],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[-0.01 0.13],'widthHeightAdjustment',...
        [-45 -250]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 5 6.5 1]);
    for cn = 1:6
        for an = 1:size(Rs,1)
            R = Rs{an,cn};
            figure(100);clf;imagesc(R.dist);colorbar;
            minDists = min(R.dist,[],2);
        end
    end
    [Y,E] = discretize(1:49,3);
    all_speeds = [];
    for cn = 1:6
        mean_speed_over_trials = [];
        aThisSpeed = [];
        for an = 1:size(Rs,1)
            thisSpeed = nanmean(Rs{an,cn}.speed);
            mean_speed_over_trials(an,:) = thisSpeed;
            for bb = 1:3
                aThisSpeed(an,bb) = nanmean(thisSpeed(Y==bb));
            end
        end
        all_speeds = [all_speeds aThisSpeed];
        axes(ff.h_axes(1,cn));
        hold on;
        xs = Rs{1,2}.xs; N = length(xs);
        mspeed = mean(mean_speed_over_trials(:,1:N)); semspeed = std(mean_speed_over_trials(:,1:N))/sqrt(5);
        plot(xs,mspeed);
        shadedErrorBar(xs,mspeed,semspeed);
        changePosition(gca,[0.1 0.15 -0.05 -0.15]);
        if cn == 1
            put_axes_labels(gca,{'',[0 0 0]},{{'Speed (cm/sec)'},[0 -1 0]});
        end
        plot([ald ald],[0 30],'m:','linewidth',0.5);
        plot([50 50],[0 30],'k--','linewidth',0.25);
        plot([100 100],[0 30],'k--','linewidth',0.25);
        xlabel('Position (cm)');
        ylim([0 30]);
        box off;
        set(gca,'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('speeds_345_IT'),600);
%%
    
    [within,dvn,xlabels] = make_within_table({'FR','Conds','Bins'},[2,3,3]);
    dataT = make_between_table({all_speeds},dvn);
    ra = RMA(dataT,within,0.05);
%%
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR_Bins','bonferroni'},[1 1 1]);
    hf = get_figure(5,[8 7 3.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;%mData.colors([5,6,7]); tcolors = repmat(tcolors,2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
%     make_bars_hollow(hbs(4:6));
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 50]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'B1','B2','B3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.05 -0.05])
    put_axes_labels(gca,{[],[0 0 0]},{'Avg. Speed (cm/sec)',[0 0 0]});
    set_title(gca,'AFoR',[1.3,45],5); set_title(gca,'BFoR',[5.3,45],5);
    save_pdf(hf,mData.pdf_folder,sprintf('avg_speed_anova_AB.pdf'),600);
    ra.ranova
    ra.mauchly
    
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR','bonferroni'},[1 1 1]);
%     s = generate_shades(3); tcolors = s.g;
    tcolors = {'k','k'};%mData.colors(1:3);
     hf = get_figure(500,[7 7 1.25 1]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
    make_bars_hollow(hbs(2));
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[0 25])
    xticks = [xdata(1:end)]; xticklabels = {'AFoR','BFoR'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    
    changePosition(gca,[0.17 -0.07 -0.5 0.1])
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{'Avg. Speed (cm/sec)',[0 -2 0]});
    format_axes(gca);
    save_pdf(hf,mData.pdf_folder,sprintf('avg_speed_anova_AB_pooled'),600);

end
%% 
an = 1; cn = 1;
plotRasters_simplest(Rs{an,cn})
%% Show sample rasters
if 1
   an = 3; cn = 1;
    % plotRasters_simplest(Rs{an,cn})
    % find(respAsB{an,cn});
    tr = respAsB{an,cn};
    ff = makeFigureRowsCols(2021,[0.5 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
        [-75 -475]);
    set(gcf,'color','w'); set(gcf,'Position',[5 4 3.25 1]);
    ff = sample_rasters(Rs{an,cn},[191 11 96 41],ff);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rastersD_A'),600);

    an = 3; cn = 4;
    % plotRasters_simplest(Rs{an,cn})
    % find(resp_valsC{an}(:,cn));
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
        [-75 -475]);
    set(gcf,'color','w'); set(gcf,'Position',[10 4 3.25 1]);
    ff = sample_rasters(Rs{an,cn},[142 217 262 302],ff);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('belt_rastersD_A'),600);
end

%% Mutual Information Time versus Distance Distributions

if 1
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            R = Rs{rr,cc};
            zMIs{rr,cc} = R.info_metrics.ShannonMI_Zsh';
        end
    end
    zMIsC = zMIs(:,1:3); zMIsA = zMIs(:,4:6);
    CN = 3;
    tcolors = {'m','c'};
    distD(:,1) = zMIsC(:,CN);
    distD(:,2) = zMIsA(:,CN);
    [distDo,allVals,allValsG] = plotDistributions(distD);
    minBin = min(allVals);
    maxBin = max(allVals);

    
    incr = 0.001; %maxBin =
    hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    %    [ha,hb,hca] = plotDistributions(distD,'colors',tcolors,'maxY',maxBin,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    [ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
    plot([1.65 1.65],[0 100],'--k');
    myxlims2 = [7 9 13];
    xlims = xlim; xlim([xlims(1) myxlims2(CN)]);
    set(gca,'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
    changePosition(gca,[0.15 0.13 -0.4 -0.13]);
    if CN == 1
        put_axes_labels(gca,{'zMI',[0 0 0]},{{'Percentage','of Neurons'},[0 0 0]});
    else
        put_axes_labels(gca,{'zMI',[0 0 0]},{'',[1 0 0]});
    end
    if CN == 3
        legs = {'Dist (D)','Time (T)'};
        legs{end+1} = [4 2 60 5];
        putLegend(gca,legs,'colors',tcolors)
    end
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_zMI_IT_%d',CN),600);
end
%% Mutual Information Time versus Distance bar graph
mzMIsC = arrayfun(@(x) nanmean(x{1}),zMIsC);
mzMIsA = arrayfun(@(x) nanmean(x{1}),zMIsA);
[within,dvn,xlabels] = make_within_table({'DT','Cond'},[2,3]);
dataT = make_between_table({mzMIsC,mzMIsA},dvn);
ra = RMA(dataT,within,0.05);
[xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'DT_by_Cond','bonferroni'},[1 0.25 1]);
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
hold on;
tcolors = [mData.shades.m;mData.shades.c];
% tcolors = repmat(tcolors,2,1)';
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
set(gca,'xlim',[0.35 xdata(end)+.65],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
xticks = xdata; xticklabels = {'Dist','Time'};%xticklabels = {'C3-D','C4-D','C3''-D','C3-T','C4-T','C3''-T'}; xticklabels = repmat(xticklabels,1,2);
set(gca,'xtick',xticks,'xticklabels',xticklabels);
xtickangle(45)
changePosition(gca,[0.05 0.02 -0.05 -0.011])
put_axes_labels(gca,{[],[0 0 0]},{{'zMI'},[0 0 0]});

save_pdf(hf,mData.pdf_folder,'zMI_bar_graph_all_cells_IT.pdf',600);


%% population vector and correlation single animal
% For belt FoR, exclude cells that are also PCs for Air FoR
if 1
    an = 1;
    ff = makeFigureRowsCols(106,[1 0.5 6 1],'RowsCols',[2 6],...
        'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.07 0.1],'widthHeightAdjustment',...
        [0.01 -60]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 5 7 2]);
    resp = get_cell_list(respAsB,[1;2;3]); resp = get_cell_list(respAsB,[4;5;6]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,respAsB,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[-0.1 1],[]);
    set_obj(ff,{'FontWeight','Normal','FontSize',6,'LineWidth',0.25});
    ht = get_obj(ff,'title'); hyl = get_obj(ff,'ylabel'); changePosition(hyl(1,1),[4 0 0]);
    set_obj(ht,'String',{'Pop. Activity','Pop. Activity','Pop. Activity','Pop. Activity','Pop. Activity','Pop. Activity';'Pop. Correlation','Pop.Correlation','Pop.Correlation','Pop. Correlation','Pop.Correlation','Pop.Correlation'});
    set_obj(ht,{'FontSize',5,'FontWeight','Normal'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_population_vector_corr_IT.pdf'),600);
end
%% average correlation of all animals
if 1
    ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 6],...
        'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.075 0.2],'widthHeightAdjustment',...
        [0.01 -250]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 5 6.99 1]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[-0.1 1],[]);
    set_obj(ff,{'FontWeight','Normal','FontSize',6,'LineWidth',0.25});
    ht = get_obj(ff,'title'); 
    set_obj(ht,'String',{'Avg. Pop. Correlation','Avg. Pop. Correlation','Avg. Pop. Correlation'});
    set_obj(ht,{'FontSize',5,'FontWeight','Normal'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_average_population_vector_corrD.pdf'),600);
end

%% population vector and correlation single animal TRIAL WISE trial wise
% For belt FoR, exclude cells that are also PCs for Air FoR
if 1
    an = 1;
    ff = makeFigureRowsCols(106,[1 0.5 6 1],'RowsCols',[2 10],...
        'spaceRowsCols',[0 -0.07],'rightUpShifts',[0.03 0.1],'widthHeightAdjustment',...
        [+60 -90]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 5 13 5]);
    tR = Rs(an,:)'; [tR,tmR] = separate_into_trials(tR,{1;2;3;4;5;6;7;8;9;10});
    respAsB_T = respAsB(an,:)';respAsB_T = repmat(respAsB_T,1,10);
    cn = 6;
    [CRc,aCRc,mRR] = find_population_vector_corr(tR,tmR,respAsB_T,3);
    ff = show_population_vector_and_corr(mData,ff,tR(cn,:),mRR(cn,:),CRc(cn,:),[],[]);
    for ii = 1:10
        axes(ff.h_axes(2,ii)); plot(nanmean(mRR{cn,ii}));
    end
    set_obj(ff,{'FontWeight','Normal','FontSize',6,'LineWidth',0.25});
    ht = get_obj(ff,'title'); hyl = get_obj(ff,'ylabel'); changePosition(hyl(1,1),[4 0 0]);
%     colormap jet
%     set_obj(ht,'String',{'Pop. Activity','Pop. Activity','Pop. Activity','Pop. Activity','Pop. Activity','Pop. Activity';'Pop. Correlation','Pop.Correlation','Pop.Correlation','Pop. Correlation','Pop.Correlation','Pop.Correlation'});
    set_obj(ht,{'FontSize',5,'FontWeight','Normal'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_population_vector_corr_trial_wise.pdf'),600);
end
%% Percentage of Responsive Cells
if 1
    respT = exec_fun_on_cell_mat(respA.vals,{'sum','length'}); respFA = respT.sum./respT.length;
    respT = exec_fun_on_cell_mat(respB.vals,{'sum','length'}); respFB = respT.sum./respT.length;
    respORA = get_cell_list(respA.vals,[1;2;3]); respORB = get_cell_list(respB.vals,[1;2;3]);
    respT = exec_fun_on_cell_mat(respORA,{'sum','length'}); respF_ORA = mean(respT.sum./respT.length,2);
    respT = exec_fun_on_cell_mat(respORB,{'sum','length'}); respF_ORB = mean(respT.sum./respT.length,2);
    [within,dvn,xlabels] = make_within_table({'FR','Cond'},[2,3]);
    dataT = make_between_table({respFA*100,respFB*100},dvn);
    ra = RMA(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR_by_Cond','bonferroni'},[1 1 1]);
%     s = generate_shades(3); tcolors = s.g;
    tcolors = colors;%mData.colors(1:3);
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
    make_bars_hollow(hbs(4:6))
    format_axes(gca);
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[0 25])
    xticks = [xdata(1:end)]; xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    any_mean = mean(100*respF_ORA);    any_sem = std(100*respF_ORA)/sqrt(5);
    pmchar=char(177); any_text = sprintf('%.0f%c%.0f%%',any_mean,pmchar,any_sem); text(0.75,20,any_text,'FontSize',6);
    any_mean = mean(100*respF_ORB);    any_sem = std(100*respF_ORB)/sqrt(5);
    pmchar=char(177); any_text = sprintf('%.0f%c%.0f%%',any_mean,pmchar,any_sem); text(4.75,20,any_text,'FontSize',6);
    changePosition(gca,[0.17 0.02 -0.1 -0.011])
    set_title(gca,'AFoR',[1.3,25],5); set_title(gca,'BFoR',[5.3,25],5,'r');
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('percentage_PCs_responsive_AB'),600);

    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR','bonferroni'},[1 1 1]);
%     s = generate_shades(3); tcolors = s.g;
    tcolors = {'k','r'};%mData.colors(1:3);
     hf = get_figure(500,[7 7 1.25 1]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[0 25])
    xticks = [xdata(1:end)]; xticklabels = {'AFoR','BFoR'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    
    changePosition(gca,[0.17 -0.07 -0.5 0.1])
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
%     put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned','Cells (%)'},[0 0 0]});
    format_axes(gca);
    save_pdf(hf,mData.pdf_folder,sprintf('percentage_PCs_responsive_AB_pooled'),600);
end

%% Mutual Information comparison
if 1
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            R = Rs{rr,cc};
            zMIsC{rr,cc} = R.info_metrics.ShannonMI_Zsh';
        end
    end
    data = arrayfun(@(x) nanmean(x{1}),zMIsC);
    dataM = arrayfun(@(x) max(x{1}),zMIsC);
    [within,dvn,xlabels] = make_within_table({'FR','Cond'},[2,3]);
    dataT = make_between_table({data},dvn);
    ra = RMA(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR_by_Cond','bonferroni'},[1 1 1]);
%     s = generate_shades(3); tcolors = s.g;
    tcolors = colors;%mData.colors(1:3);
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
    make_bars_hollow(hbs(4:6))
    format_axes(gca);
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[0 1.7])
    xticks = [xdata(1:end)]; xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.17 0.02 -0.1 -0.011])
    set_title(gca,'AFoR',[1.3,1.6],5); set_title(gca,'BFoR',[5.3,1.6],5,'r');
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'zMI'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_AB'),600);

    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR','bonferroni'},[1 1 1]);
%     s = generate_shades(3); tcolors = s.g;
    tcolors = {'k','r'};%mData.colors(1:3);
     hf = get_figure(500,[7 7 1.25 1]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[0 1.7])
    xticks = [xdata(1:end)]; xticklabels = {'AFoR','BFoR'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    
    changePosition(gca,[0.17 -0.07 -0.5 0.1])
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
%     put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned','Cells (%)'},[0 0 0]});
    format_axes(gca);
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_AB_pooled'),600);
end

%% r-squared
if 1
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            R = Rs{rr,cc};
            [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
            zMIsC{rr,cc} = rs';
        end
    end
    data = arrayfun(@(x) nanmean(x{1}),zMIsC);
    dataM = arrayfun(@(x) max(x{1}),zMIsC);
    [within,dvn,xlabels] = make_within_table({'FR','Cond'},[2,3]);
    dataT = make_between_table({data},dvn);
    ra = RMA(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR_by_Cond','bonferroni'},[1 1 1]);
    h(h==1) = 0;
%     s = generate_shades(3); tcolors = s.g;
    tcolors = colors;%mData.colors(1:3);
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
    make_bars_hollow(hbs(4:6))
    format_axes(gca);
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[0 0.8])
    xticks = [xdata(1:end)]; xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.17 0.02 -0.1 -0.011])
    set_title(gca,'AFoR',[1.3,0.8],5); set_title(gca,'BFoR',[5.3,0.8],5,'r');
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'R-Squared'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('rsquared_AB'),600);

    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR','bonferroni'},[1 1 1]);
%     s = generate_shades(3); tcolors = s.g;
    tcolors = {'k','r'};%mData.colors(1:3);
     hf = get_figure(500,[7 7 1.25 1]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[0 0.8])
    xticks = [xdata(1:end)]; xticklabels = {'AFoR','BFoR'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    
    changePosition(gca,[0.17 -0.07 -0.5 0.1])
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
%     put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned','Cells (%)'},[0 0 0]});
    format_axes(gca);
    save_pdf(hf,mData.pdf_folder,sprintf('rsquared_AB_pooled'),600);
end


%% Percentage of Responsive Cells Unique in Conditions
if 1
    percCells = [];
    percCells(:,1) = get_cell_list(resp_valsCS,[1 -2 -3],1)'; percCells(:,2) = get_cell_list(resp_valsCS,[-1 2 -3],1)'; percCells(:,3) = get_cell_list(resp_valsCS,[-1 -2 3],1)';
    within = make_within_table({'Cond'},3);
    dataT = make_between_table({percCells*100},{'C31','C4','C32'});
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
%     s = generate_shades(3); tcolors = s.g;
    tcolors = mData.colors;
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
    format_axes(gca);
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[0 25])
    xticks = [xdata(1:end)]; xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    any_mean = mean(100*resp_OR_fractionCS);    any_sem = std(100*resp_OR_fractionCS)/sqrt(5);
%     pmchar=char(177); any_text = sprintf('%.0f%c%.0f%%',any_mean,pmchar,any_sem); text(0.75,20,any_text,'FontSize',6);
    changePosition(gca,[0.17 0.02 -0.4 -0.011])
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('percentage_PCs_responsive_unique'),600);
end
%% Spatial correlation between adjacent trails
if 1
    [within,dvn,xlabels] = make_within_table({'Condition','TrialPairs'},[3,9]);
    var1 = arrayfun(@(x) mean(x{1}),out1.adj_SP_corr_diag); 
    var2 = arrayfun(@(x) mean(x{1}),out2.adj_SP_corr_diag);
    var3 = arrayfun(@(x) mean(x{1}),out3.adj_SP_corr_diag);
    between = make_between_table({var1,var2,var3},dvn);
    ra = RMA(between,within);
    %%
    hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_line_graph(mData,ra,0);
    tcolors = mData.colors;
    ii = 1; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 2; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linestyle','-.','linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 3; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    legs = {'C3','C4','C3'''};
    legs{end+1} = [1 0.2 0.25 0.25];
    putLegendH(gca,legs,tcolors)
    ylims = ylim;
    set(gca,'xlim',[0.5 xdata(end)+0.5],'ylim',[ylims(1) ylims(2)],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'T1-T2','T2-T3','T3-T4','T4-T5','T5-T6','T6-T7','T7-T8','T8-T9','T9-T10'};
    xticklabels = repmat(xticklabels,1,2);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(45);
    changePosition(gca,[0.05 0.12 0.05 -0.08]);
    put_axes_labels(gca,{'Trial-Pairs',[0 0 0]},{{'Spatial Correlation'},[0 -0.05 0]});
    save_pdf(hf,mData.pdf_folder,'spatial_correlation_trials',600);
end
%% Spatial correlation between adjacent trails (taking mean) 
if 1
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'TrialPairs','bonferroni'},[1 0.5 1]);
    hollowsep = 19;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.75 1],'color','w');
    hold on;
    tcolors = colors(1:9); tcolors = repmat(tcolors,[1 3]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
            'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
            'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    ylims = ylim;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[ylims(1) maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'T1-T2','T2-T3','T3-T4','T4-T5','T5-T6','T6-T7','T7-T8','T8-T9','T9-T10'};
    xticklabels = repmat(xticklabels,1,3);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(45);
%     changePosition(gca,[0.1 0.11 -0.2 -0.05]);
    changePosition(gca,[0.0 0.03 0.08 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatial','Correlation'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'spatial_correlation_trials_bar_graph',600);
end

%% Spatial correlation trails with mean of trials
if 1
    [within,dvn,xlabels] = make_within_table({'Condition','TrialPairs'},[3,10]);
    var1 = []; var2 = []; var3 = [];
    for ii = 1:length(out1.all_SP_corr_diag_with_mean)
        var1(ii,:) = (arrayfun(@(x) mean(x{1}),out1.all_SP_corr_diag_with_mean{ii}))'; 
        var2(ii,:) = (arrayfun(@(x) mean(x{1}),out2.all_SP_corr_diag_with_mean{ii}))'; 
        var3(ii,:) = (arrayfun(@(x) mean(x{1}),out3.all_SP_corr_diag_with_mean{ii}))'; 
    end
    between = make_between_table({var1,var2,var3},dvn);
    ra = repeatedMeasuresAnova(between,within);
    %%
    hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_line_graph(mData,ra,0);
    tcolors = mData.colors;
    ii = 1; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 2; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linestyle','-.','linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 3; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    legs = {'C3','C4','C3'''};
    legs{end+1} = [1 0.2 0.25 0.25];
    putLegendH(gca,legs,tcolors)
    ylims = ylim;
    set(gca,'xlim',[0.5 xdata(end)+0.5],'ylim',[ylims(1) ylims(2)],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'T1-T2','T2-T3','T3-T4','T4-T5','T5-T6','T6-T7','T7-T8','T8-T9','T9-T10'};
    xticklabels = repmat(xticklabels,1,2);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(45);
    changePosition(gca,[0.05 0.12 0.05 -0.08]);
    put_axes_labels(gca,{'Trial-Pairs',[0 0 0]},{{'Spatial Correlation'},[0 -0.05 0]});
    save_pdf(hf,mData.pdf_folder,'spatial_correlation_trials',600);
end
%% Spatial correlation trails with mean of trials (mean bar graph)
if 1
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0.5 1]);
    hollowsep = 19;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 3.25 1],'color','w');
    hold on;
    tcolors = colors(1:10); tcolors = repmat(tcolors,[1 3]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
            'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
            'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    ylims = ylim;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[ylims(1) maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'T1','T2','T3','T4-T5','T5-T6','T6-T7','T7-T8','T8-T9','T9-T10'};
    xticklabels = repmat(xticklabels,1,3);
    set(gca,'xtick',xticks,'xticklabels',xlabels);
    xtickangle(45);
%     changePosition(gca,[0.1 0.11 -0.2 -0.05]);
    changePosition(gca,[-0.02 0.03 0.1 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Correlation'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'spatial_correlation_trials_bar_graph',600);
end
%% Rate remap (Delta FR score) between adjacent trails
if 1
    [within,dvn,xlabels] = make_within_table({'Condition','TrialPairs'},[3,9]);
    var1 = arrayfun(@(x) nanmean(x{1}),out1.adj_RR_SP); 
    var2 = arrayfun(@(x) nanmean(x{1}),out2.adj_RR_SP);
    var3 = arrayfun(@(x) nanmean(x{1}),out3.adj_RR_SP);
    between = make_between_table({var1,var2,var3},dvn);
    ra = repeatedMeasuresAnova(between,within);
    %%
    hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_line_graph(mData,ra,0);
    tcolors = mData.colors;
    ii = 1; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 2; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linestyle','-.','linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 3; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    legs = {'C3','C4','C3'''};
    legs{end+1} = [1.5 0.3 1.02 0.25];
    putLegendH(gca,legs,tcolors)
    ylims = ylim;
    set(gca,'xlim',[0.5 xdata(end)+0.5],'ylim',[ylims(1)-0.01 ylims(2)],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'T1-T2','T2-T3','T3-T4','T4-T5','T5-T6','T6-T7','T7-T8','T8-T9','T9-T10'};
    xticklabels = repmat(xticklabels,1,2);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(45);
    changePosition(gca,[0.065 0.12 0.05 -0.08]);
    put_axes_labels(gca,{'Trial-Pairs',[0 0 0]},{{'\Delta FR Score'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'Delta_FR_score_trials',600);
end
%% Spatial correlation between conditions
if 1
    [within,dvn,xlabels] = make_within_table({'Cond'},2);
    var_CE = arrayfun(@(x) mean(x{1}),outRemap.adj_SP_corr_diag);
    dataT = make_between_table({var_CE},dvn);
    ra = repeatedMeasuresAnova(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0 1]);
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors ={colors{7};colors{8};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.02,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C3-C4','C4-C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     set(hbs(2),'facecolor','none','edgecolor',tcolors{2});
    xtickangle(45);
%     changePosition(gca,[0.2 0.03 -0.55 -0.05]);
    changePosition(gca,[0.1 0.03 -0.4 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatial Correlation'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'air_place_cell corr remap',600);
end

%% Delta FR score (RAte remap) between conditions
if 1
    [within,dvn,xlabels] = make_within_table({'Cond'},2);
    var_CE = arrayfun(@(x) nanmean(x{1}),outRemap.adj_RR_SP);
    dataT = make_between_table({var_CE},dvn);
    ra = repeatedMeasuresAnova(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0 1]);
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors ={colors{7};colors{8};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.02,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+0.1],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C3-C4','C4-C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     set(hbs(2),'facecolor','none','edgecolor',tcolors{2});
    xtickangle(45);
%     changePosition(gca,[0.2 0.03 -0.55 -0.05]);
    changePosition(gca,[0.1 0.03 -0.4 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'\Delta FR Score'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'air_place_cell rate remap (DFR)',600);
end

%% Popuation Vector correlation between conditions
if 1
    [within,dvn,xlabels] = make_within_table({'Cond'},2);
    var_CE = arrayfun(@(x) nanmean(x{1}),outRemap.adj_PV_corr_diag);
    dataT = make_between_table({var_CE},dvn);
    ra = repeatedMeasuresAnova(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0 1]);
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors ={colors{7};colors{8};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY+0.1],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C3-C4','C4-C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     set(hbs(2),'facecolor','none','edgecolor',tcolors{2});
    xtickangle(45);
%     changePosition(gca,[0.2 0.03 -0.55 -0.05]);
    changePosition(gca,[0.2 0.03 -0.4 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Pop. Vec.','Correlation'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'air_place_cell PV corr',600);
end

%%
%% overlap RM bar graph
if 1
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    mOI = NaN(6,6,5);
    mask = triu(ones(6,6),3) & ~triu(ones(6,6),4);
%     mask = logical(zeros(6,6)); mask(1,[2 3 4 5]) = 1;
    ois = [];
    OI_E = get_overlap_index(respAB.vals);
    for ii = 1:5
        C12(ii,1) = OI_E{ii}(1,2);
        C12(ii,2) = OI_E{ii}(2,3);
        mOI(:,:,ii) = OI_E{ii};
        ois(ii,:) = OI_E{ii}(mask);
    end
    figure(100);clf;imagesc(nanmean(mOI,3));colorbar
    [within,dvn,xlabels] = make_within_table({'Cond'},[3]);
    var_CE = ois;
    dataT = make_between_table({var_CE},dvn);
    ra = RMA(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','bonferroni'},[1 0 1]);
    h(h==1) = 0;
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors =colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     set(hbs(2),'facecolor','none','edgecolor',tcolors{2});
    xtickangle(45);
    set_title(gca,'AFoR-BFoR',[0.25,0.15],5)
%     changePosition(gca,[0.2 0.03 -0.55 -0.05]);
    changePosition(gca,[0.1 0.03 -0.25 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Overlap Index'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('air_overlap_stats_place_cells'),600);
end

%% Place Field Properties
props = {'Field Width (cm)','Field Center (cm)','Field FR (AU)'};
fileNames = {'PFW','PFC','PFFR'};
all_incr = [1 1 1];
% legs_v = 
for pri = 3%:3;
if 1
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            R = Rs{rr,cc};
            resp = respAsB{rr,cc};
            [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
            switch pri
                case 1
                    zMIsC{rr,cc} = PWs(resp)';
                case 2
                    zMIsC{rr,cc} = centers(resp)';
                case 3
                    zMIsC{rr,cc} = MFR(resp)';
            end
        end
    end

    tcolors = mData.colors;
    distD = zMIsC;
    [distDo,allVals] = getAveragesAndAllValues(distD);
    minBin = min(allVals);
    maxBin = max(allVals);
    incr = all_incr(pri); %maxBin =
    hf = get_figure(8,[5 7 1.25 1]);
    [ha,hb,~,bins] = plotAverageDistributions(distD,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
    format_axes(gca);
    if pri == 1
        legs = {'C3','C4','C3''',[50 3 70 5]}; putLegend(gca,legs,tcolors);
    end
    changePosition(gca,[0.1 0.13 -0.25 -0.13]);
    put_axes_labels(gca,{props{pri},[0 0 0]},{{'Neurons (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_%s',fileNames{pri}),600);
end

if 1
    mzMIsC = arrayfun(@(x) nanmean(x{1}),zMIsC);
    [within,dvn,xlabels] = make_within_table({'FR','Cond'},[2,3]);
    dataT = make_between_table({mzMIsC},dvn);
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0.25 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.07 0.02 -0.4 -0.011])
    put_axes_labels(gca,{[],[0 0 0]},{props{pri},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells.pdf',fileNames{pri}),600);
end
end

%% Place Field Properties widths vs centers scatter plot and correlation
props = {'Field Width (cm)','Field Center (cm)','Field FR (AU)'};
fileNames = {'PFW','PFC','PFFR'};
all_incr = [1 1 1];

if 1
    for rr = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            R = Rs{rr,cc};
            resp = respAnB{rr,cc};
            [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
            W{rr,cc} = PWs(resp)';
            C{rr,cc} = centers(resp)';
            CW_corr(rr,cc) = corr(C{rr,cc},W{rr,cc});
        end
    end
    tcolors = mData.colors;
    hf = get_figure(8,[5 7 1.25 1]);hold on;
    an = 3; 
    for cn = 1:3
        scatter(C{an,cn},W{an,cn},2,tcolors{cn});
        [x,p] = polyfit(C{an,cn},W{an,cn},1);
        y_hat = polyval(x,C{an,cn});
        plot(C{an,cn},y_hat,'color',tcolors{cn});
    end
    format_axes(gca);
    changePosition(gca,[0.1 0.13 -0.25 -0.13]);
    put_axes_labels(gca,{'Field Center (cm)',[0 0 0]},{'Field Width (cm)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('scatter_%s','PFW_PFC'),600);

    [within,dvn,xlabels] = make_within_table({'FR','Cond'},[2,3]);
    dataT = make_between_table({CW_corr},dvn);
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0.25 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1)-0.1 maxY]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.23 0.02 -0.4 -0.011])
    put_axes_labels(gca,{[],[0 0 0]},{{'Corr. Field Width','Vs Field Center'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells.pdf','CW_corr'),600);
end

%% Place Field Properties on the belt
% respCells1 = get_cell_list(resp_valsCS,[1 -2 -3]);     respCells2 = get_cell_list(resp_valsCS,[-1 2 -3]);     respCells3 = get_cell_list(resp_valsCS,[-1 -2 3]);
% all_respCells = {respCells1,respCells2,respCells3};
props = {'Field Width (cm)',{'Spatially Tuned','Cells (%)'},'Field FR (AU)'};
fileNames = {'PFW','PFC','PFFR'}; yspacing = [10 5 5]; 
pri = 2;
Air = 0;
if 1
    sRs = Rs;
    minBin = 0;     maxBin = 150;     incr = 50; % choosing more than 3 bins, give significant anova but not significant multcompare
    % three bins also make sense because LED in condition 2 comes ON at
    % 110cm
    bins = minBin:incr:maxBin;
    num_cells = []; PFW_belt = []; MFR_belt = [];
    for rr = 1:size(sRs,1)
        for cc = 1:size(sRs,2)
            R = sRs{rr,cc};
            [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
            resp = respAsB{rr,cc}; 
            centers = centers(resp)'; PWs = PWs(resp)'; MFR = MFR(resp)';
            for ii = 1:(length(bins)-1)
                binS = bins(ii); binE = bins(ii+1);
                inds = find(centers >= binS & centers < binE);
                num_cells(rr,cc,ii) = 100*length(inds)/sum(resp);
                PFW_belt(rr,cc,ii) = mean(PWs(inds));
                MFR_belt(rr,cc,ii) = mean(MFR(inds));
            end
        end
    end
    mzMIsC = [];
    for rr = 1:size(sRs,1)
        switch pri
            case 1
                mzMIsC(rr,:) = reshape(squeeze(PFW_belt(rr,:,:))',1,(length(bins)-1)*size(sRs,2));                
            case 2
                mzMIsC(rr,:) = reshape(squeeze(num_cells(rr,:,:))',1,(length(bins)-1)*size(sRs,2));
            case 3
                mzMIsC(rr,:) = reshape(squeeze(MFR_belt(rr,:,:))',1,(length(bins)-1)*size(sRs,2));
        end
    end
    [within,dvn,xlabels] = make_within_table({'FR','Cond','Bins'},[2,3,length(bins)-1]);
    dataT = make_between_table({mzMIsC},dvn);
    ra = RMA(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR_by_Cond_Bins','bonferroni'},[1 1 1]);
    hf = get_figure(5,[5 7 2.5 1]);
    s = generate_shades(length(bins)-1);
    tcolors = [mData.colors(1);mData.colors(1);mData.colors(1);mData.colors(2);mData.colors(2);mData.colors(2);mData.colors(3);mData.colors(3);mData.colors(3)]; tcolors = repmat(tcolors,2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',yspacing(pri),'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    if pri == 2
        set_axes_limits(gca,[0.35 xdata(end)+.65],[0 85]);
    else
        set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);
    end
    format_axes(gca);
    make_bars_hollow(hbs(10:end));
    xticks = xdata; xticklabels = {'C3-B1','C3-B2','C3-B3','C4-B1','C4-B2','C4-B3','C3''-B1','C3''-B2','C3''-B3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 25 50 75]); xtickangle(45)
    if pri == 2
        changePosition(gca,[0.02 0.02 0.01 -0.011])
    else
        changePosition(gca,[0.01 0.02 -0.03 -0.011])
    end
    put_axes_labels(gca,{[],[0 0 0]},{props{pri},[0 -2 0]});
    set_title(gca,'AFoR',[4.25,75],5); set_title(gca,'BFoR',[14.25,75],5)
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells_belt.pdf',fileNames{pri}),600);


%%
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR_by_Bins','bonferroni'},[1 1 1]);
    hf = get_figure(6,[11 7 1.25 1]);
    s = generate_shades(length(bins)-1);
    tcolors = [s.y;s.y];%;s.y;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',15,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set_axes_limits(gca,[0.25 xdata(end)+.75],[0 85]);
    format_axes(gca);
    make_bars_hollow(hbs(4:6))
    xticks = xdata; xticklabels = {'B1','B2','B3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 25 50 75]); xtickangle(45)
    changePosition(gca,[0.2 0.02 -0.2 -0.05])
%     put_axes_labels(gca,{[],[0 0 0]},{props{pri},[0 -2 0]});
    set_title(gca,'AFoR',[1.25,75],5); set_title(gca,'BFoR',[5.25,75],5)
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells_belt_pooled.pdf',fileNames{pri}),600);
    ra.ranova
    ra.mauchly
end

%% Place Field Properties on the belt
% respCells1 = get_cell_list(resp_valsCS,[1 -2 -3]);     respCells2 = get_cell_list(resp_valsCS,[-1 2 -3]);     respCells3 = get_cell_list(resp_valsCS,[-1 -2 3]);
% all_respCells = {respCells1,respCells2,respCells3};
props = {'Field Width (cm)',{'Spatially Tuned','Cells (%)'},'Field FR (AU)'};
fileNames = {'PFW','PFC','PFFR'}; yspacing = [10 5 5]; 
pri = 1;
Air = 0;
if 1
    sRs = Rs;
    minBin = 0;     maxBin = 150;     incr = 50; % choosing more than 3 bins, give significant anova but not significant multcompare
    % three bins also make sense because LED in condition 2 comes ON at
    % 110cm
    bins = minBin:incr:maxBin;
    num_cells = []; PFW_belt = []; MFR_belt = [];
    for rr = 1:size(sRs,1)
        for cc = 1:size(sRs,2)
            R = sRs{rr,cc};
            [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
            resp = respAsB{rr,cc}; 
            centers = centers(resp)'; PWs = PWs(resp)'; MFR = MFR(resp)';
            for ii = 1:(length(bins)-1)
                binS = bins(ii); binE = bins(ii+1);
                inds = find(centers >= binS & centers < binE);
                num_cells(rr,cc,ii) = 100*length(inds)/sum(resp);
                PFW_belt(rr,cc,ii) = mean(PWs(inds));
                MFR_belt(rr,cc,ii) = mean(MFR(inds));
            end
        end
    end
    mzMIsC = [];
    for rr = 1:size(sRs,1)
        switch pri
            case 1
                mzMIsC(rr,:) = reshape(squeeze(PFW_belt(rr,:,:))',1,(length(bins)-1)*size(sRs,2));                
            case 2
                mzMIsC(rr,:) = reshape(squeeze(num_cells(rr,:,:))',1,(length(bins)-1)*size(sRs,2));
            case 3
                mzMIsC(rr,:) = reshape(squeeze(MFR_belt(rr,:,:))',1,(length(bins)-1)*size(sRs,2));
        end
    end
    [within,dvn,xlabels] = make_within_table({'FR','Cond','Bins'},[2,3,length(bins)-1]);
    dataT = make_between_table({mzMIsC},dvn);
    ra = RMA(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR_by_Cond_Bins','bonferroni'},[1 1 1]);
    hf = get_figure(5,[5 7 2.5 1]);
    s = generate_shades(length(bins)-1);
    tcolors = [mData.colors(1);mData.colors(1);mData.colors(1);mData.colors(2);mData.colors(2);mData.colors(2);mData.colors(3);mData.colors(3);mData.colors(3)]; tcolors = repmat(tcolors,2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',yspacing(pri),'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);

        set_axes_limits(gca,[0.35 xdata(end)+.65],[0 45]);

    format_axes(gca);
    make_bars_hollow(hbs(10:end));
    xticks = xdata; xticklabels = {'C3-B1','C3-B2','C3-B3','C4-B1','C4-B2','C4-B3','C3''-B1','C3''-B2','C3''-B3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 25 50]); xtickangle(45)
    changePosition(gca,[0.01 0.02 -0.03 -0.011])
    put_axes_labels(gca,{[],[0 0 0]},{props{pri},[0 -2 0]});
    set_title(gca,'AFoR',[4.25,45],5); set_title(gca,'BFoR',[14.25,45],5)
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells_belt.pdf',fileNames{pri}),600);


%%
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR_by_Cond','bonferroni'},[1 1 1]);
    hf = get_figure(6,[11 7 1.25 1]);
    s = generate_shades(length(bins)-1);
    tcolors = mData.colors(1:3); tcolors = repmat(tcolors,2,1);%[s.m;s.m];%;s.y;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set_axes_limits(gca,[0.25 xdata(end)+.75],[0 45]);
    format_axes(gca);
    make_bars_hollow(hbs(4:6))
    xticks = xdata; xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 25 50]); xtickangle(45)
    changePosition(gca,[0.2 0.02 -0.2 -0.05])
%     put_axes_labels(gca,{[],[0 0 0]},{props{pri},[0 -2 0]});
    set_title(gca,'AFoR',[1.25,45],5); set_title(gca,'BFoR',[5.25,45],5)
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells_belt_pooled.pdf',fileNames{pri}),600);
    ra.ranova
%     ra.mauchly
end

%% Place Field Properties on the belt
% respCells1 = get_cell_list(resp_valsCS,[1 -2 -3]);     respCells2 = get_cell_list(resp_valsCS,[-1 2 -3]);     respCells3 = get_cell_list(resp_valsCS,[-1 -2 3]);
% all_respCells = {respCells1,respCells2,respCells3};
props = {'Field Width (cm)',{'Spatially Tuned','Cells (%)'},'Field FR (AU)'};
fileNames = {'PFW','PFC','PFFR'}; yspacing = [10 5 5]; 
pri = 3;
Air = 0;
if 1
    sRs = Rs;
    minBin = 0;     maxBin = 150;     incr = 50; % choosing more than 3 bins, give significant anova but not significant multcompare
    % three bins also make sense because LED in condition 2 comes ON at
    % 110cm
    bins = minBin:incr:maxBin;
    num_cells = []; PFW_belt = []; MFR_belt = [];
    for rr = 1:size(sRs,1)
        for cc = 1:size(sRs,2)
            R = sRs{rr,cc};
            [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
            resp = respAsB{rr,cc}; 
            centers = centers(resp)'; PWs = PWs(resp)'; MFR = MFR(resp)';
            for ii = 1:(length(bins)-1)
                binS = bins(ii); binE = bins(ii+1);
                inds = find(centers >= binS & centers < binE);
                num_cells(rr,cc,ii) = 100*length(inds)/sum(resp);
                PFW_belt(rr,cc,ii) = mean(PWs(inds));
                MFR_belt(rr,cc,ii) = mean(MFR(inds));
            end
        end
    end
    mzMIsC = [];
    for rr = 1:size(sRs,1)
        switch pri
            case 1
                mzMIsC(rr,:) = reshape(squeeze(PFW_belt(rr,:,:))',1,(length(bins)-1)*size(sRs,2));                
            case 2
                mzMIsC(rr,:) = reshape(squeeze(num_cells(rr,:,:))',1,(length(bins)-1)*size(sRs,2));
            case 3
                mzMIsC(rr,:) = reshape(squeeze(MFR_belt(rr,:,:))',1,(length(bins)-1)*size(sRs,2));
        end
    end
    [within,dvn,xlabels] = make_within_table({'FR','Cond','Bins'},[2,3,length(bins)-1]);
    dataT = make_between_table({mzMIsC},dvn);
    ra = RMA(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR_by_Cond_Bins','bonferroni'},[1 1 1]);
    hf = get_figure(5,[5 7 2.5 1]);
    s = generate_shades(length(bins)-1);
    tcolors = [mData.colors(1);mData.colors(1);mData.colors(1);mData.colors(2);mData.colors(2);mData.colors(2);mData.colors(3);mData.colors(3);mData.colors(3)]; tcolors = repmat(tcolors,2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',yspacing(pri),'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);

        set_axes_limits(gca,[0.35 xdata(end)+.65],[0 15]);

    format_axes(gca);
    make_bars_hollow(hbs(10:end));
    xticks = xdata; xticklabels = {'C3-B1','C3-B2','C3-B3','C4-B1','C4-B2','C4-B3','C3''-B1','C3''-B2','C3''-B3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 7 14]); xtickangle(45)
    changePosition(gca,[0.01 0.02 -0.03 -0.011])
    put_axes_labels(gca,{[],[0 0 0]},{props{pri},[0 0 0]});
    set_title(gca,'AFoR',[4.25,15],5); set_title(gca,'BFoR',[14.25,15],5)
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells_belt.pdf',fileNames{pri}),600);


%%
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR_by_Bins','bonferroni'},[1 1 1]);
    hf = get_figure(6,[11 7 1.25 1]);
    s = generate_shades(length(bins)-1);
    tcolors = mData.colors(1:3); tcolors = repmat(tcolors,2,1);%[s.m;s.m];%;s.y;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set_axes_limits(gca,[0.25 xdata(end)+.75],[0 15]);
    format_axes(gca);
    make_bars_hollow(hbs(4:6))
    xticks = xdata; xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 15]); xtickangle(45)
    changePosition(gca,[0.2 0.02 -0.2 -0.05])
%     put_axes_labels(gca,{[],[0 0 0]},{props{pri},[0 -2 0]});
    set_title(gca,'AFoR',[1.25,15],5); set_title(gca,'BFoR',[5.25,15],5)
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells_belt_pooled.pdf',fileNames{pri}),600);
    ra.ranova
%     ra.mauchly
end

%% Place Field Properties on the belt
% respCells1 = get_cell_list(resp_valsCS,[1 -2 -3]);     respCells2 = get_cell_list(resp_valsCS,[-1 2 -3]);     respCells3 = get_cell_list(resp_valsCS,[-1 -2 3]);
% all_respCells = {respCells1,respCells2,respCells3};
props = {'Field Width (cm)',{'Spatially Tuned','Cells (%)'},'Field FR (AU)','zMI'};
fileNames = {'PFW','PFC','PFFR','zMI'}; yspacing = [10 5 5 1]; 
pri = 4;
Air = 0;
if 1
    sRs = Rs;
    minBin = 0;     maxBin = 150;     incr = 50; % choosing more than 3 bins, give significant anova but not significant multcompare
    % three bins also make sense because LED in condition 2 comes ON at
    % 110cm
    bins = minBin:incr:maxBin;
    num_cells = []; PFW_belt = []; MFR_belt = [];
    for rr = 1:size(sRs,1)
        for cc = 1:size(sRs,2)
            R = sRs{rr,cc};
            tzMI = R.info_metrics.ShannonMI_Zsh;
            [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
            resp = respAsB{rr,cc}; 
            centers = centers(resp)'; PWs = PWs(resp)'; MFR = MFR(resp)'; tzMI = tzMI(resp)';
            for ii = 1:(length(bins)-1)
                binS = bins(ii); binE = bins(ii+1);
                inds = find(centers >= binS & centers < binE);
                num_cells(rr,cc,ii) = 100*length(inds)/sum(resp);
                PFW_belt(rr,cc,ii) = mean(PWs(inds));
                MFR_belt(rr,cc,ii) = mean(MFR(inds));
                zMI_belt(rr,cc,ii) = mean(tzMI(inds));
            end
        end
    end
    mzMIsC = [];
    for rr = 1:size(sRs,1)
        switch pri
            case 1
                mzMIsC(rr,:) = reshape(squeeze(PFW_belt(rr,:,:))',1,(length(bins)-1)*size(sRs,2));                
            case 2
                mzMIsC(rr,:) = reshape(squeeze(num_cells(rr,:,:))',1,(length(bins)-1)*size(sRs,2));
            case 3
                mzMIsC(rr,:) = reshape(squeeze(MFR_belt(rr,:,:))',1,(length(bins)-1)*size(sRs,2));
            case 4
                mzMIsC(rr,:) = reshape(squeeze(zMI_belt(rr,:,:))',1,(length(bins)-1)*size(sRs,2));
        end
    end
    [within,dvn,xlabels] = make_within_table({'FR','Cond','Bins'},[2,3,length(bins)-1]);
    dataT = make_between_table({mzMIsC},dvn);
    ra = RMA(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR_by_Cond_Bins','bonferroni'},[1 1 1]);
    hf = get_figure(5,[5 7 2.5 1]);
    s = generate_shades(length(bins)-1);
    tcolors = [mData.colors(1);mData.colors(1);mData.colors(1);mData.colors(2);mData.colors(2);mData.colors(2);mData.colors(3);mData.colors(3);mData.colors(3)]; tcolors = repmat(tcolors,2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',yspacing(pri),'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);

        set_axes_limits(gca,[0.35 xdata(end)+.65],[0 15]);

    format_axes(gca);
    make_bars_hollow(hbs(10:end));
    xticks = xdata; xticklabels = {'C3-B1','C3-B2','C3-B3','C4-B1','C4-B2','C4-B3','C3''-B1','C3''-B2','C3''-B3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 7 14]); xtickangle(45)
    changePosition(gca,[0.01 0.02 -0.03 -0.011])
    put_axes_labels(gca,{[],[0 0 0]},{props{pri},[0 0 0]});
    set_title(gca,'AFoR',[4.25,15],5); set_title(gca,'BFoR',[14.25,15],5)
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells_belt.pdf',fileNames{pri}),600);


    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR','bonferroni'},[1 1 1]);
    hf = get_figure(6,[11 7 1.25 1]);
    s = generate_shades(length(bins)-1);
    tcolors = mData.colors(1:3); tcolors = repmat(tcolors,2,1);%[s.m;s.m];%;s.y;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set_axes_limits(gca,[0.25 xdata(end)+.75],[0 15]);
    format_axes(gca);
%     make_bars_hollow(hbs(4:6))
    xticks = xdata; xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 15]); xtickangle(45)
    changePosition(gca,[0.2 0.02 -0.2 -0.05])
%     put_axes_labels(gca,{[],[0 0 0]},{props{pri},[0 -2 0]});
    set_title(gca,'AFoR',[1.25,15],5); set_title(gca,'BFoR',[5.25,15],5)
    save_pdf(hf,mData.pdf_folder,sprintf('%s_bar_graph_place_cells_belt_pooled.pdf',fileNames{pri}),600);
    ra.ranova
%     ra.mauchly
end

%% Percentage of Responsive Cells AsB ([AnotB BnotA])
if 1
    respT = exec_fun_on_cell_mat(respAsB,{'sum','length'}); respFA = respT.sum./respT.length;
    respORA = get_cell_list(respAsB(:,[1,2,3]),[1;2;3]); respORB = get_cell_list(respAsB(:,[4,5,6]),[1;2;3]);
    respT = exec_fun_on_cell_mat(respORA,{'sum','length'}); respF_ORA = mean(respT.sum./respT.length,2);
    respT = exec_fun_on_cell_mat(respORB,{'sum','length'}); respF_ORB = mean(respT.sum./respT.length,2);
    [within,dvn,xlabels] = make_within_table({'FR','Cond'},[2,3]);
    dataT = make_between_table({respFA*100},dvn);
    ra = RMA(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR_by_Cond','bonferroni'},[1 1 1]);
%     s = generate_shades(3); tcolors = s.g;
    tcolors = colors;%mData.colors(1:3);
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
    make_bars_hollow(hbs(4:6))
    format_axes(gca);
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[0 25])
    xticks = [xdata(1:end)]; xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    any_mean = mean(100*respF_ORA);    any_sem = std(100*respF_ORA)/sqrt(5);
    pmchar=char(177); any_text = sprintf('%.0f%c%.0f%%',any_mean,pmchar,any_sem); text(0.75,20,any_text,'FontSize',6);
    any_mean = mean(100*respF_ORB);    any_sem = std(100*respF_ORB)/sqrt(5);
    pmchar=char(177); any_text = sprintf('%.0f%c%.0f%%',any_mean,pmchar,any_sem); text(4.75,20,any_text,'FontSize',6);
    changePosition(gca,[0.17 0.02 -0.1 -0.011])
    set_title(gca,'AFoR',[1.3,25],5); set_title(gca,'BFoR',[5.3,25],5,'r');
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('percentage_PCs_responsive_AsB'),600);

    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR','bonferroni'},[1 1 1]);
%     s = generate_shades(3); tcolors = s.g;
    tcolors = {'k','r'};%mData.colors(1:3);
     hf = get_figure(500,[7 7 1.25 1]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[0 25])
    xticks = [xdata(1:end)]; xticklabels = {'AFoR','BFoR'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    
    changePosition(gca,[0.17 -0.07 -0.5 0.1])
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
%     put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned','Cells (%)'},[0 0 0]});
    format_axes(gca);
    save_pdf(hf,mData.pdf_folder,sprintf('percentage_PCs_responsive_AsB_pooled'),600);
end

%% Percentage of Responsive Cells Unique in Conditions
if 1
    percCells = [];
    cn = 1; cols = [1,2,3]; condOp = [1 -2 -3]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,1) = respT(:,1);
    cn = 2; cols = [1,2,3]; condOp = [-1 2 -3]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,2) = respT(:,1);
    cn = 3; cols = [1,2,3]; condOp = [-1 -2 3]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,3) = respT(:,1);
    cn = 1; cols = [4,5,6]; condOp = [1 -2 -3]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,4) = respT(:,1);
    cn = 2; cols = [4,5,6]; condOp = [-1 2 -3]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,5) = respT(:,1);
    cn = 3; cols = [4,5,6]; condOp = [-1 -2 3]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,6) = respT(:,1);
   
    [within,dvn,xlabels] = make_within_table({'FR','Cond'},[2,3]);
    dataT = make_between_table({percCells*100},dvn);
    ra = RMA(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR_by_Cond','bonferroni'},[1 1 1]);
    h(h==1) = 0;
%     s = generate_shades(3); tcolors = s.g;
    tcolors = mData.colors(1:3); tcolors = repmat(tcolors,1,2);
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',15,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    make_bars_hollow(hbs(4:6));
    format_axes(gca);
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[0 maxY])
    xticks = [xdata(1:end)]; xticklabels = {'C3','C4','C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
       changePosition(gca,[0.21 0.02 -0.3 -0.011])
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned','Unique Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('percentage_PCs_responsive_unique_AsB'),600);
    
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR','bonferroni'},[1 1 1]);
%     s = generate_shades(3); tcolors = s.g;
    tcolors = {'k','r'};%mData.colors(1:3);
     hf = get_figure(500,[7 7 1.25 1]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',20,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set(hbs(3),'EdgeColor','k');
%     hatch(hbs(3),45,'k','-',3,0.5); hatch(hbs(3),45,'k','-',3,0.5);
    set_axes_limits(gca,[0.25 xdata(end)+0.75],[0 maxY])
    xticks = [xdata(1:end)]; xticklabels = {'AFoR','BFoR'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.17 -0.07 -0.5 0.1])
%     changePosition(gca,[0.2 0.03 -0.4 -0.11]);
%     put_axes_labels(gca,{[],[0 0 0]},{{'Spatially Tuned','Cells (%)'},[0 0 0]});
    format_axes(gca);
    save_pdf(hf,mData.pdf_folder,sprintf('percentage_PCs_responsive_unique_AsB_pooled'),600);
end    
    
%% Place cell emergence disruption stability - big combined test
if 1
    
    percCells = [];
    indp = 0; cols = [1,2,3];
    % new cells
    cn = 2; condOp = [-1 2]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,indp+1) = respT(:,1);
    cn = 3; condOp = [-2 3]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,indp+2) = respT(:,1);
    cn = 3; condOp = [-1 3]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,indp+3) = respT(:,1);
    
    % disrupted cells
    cn = 1; condOp = [1 -2]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,indp+4) = respT(:,1);
    cn = 2; condOp = [2 -3]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,indp+5) = respT(:,1);
    cn = 1; condOp = [1 -3]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,indp+6) = respT(:,1);
    
    % remained cells
    cn = 1; condOp = [1 2]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,indp+7) = respT(:,1);
    cn = 2; condOp = [2 3]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,indp+8) = respT(:,1);
    cn = 1; cols = [1,2,3]; condOp = [1 3]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,indp+9) = respT(:,1);
    
    indp = 9; cols = [4,5,6];
    % new cells
    cn = 2; condOp = [-1 2]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,indp+1) = respT(:,1);
    cn = 3; condOp = [-2 3]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,indp+2) = respT(:,1);
    cn = 3; condOp = [-1 3]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,indp+3) = respT(:,1);
    
    % disrupted cells
    cn = 1; condOp = [1 -2]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,indp+4) = respT(:,1);
    cn = 2; condOp = [2 -3]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,indp+5) = respT(:,1);
    cn = 1; condOp = [1 -3]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,indp+6) = respT(:,1);
    
    % remained cells
    cn = 1; condOp = [1 2]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,indp+7) = respT(:,1);
    cn = 2; condOp = [2 3]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,indp+8) = respT(:,1);
    cn = 1; cols = [1,2,3]; condOp = [1 3]; respT = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),condOp),{'sum'}); respT1 = exec_fun_on_cell_mat(get_cell_list(respAsB(:,cols),cn),{'sum'}); respT = respT.sum./respT1.sum; percCells(:,indp+9) = respT(:,1);
   
 
    [within,dvn,xlabels] = make_within_table({'FR','Conds','Type'},[2,3 3]);
    dataT = make_between_table({percCells},dvn);
    ra = RMA(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'FR_by_Conds_Type','bonferroni'},[1 1 1]);
    hf = get_figure(5,[8 7 6.99 ]);
    tcolors = mData.colors(1:3); tcolors = repmat(tcolors,1,6);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);
    format_axes(gca);
    xticks = xdata; xticklabels = xlabels;
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0 0.02 -0.15 -0.05])
    put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('dynamics_cells.pdf'),600);
    ra.ranova
    ra.mauchly
end



%%