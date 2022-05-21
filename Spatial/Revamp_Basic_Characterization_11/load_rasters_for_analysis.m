function load_rasters_for_analysis

%% Load Data
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei = evalin('base','ei'); 
    
    selContexts = [3 4 5];
    rasterNames = {'air77T','air77T','air77T'};
    o = get_data(ei,selContexts,rasterNames);
    break
end
n = 0;
%%
%% Speed Figure
while 1
    Rs = o.Rs;
    for ii = 1:length(ei)
        b1 = ei{ii}.b;
        for jj = 1:10
            alds(ii,jj) = b1.dist(b1.stim_r(jj+10)) - b1.dist(b1.air_puff_r(jj+20));
        end
    end
    ald = round(mean(alds(:)));
     ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 6],'spaceRowsCols',[0 0],'rightUpShifts',[-0.05 0.13],'widthHeightAdjustment',[-3 -350]);
    set(gcf,'color','w'); set(gcf,'Position',[5 5 5.60 1]);
    [Y,E] = discretize(1:49,3);
    all_speeds = []; cTxt = {'Ar-Trials','Ar-Inter-Trials','ArL-Trials','ArL-Inter-Trials','Ar*-Trials','Ar*-Inter-Trials'}; 
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
        xs = Rs{1,cn}.xs; N = length(xs);
        mspeed = mean(mean_speed_over_trials(:,1:N)); semspeed = std(mean_speed_over_trials(:,1:N))/sqrt(5);
        plot(xs,mspeed);
        shadedErrorBar(xs,mspeed,semspeed);
        changePosition(gca,[0.1 0.15 -0.05 -0.15]);
        if cn == 1
            put_axes_labels(gca,{'',[0 0 0]},{{'Speed (cm/sec)'},[0 -1 0]});
        end
        if ismember(cn,[1 3 5])
            if cn == 3
                plot([ald ald],[0 30],'m:','linewidth',0.5);
            end
            plot([50 50],[0 30],'k--','linewidth',0.25);
            plot([100 100],[0 30],'k--','linewidth',0.25);
            bTxt = {'dB1','dB2','dB3'}; 
            xbTxt = [25 75 125]-7; ybTxt = 31;
            for ii = 1:length(bTxt)
                text(xbTxt(ii),ybTxt,bTxt{ii},'FontSize',5);
            end
            xlabel('Distance (cm)');
        else
            plot([5 5],[0 30],'k--','linewidth',0.25);
            plot([10 10],[0 30],'k--','linewidth',0.25);
            bTxt = {'tB1','tB2','tB3'}; xbTxt = [2.5 7.5 12.5]-1; ybTxt = 31;
            for ii = 1:length(bTxt)
                text(xbTxt(ii),ybTxt,bTxt{ii},'FontSize',5);
            end
            xlabel('Time (sec)');
        end
        text(xbTxt(1),ybTxt+10,cTxt{cn},'FontSize',5);
        ylim([0 30]);
        box off;
        format_axes(gca);
    end
    
    
    
    save_pdf(ff.hf,mData.pdf_folder,sprintf('speeds_345'),600);
%%
    all_speeds_Trials = all_speeds(:,[1:3 7:9 13:15]+0);
    [within,dvn,xlabels] = make_within_table({'Conds','Bins'},[3,3]);
    dataT = make_between_table({all_speeds_Trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Conds_by_Bins','hsd'},[1 0.5 1]);
    h(h==1) = 0;
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors(7:9); tcolors = repmat(tcolors,3,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 37]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[-0.01 0.02 0.05 -0.05])
    put_axes_labels(gca,{[],[0 0 0]},{'Avg. Speed (cm/sec)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('avg_speed_anova.pdf'),600);
    ra.ranova
    ra.mauchly
    %%
    all_speeds_Trials = all_speeds(:,[1:3 7:9 13:15]+0);
    [within,dvn,xlabels] = make_within_table({'Conds','Bins'},[3,3]);
    dataT = make_between_table({all_speeds_Trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Bins','hsd'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors(7:9);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 37]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'yticklabels',[]); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.35 -0.05])
%     put_axes_labels(gca,{[],[0 0 0]},{'Avg. Speed (cm/sec)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('avg_speed_anova_pooled.pdf'),600);
    %%
    all_speeds_Trials = all_speeds(:,[7:9]+0);
    [within,dvn,xlabels] = make_within_table({'Bins'},[3]);
    dataT = make_between_table({all_speeds_Trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Bins','hsd'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors(7:9);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 37]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'yticklabels',[]); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.35 -0.05])
%     put_axes_labels(gca,{[],[0 0 0]},{'Avg. Speed (cm/sec)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('avg_speed_anova_pooled.pdf'),600);
    %%
    all_speeds_Trials = all_speeds(:,[1:3 7:9 13:15]+3);
    [within,dvn,xlabels] = make_within_table({'Conds','Bins'},[3,3]);
    dataT = make_between_table({all_speeds_Trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Bins','hsd'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors([10 11 12]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 37]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'tB1','tB2','tB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.35 -0.05])
    put_axes_labels(gca,{[],[0 0 0]},{'Avg. Speed (cm/sec)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('avg_speed_anova_pooled_IT.pdf'),600);
    %%
    break;
end