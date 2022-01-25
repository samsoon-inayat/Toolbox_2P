function air_time_cells_figure_properties

%% load data tuned versus weakly tuned cells
while 1
    si = [Ar_i_T ArL_i_T Ars_i_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    gauss = props1.good_FR_and_Gauss_loose; n_gauss = props1.good_FR_and_notGauss_loose;
%     gauss = cell_list_op(props1.good_FR,gFR_D_g_T,'and'); n_gauss = cell_list_op(props1.good_FR,gFR_T_g_D,'and');
    break;
end
disp('Done')

%% look at the distribution of peak locations for tuned and weakly tuned cells
while 1
    minBin = 0;     maxBin = 15;     incr = 5; % choosing more than 3 bins, give significant anova but not significant multcompare
    % three bins also make sense because LED in condition 2 comes ON at
    % 110cm
    bins = minBin:incr:maxBin;
    allG = []; allnG = []; allPWs = []; allzMIsG = []; allzMIsnG = [];
    for rr = 1:size(gauss,1)
        pG = []; pnG = []; pGPWs = []; pGzMIs = []; pnGzMIs = [];
        for cc = 1:size(gauss,2)
            R = Rs{rr,cc};
            tPL = props1.peak_locations{rr,cc};            tzMIs = props1.zMI{rr,cc};
            tGauss = gauss{rr,cc};
            tnGauss = n_gauss{rr,cc};
            PLG = NaN(size(tPL)); PLnG = PLG;
            PLG(tGauss) = tPL(tGauss); PLnG(tnGauss) = tPL(tnGauss);
            [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);             
            bData = []; bData1 = []; bDataPWs = []; bDatazMIs = []; bDatazMIs1 = [];
            for ii = 1:(length(bins)-1)
                binS = bins(ii); binE = bins(ii+1);
                inds = PLG >= binS & PLG < binE;
                bData(ii) = 100*sum(inds)/length(tPL);
                bDataPWs(ii) = nanmean(PWs(inds));
                bDatazMIs(ii) = nanmean(tzMIs(inds));
                inds = PLnG >= binS & PLnG < binE;
                bData1(ii) = 100*sum(inds)/length(tPL);
                bDatazMIs1(ii) = nanmean(tzMIs(inds));
            end
            pG = [pG bData]; pnG = [pnG bData1]; pGPWs = [pGPWs bDataPWs];  pGzMIs = [pGzMIs bDatazMIs]; pnGzMIs = [pnGzMIs bDatazMIs1];
        end
        allG(rr,:) = pG; allnG(rr,:) = pnG; allPWs(rr,:) = pGPWs; allzMIsG(rr,:) = pGzMIs; allzMIsnG(rr,:) = pnGzMIs;
    end
    disp('Done')
    %% For percentage of cells over belt
    [within,dvn,xlabels] = make_within_table({'CT','Cond','Bin'},[2,3,3]);
    dataT = make_between_table({allG,allnG},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT_Cond_Bin','hsd'},[1 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3 3 3 3],[1 2]); xdata(10:end) = xdata(10:end) + 1;
    hf = get_figure(5,[8 7 3 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.colors(1:9),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'tB1','tB2','tB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[-0.041 0.0 0.12 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('percent_cells_on_IT.pdf'),600);
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(1:2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'gTt','gUt'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.5 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('percent_cells_on_IT_pooled_CT.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT_by_Bin','hsd'},[1 1 1]);
    xdata = make_xdata([3 3],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.dcolors(3:5),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     maxY = maxY;
    make_bars_hollow(hbs(4:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 12]); format_axes(gca);
    xticks = xdata; xticklabels = {'tB1','tB2','tB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    changePosition(gca,[0.05 0.0 -0.1 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
    %% For place widths over belt of tuned cells
    [within,dvn,xlabels] = make_within_table({'Cond','Bin'},[3,3]);
    dataT = make_between_table({allPWs},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_Bin','hsd'},[1 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3],[1 2]);
    hf = get_figure(5,[8 7 1.7 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'tB1','tB2','tB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.0 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{'Tuning Width (sec)'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('PWs_on_IT.pdf'),600);
    %% For zMIs of cells over belt
    [within,dvn,xlabels] = make_within_table({'CT','Cond','Bin'},[2,3,3]);
    dataT = make_between_table({allzMIsG,allzMIsnG},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT_Cond_Bin','hsd'},[1 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3 3 3 3],[1 2]); xdata(10:end) = xdata(10:end) + 1;
    hf = get_figure(5,[8 7 3 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.colors(1:9),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 1.5]); format_axes(gca);
    xticks = xdata; xticklabels = {'tB1','tB2','tB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[-0.041 0.0 0.12 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'zMI'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_cells_on_IT.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(1:2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 1.5]); format_axes(gca);
    xticks = xdata; xticklabels = {'gTt','gUt'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.5 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_cells_on_IT_pooled_CT.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Bin','hsd'},[1 1 1]);
    xdata = make_xdata([3],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(3:5);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.25,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 1.5]); format_axes(gca);
    xticks = xdata; xticklabels = {'tB1','tB2','tB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.4 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_cells_on_IT_pooled_Bin.pdf'),600);
%%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT_by_Bin','hsd'},[1 1 1]);
    xdata = make_xdata([3 3],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.dcolors(3:5),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     maxY = maxY;
    make_bars_hollow(hbs(4:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    changePosition(gca,[0.05 0.0 -0.1 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
    %%
    break;
end


%% trial to trial peak location difference
while 1
    for rr = 1:size(gauss,1)
        dPLG = []; dPLnG = [];
        for cc = 1:size(gauss,2)
            R = Rs{rr,cc};
            tPL = props1.peak_locations_trials{rr,cc};
            tGauss = gauss{rr,cc};
            tnGauss = n_gauss{rr,cc};
            PLG = tPL(tGauss,:); dPLG = [dPLG mean(diff(PLG,[],2))];
            PLnG = tPL(tnGauss,:); dPLnG = [dPLnG mean(diff(PLnG,[],2))];
        end
        alldPLG(rr,:) = dPLG;
        alldPLnG(rr,:) = dPLnG;
    end
    disp('Done')
    %% For anovarm
    [within,dvn,xlabels] = make_within_table({'CT','Cond','TP'},[2,3,9]);
    dataT = make_between_table({alldPLG,alldPLnG},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT_Cond_TP','hsd'},[1 1 1]);
    xdata = make_xdata([3],[1]);
    hf = get_figure(5,[8 7 1.5 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.5 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('percent_cells_on_belt.pdf'),600);
    %%
    break;
end

%% Place cell emergence disruption stability
props1 = get_props_Rs(Rs,50);
respAnB = props1.good_FR;

for tC = 4%:4
while 1
    txtT = {'Unique','New','Disrupted','Common'};
    or_cells = get_cell_list(respAnB,[1;2;3]);
    percCells = []; 
    switch tC
        case 1 % unique place cells
            respCells1 = get_cell_list(respAnB,[1 -2 -3]);    respCells2 = get_cell_list(respAnB,[-1 2 -3]);    respCells3 = get_cell_list(respAnB,[-1 -2 3]);
        case 2 % new place cells
            respCells1 = get_cell_list(respAnB,[-1 2]);    respCells2 = get_cell_list(respAnB,[-2 3]);    respCells3 = get_cell_list(respAnB,[-1 3]);
        case 3 % disrupted place cells
            respCells1 = get_cell_list(respAnB,[1 -2]);    respCells2 = get_cell_list(respAnB,[2 -3]);    respCells3 = get_cell_list(respAnB,[1 -3]);
        case 4 % remained place cells
            respCells1 = get_cell_list(respAnB,[1 2]);    respCells2 = get_cell_list(respAnB,[2 3]);    respCells3 = get_cell_list(respAnB,[1 3]);
    end
    
    respCells = [respCells1(:,1) respCells2(:,1) respCells3(:,1)];
    percCells = 100*exec_fun_on_cell_mat(respCells,'sum')./exec_fun_on_cell_mat(or_cells,'length');

    [within,dvn,xlabels] = make_within_table({'Conds'},[3]);
    if tC > 1
       xlabels = {'Ar-ArL','ArL-Ar*','Ar-Ar*'};
    else
        xlabels = {'Ar','ArL','Ar*'};
    end
    dataT = make_between_table({percCells},dvn);
    ra = RMA(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Conds','hsd'},[1 0.25 1]);
    [om,osem] = findMeanAndStandardError(mVar)
    hf = get_figure(5,[8 7 1.25 1]);
    if tC > 1
        tcolors = mData.dcolors(6:8);%[s.m;s.c;s.y];
    else
        tcolors = mData.colors(6:8);%[s.m;s.c;s.y];
    end
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0 maxY]);
    format_axes(gca);
    xticks = xdata; xticklabels = xlabels;
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.1 0.02 -0.45 -0.05])
    yshifts = [0 0 -2 0]
    put_axes_labels(gca,{[],[0 0 0]},{sprintf('%s Cells (%%)',txtT{tC}),[0 yshifts(tC) 0]});
%     text(0.75,maxY+3,txtT{tC},'FontSize',6)
    save_pdf(hf,mData.pdf_folder,sprintf('%s_cells.pdf',txtT{tC}),600);
    ra.ranova
%     ra.mauchly
    break;
end
end



%% look at the distribution of peak locations for All Resp > 50% cells
while 1
    minBin = 0;     maxBin = 15;     incr = 5; % choosing more than 3 bins, give significant anova but not significant multcompare
    % three bins also make sense because LED in condition 2 comes ON at
    % 110cm
    bins = minBin:incr:maxBin;
    allG = []; allnG = []; allPWs = []; allzMIsG = []; allzMIsnG = [];
    allA = []; allzMIsA = [];
    for rr = 1:size(gauss,1)
        pG = []; pnG = []; pGPWs = []; pGzMIs = []; pnGzMIs = []; pA = []; pAzMIs = [];
        for cc = 1:size(gauss,2)
            R = Rs{rr,cc};
            tPL = props1.peak_locations{rr,cc};            tzMIs = props1.zMI{rr,cc};
            tGauss = gauss{rr,cc};
            tnGauss = n_gauss{rr,cc};
            PLG = NaN(size(tPL)); PLnG = PLG;
            PLG(tGauss) = tPL(tGauss); PLnG(tnGauss) = tPL(tnGauss);
            [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);             
            bData = []; bData1 = []; bDataPWs = []; bDatazMIs = []; bDatazMIs1 = []; bData2 = []; bDatazMIs2 = []; 
            for ii = 1:(length(bins)-1)
                binS = bins(ii); binE = bins(ii+1);
                inds = PLG >= binS & PLG < binE;
                bData(ii) = 100*sum(inds)/length(tPL);
                bDataPWs(ii) = nanmean(PWs(inds));
                bDatazMIs(ii) = nanmean(tzMIs(inds));
                inds = PLnG >= binS & PLnG < binE;
                bData1(ii) = 100*sum(inds)/length(tPL);
                bDatazMIs1(ii) = nanmean(tzMIs(inds));
                
                inds = tPL >= binS & tPL < binE;
                bData2(ii) = 100*sum(inds)/length(tPL);
                bDatazMIs2(ii) = nanmean(tzMIs(inds));
            end
            pG = [pG bData]; pnG = [pnG bData1]; pGPWs = [pGPWs bDataPWs];  pGzMIs = [pGzMIs bDatazMIs]; pnGzMIs = [pnGzMIs bDatazMIs1];
            pA = [pA bData2]; pAzMIs = [pAzMIs bDatazMIs2];
        end
        allG(rr,:) = pG; allnG(rr,:) = pnG; allPWs(rr,:) = pGPWs; allzMIsG(rr,:) = pGzMIs; allzMIsnG(rr,:) = pnGzMIs;
        allA(rr,:) = pA; allzMIsA(rr,:) = pAzMIs;
    end
    disp('Done')
    %% For percentage of cells over belt
    [within,dvn,xlabels] = make_within_table({'Cond','Bin'},[3,3]);
    dataT = make_between_table({allG},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_Bin','hsd'},[1 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3],[1 2]); %xdata(10:end) = xdata(10:end) + 1;
    hf = get_figure(5,[8 7 3 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.colors(1:9),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 12]); format_axes(gca);
    xticks = xdata; xticklabels = {'tB1','tB2','tB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[-0.041 0.0 0.12 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('percent_cells_on_IT.pdf'),600);
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(1:2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'gTt','gUt'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.5 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('percent_cells_on_IT_pooled_CT.pdf'),600);
    
    %% For place widths over belt of tuned cells
    [within,dvn,xlabels] = make_within_table({'Cond','Bin'},[3,3]);
    dataT = make_between_table({allPWs},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_Bin','hsd'},[1 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3],[1 2]);
    hf = get_figure(5,[8 7 1.7 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'tB1','tB2','tB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.0 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{'Tuning Width (sec)'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('PWs_on_IT.pdf'),600);
    %% For zMIs of cells over belt
    [within,dvn,xlabels] = make_within_table({'CT','Cond','Bin'},[2,3,3]);
    dataT = make_between_table({allzMIsG,allzMIsnG},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT_Cond_Bin','hsd'},[1 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3 3 3 3],[1 2]); xdata(10:end) = xdata(10:end) + 1;
    hf = get_figure(5,[8 7 3 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.colors(1:9),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 1.5]); format_axes(gca);
    xticks = xdata; xticklabels = {'tB1','tB2','tB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[-0.041 0.0 0.12 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'zMI'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_cells_on_IT.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(1:2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 1.5]); format_axes(gca);
    xticks = xdata; xticklabels = {'gTt','gUt'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.5 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_cells_on_IT_pooled_CT.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Bin','hsd'},[1 1 1]);
    xdata = make_xdata([3],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(3:5);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.25,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 1.5]); format_axes(gca);
    xticks = xdata; xticklabels = {'tB1','tB2','tB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.4 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_cells_on_IT_pooled_Bin.pdf'),600);
    %%
    break;
end




%% Overlap Indices ImageSC
while 1
    ntrials = 50;
    si = [Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
    props1 = get_props_Rs(o.Rs(:,si),[50,100]);
%     props1 = get_props_Rs(o.Rs(:,si),[40,70]);
%     props1 = get_props_Rs(o.Rs(:,si),[10,40]);
    resp = [props1.good_FR];% resp_speed];
%     resp(:,3:5) = props1.good_FR_T1(:,3:5);
    resp(:,6:8) = cell_list_op(props1.good_FR_IT1(:,6:8),cell_list_op(props1.good_FR_IT(:,6:8),[],'not'),'and');
%     resp = [resp props1.good_FR_IT1(:,6:8)];
%     resp = [resp(:,1:8) props1.good_FR_and_tuned(:,9:11)];
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = rasterNamesTxt(si); 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.15 1.15]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
%     set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
%     text(-0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); 
    set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.15 0.05 0.08 0.09]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',ntrials),600);
    %%
    break;
end

%% agglomerative hierarchical clustering
while 1
    mOI1 = mOI;
    mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1);
    tree = linkage(mOI1,'average','euclidean');
%     tree = linkage(Di,'average');
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 1.35 0.95]);
    set(H,'linewidth',1);
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel('Eucledian Distance');changePosition(hx,[0 -0.05 0]);
    changePosition(gca,[0.07 0.0 0.05 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);
    %%
    break;
end
