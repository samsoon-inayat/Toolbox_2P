function air_place_cells_figure_properties

%% load data tuned versus weakly tuned cells
while 1
    si = [Ar_t_D ArL_t_D Ars_t_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    gauss = props1.vals; 

    break;
end
disp('Done')

%% look at the distribution of peak locations for tuned and weakly tuned cells
while 1
    minBin = 0;     maxBin = 150;     incr = 50; % choosing more than 3 bins, give significant anova but not significant multcompare
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
            PLG = NaN(size(tPL)); PLnG = PLG;
            PLG(tGauss) = tPL(tGauss); 
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
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 12]); format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[-0.041 0.0 0.12 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('percent_cells_on_belt.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(1:2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.85,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 12]); format_axes(gca);
    xticks = xdata; xticklabels = {'gT','gU'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.5 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('percent_cells_on_belt_pooled_CT.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Bin','hsd'},[1 1 1]);
    xdata = make_xdata([3],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(3:5);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.4 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('percent_cells_on_belt_pooled_Bin.pdf'),600);
    
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
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    changePosition(gca,[0.05 0.0 -0.1 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('percent_cells_on_belt_pooled_CT_by_Bin.pdf'),600);
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
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.0 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{'Tuning Width (cm)'},[0 -5 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('PWs_on_belt.pdf'),600);
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
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 6]); format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[-0.041 0.0 0.12 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'zMI'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_cells_on_belt.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(1:2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.85,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 6]); format_axes(gca);
    xticks = xdata; xticklabels = {'gT','gU'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.5 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_cells_on_belt_pooled_CT.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Bin','hsd'},[1 1 1]);
    xdata = make_xdata([3],[1 2]); 
    hf = get_figure(6,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(3:5);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 6]); format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.4 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_cells_on_belt_pooled_Bin.pdf'),600);
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
respAnB = props1.good_FR_and_Gauss_loose;

for tC = 2%:4
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
    ra.ranova
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Conds','hsd'},[1 0.25 1]);
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
    yshifts = [0 0 -2 0];
    put_axes_labels(gca,{[],[0 0 0]},{sprintf('%s Cells (%%)',txtT{tC}),[0 yshifts(tC) 0]});
%     text(0.75,maxY+3,txtT{tC},'FontSize',6)
    save_pdf(hf,mData.pdf_folder,sprintf('%s_cells.pdf',txtT{tC}),600);
    break;
end
end



%% look at the distribution of peak locations for tuned and weakly tuned cells
while 1
    minBin = 0;     maxBin = 150;     incr = 50; % choosing more than 3 bins, give significant anova but not significant multcompare
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
    xdata = make_xdata([3 3 3],[1 2]); xdata(10:end) = xdata(10:end) + 1;
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
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 66]); format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[-0.041 0.0 0.12 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('percent_cells_on_belt.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Bin','hsd'},[1 1 1]);
    xdata = make_xdata([3],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(3:5);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 96]); format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.4 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('percent_cells_on_belt_pooled_Bin.pdf'),600);
    
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
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.0 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{'Tuning Width (cm)'},[0 -5 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('PWs_on_belt.pdf'),600);
    %% For zMIs of cells over belt
    [within,dvn,xlabels] = make_within_table({'Cond','Bin'},[3,3]);
    dataT = make_between_table({allzMIsA},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_Bin','hsd'},[1 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3],[1 2]); 
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
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 2.6]); format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[-0.041 0.0 0.12 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'zMI'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_cells_on_belt.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Bin','hsd'},[1 1 1]);
    xdata = make_xdata([3],[1 2]); 
    hf = get_figure(6,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(3:5);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 2.6]); format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.4 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_cells_on_belt_pooled_Bin.pdf'),600);
    %%
    break;
end




