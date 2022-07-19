

%% one graph - 2x2x2
while 1
    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[10 3 1.75 1.25],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
    switch varT
        case 1 % responsive cells 
            MY = 50; ysp = 3; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'Percent of Cells'};
        case 2
            MY = 70; ysp = 3; mY = 0; titletxt = 'Response Fidelity'; ylabeltxt = {'Percent of Trials'};
        case 3
            MY = 2; ysp = 0.01; mY = 0; titletxt = 'Mutual Information'; ylabeltxt = {'Z-Score'};
        case 4
            MY = 1.7; ysp = 3; mY = 0; titletxt = 'R-squared'; ylabeltxt = {'A.U.'};
        case 7
            MY = 1.7; ysp = 0.05; mY = 1; titletxt = 'Hausdorff Frac. Dim'; ylabeltxt = {'A.U.'};
        case 10
            MY = 80; ysp = 1; mY = 0; titletxt = 'Peak Locations'; ylabeltxt = {'cm'};
        case 11
            MY = 80; ysp = 1; mY = 0; titletxt = 'Peak Locations'; ylabeltxt = {'cm'};
    end
    stp = 0.25*magfac; widths = ([1.2 1.3 1.3 1.3 1.3 0.5 0.5 0.5]+0.25)*magfac; gap = 0.16*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{''});
    tcolors = {colors{1};colors{2};colors{1};colors{2};colors{1};colors{2};colors{1};colors{2}};

    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_FoF_Cond','hsd'},[1.5 1 1]);
        xdata = make_xdata([2 2 2 2],[1 2 3]);   
%         combs = [[1:2:12]' [2:2:12]']; p = ra.MC.bonferroni.Group_by_Cond{1:2:12,6}; h = p<0.05;
    h(h==1) = 0;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
    xticklabels = {'C1','C2'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(0);
    make_bars_hollow(hbs(5:end));
    put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Air','Belt','Air','Belt'});
    ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.051 0 0]);
    format_axes_b(gca);
    set(ht,'FontWeight','Bold');
    save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
    %%
    break;
end


%% two graphs, 2x2x2 all and main effect
while 1
    ff = makeFigureRowsCols(107,[10 3 2.62 1.25],'RowsCols',[1 2],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
    switch varT
        case 1 % responsive cells 
            MY = 70; ysp = 3; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'Percent of Cells'}; % for all cells (vals) MY = 80
        case 2
            MY = 70; ysp = 4; mY = 0; titletxt = 'Response Fidelity'; ylabeltxt = {'Percent of Trials'};% for all cells (vals) MY = 70
        case 3
            MY = 0.1; ysp = 0.01; mY = -0.15; titletxt = 'Mutual Information'; ylabeltxt = {'Z-Score'};
            MY = 3.5; ysp = 0.3; mY = 0; titletxt = 'Mutual Information'; ylabeltxt = {'Z-Score'};
        case 9
            MY = 20; ysp = 1; mY = -0.15; titletxt = 'Spatial Tuning'; ylabeltxt = {'Width (cm)'};
    end
    stp = 0.28*magfac; widths = ([1.2 0.5 1.3 1.3 1.3 0.5 0.5 0.5]+0.25)*magfac; gap = 0.1*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{''});
    
    axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_FoF_Cond','hsd'},[1.5 1 1]);
        xdata = make_xdata([2 2 2 2],[1 1.5]);   
%         combs = [[1:2:12]' [2:2:12]']; p = ra.MC.bonferroni.Group_by_Cond{1:2:12,6}; h = p<0.05;
    h(h==1) = 0;
    tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
    xticklabels = {'C1','C2'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(0);
    make_bars_hollow(hbs(5:end));
    put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Air','Belt','Air','Belt'});
    ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.051 0 0]);
    format_axes_b(gca);
    set(ht,'FontWeight','Bold');
    
    axes(ff.h_axes(1,2));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'FoF','hsd'},[1.5 1 1]);
        xdata = make_xdata([2],[1 1.5]);   
    %     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
%     h(h==1) = 0;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
    xticklabels = {'Air','Belt'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(0);
    make_bars_transparent(hbs,0.5);
%     hatch(hbs,0,'w','-',2,0.25); %hatch(obj,angle,color,style,step,width)
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Pooled'});
    format_axes(gca);
    save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
    %%
    break;
end



%% one graph
while 1
    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[10 3 1.75 1.25],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
    switch varT
        case 1 % responsive cells 
            MY = 50; ysp = 3; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'Percent of Cells'};
        case 2
            MY = 70; ysp = 3; mY = 0; titletxt = 'Response Fidelity'; ylabeltxt = {'Percent of Trials'};
        case 4
            MY = 1.7; ysp = 3; mY = 0; titletxt = 'R-squared'; ylabeltxt = {'A.U.'};
        case 7
            MY = 1.7; ysp = 0.05; mY = 1; titletxt = 'Hausdorff Frac. Dim'; ylabeltxt = {'A.U.'};
        case 10
            MY = 80; ysp = 1; mY = 0; titletxt = 'Peak Locations'; ylabeltxt = {'cm'};
        case 11
            MY = 80; ysp = 1; mY = 0; titletxt = 'Peak Locations'; ylabeltxt = {'cm'};
    end
    stp = 0.25*magfac; widths = ([1.2 1.3 1.3 1.3 1.3 0.5 0.5 0.5]+0.25)*magfac; gap = 0.16*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{''});
    tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};

    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_Cond','hsd'},[1.5 1 1]);
        xdata = make_xdata([4 4],[1 1.5]);   
%         combs = [[1:2:12]' [2:2:12]']; p = ra.MC.bonferroni.Group_by_Cond{1:2:12,6}; h = p<0.05;
    h(h==1) = 0;
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
    %%
    break;
end


%% two graphs, all and cond
while 1
    ff = makeFigureRowsCols(107,[10 3 2.62 1.25],'RowsCols',[1 2],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
    switch varT
        case 1 % responsive cells 
            MY = 70; ysp = 3; mY = 0; titletxt = 'Responsivity'; ylabeltxt = {'Percent of Cells'}; % for all cells (vals) MY = 80
        case 2
            MY = 70; ysp = 4; mY = 0; titletxt = 'Response Fidelity'; ylabeltxt = {'Percent of Trials'};% for all cells (vals) MY = 70
        case 3
            MY = 0.1; ysp = 0.01; mY = -0.15; titletxt = 'Mutual Information'; ylabeltxt = {'Z-Score'};
            MY = 3.5; ysp = 0.3; mY = 0; titletxt = 'Mutual Information'; ylabeltxt = {'Z-Score'};
        case 9
            MY = 20; ysp = 1; mY = -0.15; titletxt = 'Spatial Tuning'; ylabeltxt = {'Width (cm)'};
    end
    stp = 0.28*magfac; widths = ([1.2 0.5 1.3 1.3 1.3 0.5 0.5 0.5]+0.25)*magfac; gap = 0.1*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{''});
    tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};
    axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Group_by_Cond','hsd'},[1.5 1 1]);
        xdata = make_xdata([4 4],[1 1.5]);   
    %     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
    h(h==1) = 0;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
    xticklabels = {'C1','C2','C3','C4'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(0);
    make_bars_hollow(hbs(5:end));
    put_axes_labels(gca,{'',[]},{ylabeltxt,[]});
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Control','APP'});
    ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.051 0 0]);
    set(ht,'FontWeight','Bold');
    format_axes(gca);
    
    axes(ff.h_axes(1,2));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
        xdata = make_xdata([4],[1 1.5]);   
    %     combs = [[1:2:12]' [2:2:12]']; p = ra.MC.hsd.Cond_by_CT_ET{1:2:12,6}; h = p<0.05;
%     h(h==1) = 0;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes_b(gca); xticks = xdata; 
    xticklabels = {'C1','C2','C3','C4'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(0);
    make_bars_transparent(hbs,0.5);
%     hatch(hbs,0,'w','-',2,0.25); %hatch(obj,angle,color,style,step,width)
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Pooled'});
    format_axes(gca);
    save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
    %%
    break;
end
