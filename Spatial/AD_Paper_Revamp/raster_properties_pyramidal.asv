function raster_properties

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;

prop_names = {'resp','N_Resp_Trials','zMI','zMINaN','HaFD','HiFD','cells_pooled'};
event_type = {'1-D','2-D','3-D','4-D','1-T','2-T','3-T','4-T'};
sic = {[C1_t_D];[C2_t_D];[C3_t_D];[C4_t_D];[C1_i_T];[C2_i_T];[C3_i_T];[C4_i_T]};
event_type = {'1-D','2-D','3-D','4-D'};
sic = {[C1_t_D];[C2_t_D];[C3_t_D];[C4_t_D]};
pni = 7;
[all_gFR_C,all_gV_C,good_zMI_C,good_zMI_MFR_C,good_zMI_MFR_Gauss_C,nan_zMI_C,all_C] = return_values_props(oC,sic,pni);
[all_gFR_A,all_gV_A,good_zMI_A,good_zMI_MFR_A,good_zMI_MFR_Gauss_A,nan_zMI_A,all_A] = return_values_props(oA,sic,pni);
%% general for all properties including responsivity, response fidelity, zMI, Rs
while 1
    ntrials = 50; 
    si = [C1_t_D C2_t_D C3_t_D C4_t_D];
    props_C = get_props_Rs(oC.Rs(:,si),ntrials); props_A = get_props_Rs(oA.Rs(:,si),ntrials);
    sel_pop_C = props_C.vals; sel_pop_A = props_A.vals;
%     sel_pop_C = props_C.good_FR_and_tuned; sel_pop_A = props_A.good_FR_and_tuned;
%     sel_pop_C = cell_list_op(props_C.good_FR_and_tuned,props_C.good_zMI,'and'); sel_pop_A = cell_list_op(props_A.good_FR_and_tuned,props_A.good_zMI,'and');
%     sel_pop_C = cell_list_op(props_C.good_Gauss,props_C.good_zMI,'and'); sel_pop_A = cell_list_op(props_A.good_Gauss,props_A.good_zMI,'and');
    varT = 1;
    switch varT
        case 1
            var_C = sel_pop_C; mean_var_C = find_percent(var_C); 
            var_A = sel_pop_A; mean_var_A = find_percent(var_A);
        case 2
            var_C = get_vals(props_C.N_Resp_Trials,sel_pop_C); mean_var_C = exec_fun_on_cell_mat(var_C,'nanmean'); 
            var_A = get_vals(props_A.N_Resp_Trials,sel_pop_A); mean_var_A = exec_fun_on_cell_mat(var_A,'nanmean');
        case 3
            var_C = get_vals(props_C.zMI,sel_pop_C); mean_var_C = exec_fun_on_cell_mat(var_C,'nanmean');
            var_A = get_vals(props_A.zMI,sel_pop_A); mean_var_A = exec_fun_on_cell_mat(var_A,'nanmean');
        case 4
            var_C = get_vals(props_C.rs,sel_pop_C); mean_var_C = exec_fun_on_cell_mat(var_C,'nanmean'); 
            var_A = get_vals(props_A.rs,sel_pop_A); mean_var_A = exec_fun_on_cell_mat(var_A,'nanmean');
        case 5
            var_C = get_vals(props_C.nan_zMI,sel_pop_C); mean_var_C = find_percent(var_C); 
            var_A = get_vals(props_A.nan_zMI,sel_pop_A); mean_var_A = find_percent(var_A);
        case 6
            var_C = get_vals(props_C.nan_rs,sel_pop_C); mean_var_C = find_percent(var_C); 
            var_A = get_vals(props_A.nan_rs,sel_pop_A); mean_var_A = find_percent(var_A);
        case 7
            var_C = get_vals(props_C.HaFD,sel_pop_C); mean_var_C = exec_fun_on_cell_mat(var_C,'nanmean');
            var_A = get_vals(props_A.HaFD,sel_pop_A); mean_var_A = exec_fun_on_cell_mat(var_A,'nanmean');
        case 8
            var_C = get_vals(props_C.HiFD,sel_pop_C); mean_var_C = exec_fun_on_cell_mat(var_C,'nanmean');
            var_A = get_vals(props_A.HiFD,sel_pop_A); mean_var_A = exec_fun_on_cell_mat(var_A,'nanmean');
        case 9
            var_C = get_vals(props_C.PWs,sel_pop_C); mean_var_C = exec_fun_on_cell_mat(var_C,'nanmean');
            var_A = get_vals(props_A.PWs,sel_pop_A); mean_var_A = exec_fun_on_cell_mat(var_A,'nanmean');
        case 10
            var_C = get_vals(props_C.centers,sel_pop_C); mean_var_C = exec_fun_on_cell_mat(var_C,'nanmean');
            var_A = get_vals(props_A.centers,sel_pop_A); mean_var_A = exec_fun_on_cell_mat(var_A,'nanmean');
    end

    varC = mean_var_C;
    varA = mean_var_A;
    [within,dvn,xlabels] = make_within_table({'Cond'},[4]);
    dataT = make_between_table({varC;varA},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
    ra.ranova
break;
end


%% one graph
while 1
    ff = makeFigureRowsCols(107,[10 3 1.75 1.25],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
    switch varT
        case 1 % responsive cells 
            MY = 80; ysp = 3; mY = 0; titletxt = 'Responsivity (Percent of Cells)';
        case 2
            MY = 70; ysp = 3; mY = 0; titletxt = 'Response Fidelity (Percent of Trials)';
        case 4
            MY = 0.7; ysp = 3; mY = 0; titletxt = 'R-squared (Arb)';
    end
    stp = 0.25; widths = [1.2 1.3 1.3 1.3 1.3 0.5 0.5 0.5]+0.25; gap = 0.16;
    adjust_axes(ff,[mY MY],stp,widths,gap,{'R-squared'});
    tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4}};

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
    put_axes_labels(gca,{'',[0 0 0]},{{'EMM'},[0 -0.1 0]});
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Control','APP'});
    ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.051 0 0]);
    
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
            MY = 45; ysp = 3; mY = 0; titletxt = 'Responsivity (Percent of Cells)'; % for all cells (vals) MY = 80
        case 2
            MY = 90; ysp = 4; mY = 0; titletxt = 'Response Fidelity (Percent of Trials)'; % for all cells (vals) MY = 70
        case 3
            MY = 3.9; ysp = 0.3; mY = 0; titletxt = 'Mutual Information (Z-Score)';
    end
    stp = 0.25; widths = [1.2 0.5 1.3 1.3 1.3 0.5 0.5 0.5]+0.25; gap = 0.1;
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
    put_axes_labels(gca,{'',[0 0 0]},{{'EMM'},[0 -0.1 0]});
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Control','APP'});
    ht = set_axes_top_text_no_line(gcf,gca,titletxt,[0 -0.051 0 0]);
    set(ht,'FontWeight','Bold');
    
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
    make_bars_hollow(hbs(5:end));
    set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,4,{'Pooled'});
    save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
    %%
    break;
end

