
%% compare percent responsive cells pooled air on or off brake vs no-brake exc and inh
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    event_type = {'Air ON','Air OFF','Light ON'};
%     sic = {[Lb Lbs Ab_On Abs_On];[Ab_Off Abs_Off];ArL_L;[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off]};%;M_On;M_Off};
    allsic = {{[Ab_On Abs_On];[Ar_On ArL_On Ars_On]};
    };
%     allsic = {{[Ab_Off Abs_Off];[Ar_Off ArL_Off Ars_Off]}};
%     ind = 3;
%     sic = allsic{1:6};
    all_resp = []; all_resp_exc = []; all_resp_inh = [];
    for jj = 1:length(allsic)
        sic = allsic{jj};
        clear all_gFR all_exc all_inh
        for ii = 1:length(sic)
            sit = sic{ii};
            tRs = o.Rs(:,sit);
            props1 = get_props_Rs(tRs,ntrials);
            gFR = props1.vals;
            gFR = cell_list_op(gFR,[],'or',1);
            all_gFR(:,ii) = gFR;
            gFR = props1.exc;
            gFR = cell_list_op(gFR,[],'or',1);
            all_exc(:,ii) = gFR;
            gFR = props1.inh;
            gFR = cell_list_op(gFR,[],'or',1);
            all_inh(:,ii) = gFR;
        end
        all_resp = [all_resp all_gFR]; all_resp_exc = [all_resp_exc all_exc]; all_resp_inh = [all_resp_inh all_inh];
    end
    good_FR = [all_resp_exc all_resp_inh];
    good_FR_any = cell_list_op(good_FR,[],'or',1);
    good_FR_all = cell_list_op(good_FR,[],'and',1);
    per_active =[]; 
    for rr = 1:size(good_FR,1)
        for cc = 1:size(good_FR,2)
            tts = good_FR{rr,cc};
            per_active(rr,cc) = 100*sum(tts)/length(tts);
        end
        tts = good_FR_any{rr,1};        per_active_any(rr) = 100*sum(tts)/length(tts);
        tts = good_FR_all{rr,1};        per_active_all(rr) = 100*sum(tts)/length(tts);
    end
    [within,dvn,xlabels] = make_within_table({'Type','Cond'},[2,2]);
    dataT = make_between_table({per_active},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
  %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Type_by_Cond','hsd'},[1.5 1 1]);
    h(h==1) = 0;
	xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors(3:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 30]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-Exc','NB-Exc','B-Inh','NB-Inh'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 15 30]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.5 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{'Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all1off.pdf',ntrials),600);
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
    mVar(3) = NaN; mVar(4) = NaN; semVar(3:4) = NaN;
	xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 30]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-Exc','NB-Exc','B-Inh','NB-Inh'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 15 30]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.5 -0.04]); %put_axes_labels(gca,{[],[0 0 0]},{{'Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all2off.pdf',ntrials),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Type','hsd'},[1.5 1 1]);
    mVar(3) = NaN; mVar(4) = NaN; semVar(3:4) = NaN;
%     h(h==1) = 0;
	xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors(7:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 30]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-Exc','NB-Exc','B-Inh','NB-Inh'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 15 30]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.5 -0.04]); %put_axes_labels(gca,{[],[0 0 0]},{{'Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all3off.pdf',ntrials),600);
    %%
    break;
end

%% compare percent response fidelity pooled air on or off brake vs no-brake
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    event_type = {'Air ON','Air OFF','Light ON'};
%     sic = {[Lb Lbs Ab_On Abs_On];[Ab_Off Abs_Off];ArL_L;[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off]};%;M_On;M_Off};
    allsic = {{[Ab_On Abs_On];[Ar_On ArL_On Ars_On]};
    };
%     allsic = {{[Ab_Off Abs_Off];[Ar_Off ArL_Off Ars_Off]};
%         };
%     ind = 3;
%     sic = allsic{1:6};
    all_resp = []; all_resp_exc = []; all_resp_inh = [];
    for jj = 1:length(allsic)
        sic = allsic{jj};
        clear all_gFR all_exc all_inh
        for ii = 1:length(sic)
            sit = sic{ii};
            tRs = o.Rs(:,sit);
            props1 = get_props_Rs(tRs,ntrials);
            gFR = props1.vals; rf = props1.N_Resp_Trials;
            all_gFR(:,ii) = mean(exec_fun_on_cell_mat(rf,'mean',gFR),2);
            gFR = props1.exc; rf = props1.N_Resp_Trials;
            all_exc(:,ii) = mean(exec_fun_on_cell_mat(rf,'mean',gFR),2);
            gFR = props1.inh;rf = props1.N_Resp_Trials;
            all_inh(:,ii) = mean(exec_fun_on_cell_mat(rf,'mean',gFR),2);
        end
        all_resp = [all_resp all_gFR]; all_resp_exc = [all_resp_exc all_exc]; all_resp_inh = [all_resp_inh all_inh];
    end
    good_FR = [all_resp_exc all_resp_inh];

    [within,dvn,xlabels] = make_within_table({'Type','Cond'},[2,2]);
    dataT = make_between_table({good_FR},dvn);
    ra = RMA(dataT,within);
    ra.ranova
  
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Type_by_Cond','hsd'},[1.5 1 1]);
	xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors(3:6);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',20,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 90]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-Exc','NB-Exc','B-Inh','NB-Inh'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 35 70]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.5 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{'Trials (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_allf1.pdf',ntrials),600);
     %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
    mVar(3) = NaN; mVar(4) = NaN; semVar(3:4) = NaN;
	xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',20,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 90]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-Exc','NB-Exc','B-Inh','NB-Inh'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 35 70]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.5 -0.04]); %put_axes_labels(gca,{[],[0 0 0]},{{'Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_allf2.pdf',ntrials),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Type','hsd'},[1.5 1 1]);
    mVar(3) = NaN; mVar(4) = NaN; semVar(3:4) = NaN;
%     h(h==1) = 0;
	xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors(7:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',20,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 90]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-Exc','NB-Exc','B-Inh','NB-Inh'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 35 70]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.5 -0.04]); %put_axes_labels(gca,{[],[0 0 0]},{{'Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_allf3.pdf',ntrials),600);
    %%
    break;
end
