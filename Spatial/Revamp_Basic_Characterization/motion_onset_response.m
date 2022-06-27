function firing_rate_motion_vs_rest

%% Motion onset and offset events
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    sic = {[M_On];[M_Off]};
    clear all_gFR all_exc all_inh all_gV 
    prop_names = {'resp','N_Resp_Trials','zMI','zMINaN','HaFD','HiFD','cells_pooled'};
    cell_sel = {'good_FR','vals','exc','inh'};
    varName = {'all_gFR','all_gV','all_exc','all_inh'};
    for cii = 1:length(cell_sel);
        cmdTxt = sprintf('clear %s',varName{cii});eval(cmdTxt);
    end
    pni = 7;
    all_exc_inh = [];
    for ii = 1:length(sic)
        sit = sic{ii};
        tRs = o.Rs(:,sit);
        props1 = get_props_Rs(tRs,ntrials);
        for cii = 1:length(cell_sel)
            cmdTxt = sprintf('gFR = props1.%s;',cell_sel{cii});eval(cmdTxt);
            if pni == 1
                rf = find_percent(gFR);
                cmdTxt = sprintf('%s(:,ii) = mean(rf,2);',varName{cii}); eval(cmdTxt)
            else
                if pni == 7
                    cmdTxt = sprintf('%s(:,ii) = cell_list_op(gFR,[],''or'',1);',varName{cii});eval(cmdTxt);
                else
                    if pni == 3 || pni == 4
                        cmdtxt = sprintf('rf = props1.%s;',prop_names{3}); eval(cmdtxt);
                        if pni == 3
                            cmdTxt = sprintf('%s(:,ii) = mean(exec_fun_on_cell_mat(rf,''nanmean'',gFR),2);',varName{cii});eval(cmdTxt);
                        else
                            for rrr = 1:size(rf,1)
                                for ccc = 1:size(rf,2)
                                    temp = isnan(rf{rrr,ccc});
                                    rf1(rrr,ccc) = 100*sum(temp(gFR{rrr,ccc}))/size(rf{rrr,ccc},1);
                                end
                            end
                            cmdTxt = sprintf('%s(:,ii) = mean(rf1,2);',varName{cii}); eval(cmdTxt)
                        end
                    else
                        cmdtxt = sprintf('rf = props1.%s;',prop_names{pni}); eval(cmdtxt);
                        cmdTxt = sprintf('%s(:,ii) = mean(exec_fun_on_cell_mat(rf,''mean'',gFR),2);',varName{cii});eval(cmdTxt);
                    end
                end
            end
        end
        all_exc_inh = [all_exc_inh all_exc(:,ii) all_inh(:,ii)];
    end
    
     %
%     avar = find_percent(all_exc_inh);
    avar = (all_exc_inh);
    [within,dvn,xlabels,withinD] = make_within_table({'ET','CT'},[2,2]); withinD3 = withinD;
    dataT = make_between_table({avar},dvn);
    ra = RMA(dataT,within,{0.05,{'lsd','hsd','bonferroni'}});
    ra.ranova
    print_for_manuscript(ra)
    %%
    avar1 = (all_gV);
    [within,dvn,xlabels,withinD] = make_within_table({'ET'},[2]); withinD3 = withinD;
    dataT = make_between_table({avar1},dvn);
    ra_cond_et = RMA(dataT,within,{0.05,{'lsd','hsd','bonferroni'}});
    ra_cond_et.ranova
    print_for_manuscript(ra_cond_et)
    %%
    break;
end


%% exc and inh  percentage of cells and zMI
while 1
    clc
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 1],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',...
        [10 -410]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 1.25 1.25]);
    MY = 18; ysp = 1; sigLinesStartYFactor = 1.5; mY = 0; % responsive cells
%     MY = 70; ysp = 5; sigLinesStartYFactor = 1.5; mY = 0; % responsive fidelity
%     MY = 1; ysp = 0.15; sigLinesStartYFactor = 0.1; mY = -0.35; % responsive fidelity
    stp = 0.43; widths = [0.55 0.4 0.4 0.4 0.4 0.4]+0.1; gap = 0.09;
    set(ff.hf,'Units','inches');
    for aii = 1:length(ff.h_axes)
        sel_ax = ff.h_axes(1,aii);
        set(sel_ax,'Units','inches');
        axesPosd = get(sel_ax,'Position');
        if aii == 1
            rt = stp; 
            ylabel(sel_ax,{'Estimated','Marginal Means'});
        else
            rt = rt + widths(aii-1) + gap; set(sel_ax,'yticklabels',[]);
        end
        set(sel_ax,'Position',[rt axesPosd(2) widths(aii) axesPosd(4)]);
        ylim(sel_ax,[mY MY]);
    end

    %++++++++++++++++++++++++
    axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'ET_by_CT','bonferroni'},[1.5 1 1]); 
	xdata = make_xdata([2 2],[1 1.5]);    tcolors = repmat(mData.dcolors([3 5]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'Exc','Inh'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,2,{'MOn','MOff'});
%     %+++++++++++++++++++++++
%      axes(ff.h_axes(1,2));
%     [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra_cond_et,{'ET','bonferroni'},[1.5 1 1]); 
% 	xdata = make_xdata([4],[1 1.5]);    tcolors = repmat(mData.colors([12 4 9]),1,4);
%     [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
%         'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
%         'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
%     set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'AOn','AOff','Arb'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
%     set_axes_top_text(ff.hf,ff.h_axes(1:2),'All Responsive Cells');
 
     save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs_3_factors.pdf',ntrials),600);
    %%
    break;
end 

%% all exc and inh (all_gV variable) percentage of cells, response fidelity, and zMI
while 1
    clc
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 2],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',...
        [10 -410]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 1.55 1.25]);
    MY = 50; ysp = 5; sigLinesStartYFactor = 1.5; mY = 0; % responsive fidelity
%     MY = 1.5; ysp = 0.15; sigLinesStartYFactor = 0.1; mY = -0.25; % responsive fidelity
    stp = 0.43; widths = [0.55 0.25 0.4 0.4 0.4 0.4]+0.1; gap = 0.09;
    set(ff.hf,'Units','inches');
    for aii = 1:length(ff.h_axes)
        sel_ax = ff.h_axes(1,aii);
        set(sel_ax,'Units','inches');
        axesPosd = get(sel_ax,'Position');
        if aii == 1
            rt = stp; 
            ylabel(sel_ax,{'Estimated','Marginal Means'});
        else
            rt = rt + widths(aii-1) + gap; set(sel_ax,'yticklabels',[]);
        end
        set(sel_ax,'Position',[rt axesPosd(2) widths(aii) axesPosd(4)]);
        ylim(sel_ax,[mY MY]);
    end

    %++++++++++++++++++++++++
    axes(ff.h_axes(1,1));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'ET_by_CT','bonferroni'},[1.5 1 1]); 
	xdata = make_xdata([2 2],[1 1.5]);    tcolors = repmat(mData.dcolors([3 5]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,[],[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'Exc','Inh'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,2,{'MOn','MOff'});
    %+++++++++++++++++++++++
     axes(ff.h_axes(1,2));
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'ET','bonferroni'},[1.5 1 1]); 
	xdata = make_xdata([2],[1 1.5]);    tcolors = repmat(mData.dcolors([6 7]),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'MOn','MOff'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,2,{'Pooled'});
%     set_axes_top_text(ff.hf,ff.h_axes(1:2),'All Responsive Cells');
 
     save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs_3_factors.pdf',ntrials),600);
    %%
    break;
end 


%% Brake vs No-Brake and Voluntary Motion On Off
% three circle Venn diagram, percentages of complementary and conjunctive cells compared among Lb/Lbs, Ar_L/Ars_L, and ArL_L
while 1
    all_resp = [];
    ntrials = 50;
    sic = {[Ab_On Abs_On Ab_Off Abs_Off Ab_Offc Abs_Offc Lb Lbs];[Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off Ar_Offc ArL_Offc Ars_Offc Ar_L ArL_L Ars_L];[M_On M_Off]};
    clear all_gFR
    for ii = 1:length(sic)
        sit = sic{ii};
        tRs = o.Rs(:,sit);
        props1 = get_props_Rs(tRs,ntrials);
        gFR = props1.vals;
        gFR = cell_list_op(gFR,[],'or',1);
        all_gFR(:,ii) = gFR;
    end
    all_resp = [all_resp all_gFR];
    good_FR = all_gFR;
    per_unique =[];
    condMat = [];
    for rr = 1:size(good_FR,2)
        for cc = 1:size(good_FR,2)
            if rr == cc
                condMat(cc) = cc;
            else
                condMat(cc) = -cc;
            end
        end
        temp_unique = get_cell_list(good_FR,condMat,1);
        unique_gFR(:,rr) = get_cell_list(good_FR,condMat,0);
        per_unique(:,rr) = temp_unique(:,1);
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},size(good_FR,2));
    dataT = make_between_table({per_unique},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    
    any_gFR = cell_list_op(all_gFR,[],'or',1);
    pangfr = 100*exec_fun_on_cell_mat(any_gFR,'sum')./exec_fun_on_cell_mat(any_gFR,'length');
    [mpangfr,sempangfr] = findMeanAndStandardError(pangfr);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = make_xdata([size(good_FR,2)],[1]);
%     h(h==1) = 0;
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
    set(gcf,'color','w');    set(gcf,'Position',[10 3 1.6 1.25]);
    MY = 39; ysp = 1.5; mY = 0; % responsive cells
    stp = 0.43; widths = [0.4 0.4 0.4 0.4 0.4 0.4]+0.1; gap = 0.09;
    adjust_axes(ff,[mY MY],stp,widths,gap,{'Estimated','Marginal Means'});
    % s = generate_shades(length(bins)-1);
    tcolors = [mData.colors(1:2);mData.dcolors(8)];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'B','NB','M'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_axes_top_text(ff.hf,ff.h_axes(1),'Complementary',[-0.04 0 0 0]);
    
    save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs_3_factors.pdf',ntrials),600);

    pos = get(gca,'Position');
%     text(1,ylims(2),sprintf('%s',event_type{ind}),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('unique_cells_across_conditions_%s_volunM.pdf',event_type{ind}),600);
    
    hf = get_figure(6,[10 7 1.5 1]);
    good_FRV = all_gFR;
    p_gFR_C = [];
    for rr = 1:size(good_FRV,2)
        for cc = 1:size(good_FRV,2)
            gFR_C = cell_list_op(good_FRV(:,rr),good_FRV(:,cc),'and',1);
            temp_gFR_C = 100*exec_fun_on_cell_mat(gFR_C,'sum')./exec_fun_on_cell_mat(gFR_C,'length');
            [m_p_gFR_C(rr,cc),sem_p_gFR_C(rr,cc)] = findMeanAndStandardError(temp_gFR_C);
            p_gFR_C(rr,cc) = mean(100*exec_fun_on_cell_mat(gFR_C,'sum')./exec_fun_on_cell_mat(gFR_C,'length'));
        end
    end
    
    gFR_C_all = cell_list_op(good_FRV(:,1),good_FRV(:,2),'and',1);    gFR_C_all = cell_list_op(gFR_C_all,good_FRV(:,3),'and',1);
    p_gFR_C_all = mean(100*exec_fun_on_cell_mat(gFR_C_all,'sum')./exec_fun_on_cell_mat(gFR_C_all,'length'));
    [mp_gFR_C_all,semp_gFR_C_all] = findMeanAndStandardError(100*exec_fun_on_cell_mat(gFR_C_all,'sum')./exec_fun_on_cell_mat(gFR_C_all,'length'));
%     AVenn = [p_gFR_C(1,1) p_gFR_C(2,2)]; IVenn = [p_gFR_C(1,2)];
%     [HVenn] = venn(AVenn,IVenn,'ErrMinMode','None');
    [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(good_FRV,0.5,0.05);
    AVenn = [p_gFR_C(1,1) p_gFR_C(2,2) p_gFR_C(3,3)]; IVenn = [p_gFR_C(1,2) p_gFR_C(1,3) p_gFR_C(2,3) p_gFR_C_all];
    [HVenn] = venn(AVenn,IVenn,'ErrMinMode','None');
    format_axes(gca);
    axis equal; axis off;
    changePosition(gca,[0.0 -0.0 -0.05 -0.05]); %put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -5 0]});
    pmchar=char(177);
%     text(-3.65,2.5,{'Brake',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(1,1),pmchar,sem_p_gFR_C(1,1))},'FontSize',6,'rotation',0);
%     text(4.5,-2.65,{'No-Brake',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(2,2),pmchar,sem_p_gFR_C(2,2))},'FontSize',6,'rotation',0);
%     text(0,4,{sprintf('Volun (%0.0f%c%0.0f%%)',m_p_gFR_C(3,3),pmchar,sem_p_gFR_C(3,3))},'FontSize',6,'rotation',0);
    set(HVenn(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); 
    set(HVenn(2),'FaceColor',tcolors{2},'FaceAlpha',0.75);
    set(HVenn(3),'FaceColor',tcolors{3},'FaceAlpha',0.75);
    ylims = ylim;
%     text(1,ylims(2)+0.51,sprintf('%s',event_type{ind}),'FontSize',6);
%     text(5,ylims(2)-0.51,sprintf('%.0f%%',IVenn(end)),'FontSize',6);
    btxt = (sprintf('%0.1f%c%0.1f%%',m_p_gFR_C(1,1),pmchar,sem_p_gFR_C(1,1)));
    nbtxt = (sprintf('%0.1f%c%0.1f%%',m_p_gFR_C(2,2),pmchar,sem_p_gFR_C(2,2)));
    inttxt = (sprintf('%0.1f%c%0.1f%%',m_p_gFR_C(1,2),pmchar,sem_p_gFR_C(1,2)));
    clc
    disp(btxt);
    disp(inttxt);
    disp(nbtxt);
    save_pdf(hf,mData.pdf_folder,sprintf('conjunctive_cells_venn_diagram_%s_volunM.pdf',event_type{ind}),600);

    %%
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(all_gFR,0.5,0.05);
    
    allrespcells = cell_list_op(all_gFR,[],'or',1);
    pallrespcells = 100*exec_fun_on_cell_mat(allrespcells,'sum')./exec_fun_on_cell_mat(allrespcells,'length');
    [mpar,sempar] = findMeanAndStandardError(pallrespcells);
    
    per_conj(:,1) = squeeze(all_CI_mat(1,2,:));
    per_conj(:,2) = squeeze(all_CI_mat(1,3,:));
    per_conj(:,3) = squeeze(all_CI_mat(2,3,:));
    
    [within,dvn,xlabels] = make_within_table({'Cond'},size(good_FR,2));
    dataT = make_between_table({per_conj},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = make_xdata([size(good_FR,2)],[1]);

    
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
    set(gcf,'color','w');    set(gcf,'Position',[10 3 1.6 1.25]);
    MY = 25; ysp = 1; mY = 0; % responsive cells
    stp = 0.43; widths = [0.4 0.4 0.4 0.4 0.4 0.4]+0.1; gap = 0.09;
    adjust_axes(ff,[mY MY],stp,widths,gap,{'Estimated','Marginal Means'});
    % s = generate_shades(length(bins)-1);
    tcolors = [mData.dcolors(1:2);mData.dcolors(9)];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'B-NB','B-M','NB-M'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_axes_top_text(ff.hf,ff.h_axes(1),'Conjunctive',[0.0 0 0 0]);

%     text(1,ylims(2),sprintf('%s',event_type{ind}),'FontSize',6);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('unique_cells_across_conditions_%s_volunM.pdf','L'),600);
    
    %%
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(all_gFR,0.5,0.05);
    mOI = mCI; semOI = semCI;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
%     txl = {'B-A-On-Exc','B-A-On-Inh','B-A-On-Unt','B-A-Off-Exc','B-A-Off-Inh','B-A-Off-Unt','NB-A-On-Exc','NB-A-On-Inh','NB-A-On-Unt','NB-A-Off-Exc','NB-A-Off-Inh','NB-A-Off-Unt'};%,'M-On','M-Off'};
    txl = {'Brake','No-Brake','Motion'};
%     txl = {'Air-On','Air-Off'};txl = repmat(txl,1,2);
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.25 1.25]);
%     hf = get_figure(5,[8 7 1.25 1.75]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
%     set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    for rr = 1:size(mCI,1)
        for cc = 1:size(mCI,1)
            if rr == cc
                text(rr-0.25,cc,sprintf('%.0f',mCI(rr,cc)),'FontSize',5);
            end
        end
    end
%     text(-0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); 
    set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[-0.01 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',ntrials),600);
    
    
    
    %% agglomerative hierarchical clustering

    mOI1 = mOI;
%     mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1,@naneucdist);
%     tree = linkage(mOI1,'average','euclidean');
    tree = linkage(Di,'average');
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 1.5 1]);
%     set(hf,'Position',[7 3 2 1]);
    set(H,'linewidth',0.5);
    txl = {'Brake','No-Brake','Motion'};
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel({'Eucledian','Distance'});changePosition(hx,[0 -0.3 0]);
    changePosition(gca,[0.15 0.0 -0.1 -0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);

    
    %%
    tcolors = mData.colors;
    hf = get_figure(6,[10 7 1.5 1]);
    good_FRV = all_gFR(:,[1 3]); 
    [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(good_FRV,0.5,0.05);
    p_gFR_C = [];
    for rr = 1:size(good_FRV,2)
        for cc = 1:size(good_FRV,2)
            gFR_C = cell_list_op(good_FRV(:,rr),good_FRV(:,cc),'and',1);
            temp_gFR_C = 100*exec_fun_on_cell_mat(gFR_C,'sum')./exec_fun_on_cell_mat(gFR_C,'length');
            [m_p_gFR_C(rr,cc),sem_p_gFR_C(rr,cc)] = findMeanAndStandardError(temp_gFR_C);
            p_gFR_C(rr,cc) = mean(100*exec_fun_on_cell_mat(gFR_C,'sum')./exec_fun_on_cell_mat(gFR_C,'length'));
        end
    end
    
    gFR_C_all = cell_list_op(good_FRV(:,1),good_FRV(:,2),'and',1);%    gFR_C_all = cell_list_op(gFR_C_all,good_FRV(:,3),'and',1);
    p_gFR_C_all = mean(100*exec_fun_on_cell_mat(gFR_C_all,'sum')./exec_fun_on_cell_mat(gFR_C_all,'length'));
    [mp_gFR_C_all,semp_gFR_C_all] = findMeanAndStandardError(100*exec_fun_on_cell_mat(gFR_C_all,'sum')./exec_fun_on_cell_mat(gFR_C_all,'length'));
    AVenn = [p_gFR_C(1,1) p_gFR_C(2,2)]; IVenn = [p_gFR_C(1,2)];
    [HVenn] = venn(AVenn,IVenn,'ErrMinMode','None');
%     AVenn = [p_gFR_C(1,1) p_gFR_C(2,2) p_gFR_C(3,3)]; IVenn = [p_gFR_C(1,2) p_gFR_C(1,3) p_gFR_C(2,3) p_gFR_C_all];
    [HVenn] = venn(AVenn,IVenn,'ErrMinMode','None');
    format_axes(gca);
    axis equal; axis off;
    changePosition(gca,[0.0 -0.0 -0.05 -0.05]); %put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -5 0]});
    pmchar=char(177);
%     text(-3.65,2.5,{'Brake',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(1,1),pmchar,sem_p_gFR_C(1,1))},'FontSize',6,'rotation',0);
%     text(4.5,-2.65,{'No-Brake',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(2,2),pmchar,sem_p_gFR_C(2,2))},'FontSize',6,'rotation',0);
    set(HVenn(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); set(HVenn(2),'FaceColor',tcolors{2},'FaceAlpha',0.75);
%     text(-6,6,{'Volunteer',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(3,3),pmchar,sem_p_gFR_C(3,3))},'FontSize',6,'rotation',0);
    ylims = ylim;
%     text(1,ylims(2)+0.51,sprintf('%s',event_type{ind}),'FontSize',6);
%     text(5,ylims(2)-0.51,sprintf('%.0f%%',IVenn),'FontSize',6);
    btxt = (sprintf('%0.1f%c%0.1f%%',m_p_gFR_C(1,1),pmchar,sem_p_gFR_C(1,1)));
    nbtxt = (sprintf('%0.1f%c%0.1f%%',m_p_gFR_C(2,2),pmchar,sem_p_gFR_C(2,2)));
    inttxt = (sprintf('%0.1f%c%0.1f%%',m_p_gFR_C(1,2),pmchar,sem_p_gFR_C(1,2)));
%     ht = title(sprintf('%s   %s   %s',btxt,inttxt,nbtxt));set(ht,'FontWeight','normal','FontSize',5);
    clc
    disp(btxt);
    disp(inttxt);
    disp(nbtxt);
    save_pdf(hf,mData.pdf_folder,sprintf('conjunctive_cells_venn_diagram_%s.pdf',''),600);
    
    
    intcells = squeeze(all_CI_mat(1,2,:));    bcells = squeeze(all_CI_mat(1,1,:)); nbcells = squeeze(all_CI_mat(2,2,:));
    bcellsu = squeeze(uni(1,2,:)); nbcellsu = squeeze(uni(2,1,:));
    [within,dvn,xlabels] = make_within_table({'Type'},[3]);
    dataT = make_between_table({[bcellsu intcells nbcellsu]},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Type','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = [mData.colors(1) mean([mData.colors{1};mData.colors{2}]) mData.colors(2)];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',7,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
%     set(hbs(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); set(hbs(3),'FaceColor',tcolors{2},'FaceAlpha',0.75);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'No-Brake','Conjunc','Motion'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.1 0.01 -0.3 0]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all3.pdf',ntrials),600);
    %%
    break;
end


%% for heat map and clustering
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    event_type = {'1-L','2-AOn','2-AOff','2-Arb','3-AOn','3-A','3-AOff','3-Arb','4-AOn','4-AL','4-AOff','4-Arb','5-AOn','5-A','5-AOff','5-Arb','6-L','7-AOn','7-AOff','7-Arb','M-On','M-Off'};
    sic = {Lb;Ab_On;Ab_Off;Ab_Offc;Ar_On;Ar_L;Ar_Off;Ar_Offc;ArL_On;ArL_L;ArL_Off;ArL_Offc;Ars_On;Ars_L;Ars_Off;Ars_Offc;Lbs;Abs_On;Abs_Off;Abs_Offc;M_On;M_Off};
%     event_type = rasterNamesTxt([Lb Ab_On Ab_Off Ab_Offc Ar_On Ar_L Ar_Off Ar_Offc ArL_On ArL_L ArL_Off ArL_Offc Ars_On Ars_L Ars_Off Ars_Offc Lbs Abs_On Abs_Off Abs_Offc M_On M_Off]);
%     sic = {Lb;Ab_On;Ab_Off;Ab_Offc;Ar_On;Ar_L;Ar_Off;Ar_Offc;ArL_On;ArL_L;ArL_Off;ArL_Offc;Ars_On;Ars_L;Ars_Off;Ars_Offc;Lbs;Abs_On;Abs_Off;Abs_Offc;M_On;M_Off};
%     
%     event_type = {'B-AOn','B-AOff','B-Arb','NB-AOn','NB-AOff','NB-Arb','B-L','NB-AL','NB-A','M-On','M-Off'};
%     sic = {[Ab_On Abs_On];[Ab_Off Abs_Off];[Ab_Offc Abs_Offc];[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off];[Ar_Offc ArL_Offc Ars_Offc];[Lb Lbs];[ArL_L];[Ar_L Ars_L];[M_On];[M_Off]}; % for heat map after light sitmulus
%     
%     event_type = {'2-B-AOn','7-B-AOn','2-B-AOff','7-B-AOff','3-NB-AOn','4-NB-AOn','5-NB-AOn','3-NB-AOff',...
%         '4-NB-AOff','5-NB-AOff','1-B-L','6-B-L','M-On','M-Off'};
%     sic = {[Ab_On];[Abs_On];[Ab_Off];[Abs_Off];[Ar_On];[ArL_On];[Ars_On];[Ar_Off];[ArL_Off];[Ars_Off];[Lb];[Lbs];[M_On];[M_Off]}; 
%     
    
    clear all_gFR all_exc all_inh all_gV 
    prop_names = {'resp','N_Resp_Trials','zMI','zMINaN','HaFD','HiFD','cells_pooled'};
    cell_sel = {'good_FR','vals','exc','inh'};
    varName = {'all_gFR','all_gV','all_exc','all_inh'};
    for cii = 1:length(cell_sel);
        cmdTxt = sprintf('clear %s',varName{cii});eval(cmdTxt);
    end
    pni = 7;
    all_exc_inh = [];
    for ii = 1:length(sic)
        sit = sic{ii};
        tRs = o.Rs(:,sit);
        props1 = get_props_Rs(tRs,ntrials);
        for cii = 1:length(cell_sel)
            cmdTxt = sprintf('gFR = props1.%s;',cell_sel{cii});eval(cmdTxt);
            if pni == 1
                rf = find_percent(gFR);
                cmdTxt = sprintf('%s(:,ii) = mean(rf,2);',varName{cii}); eval(cmdTxt)
            else
                if pni == 7
                    cmdTxt = sprintf('%s(:,ii) = cell_list_op(gFR,[],''or'',1);',varName{cii});eval(cmdTxt);
                else
                    if pni == 3 || pni == 4
                        cmdtxt = sprintf('rf = props1.%s;',prop_names{3}); eval(cmdtxt);
                        if pni == 3
                            cmdTxt = sprintf('%s(:,ii) = mean(exec_fun_on_cell_mat(rf,''nanmean'',gFR),2);',varName{cii});eval(cmdTxt);
                        else
                            for rrr = 1:size(rf,1)
                                for ccc = 1:size(rf,2)
                                    temp = isnan(rf{rrr,ccc});
                                    rf1(rrr,ccc) = 100*sum(temp(gFR{rrr,ccc}))/size(rf{rrr,ccc},1);
                                end
                            end
                            cmdTxt = sprintf('%s(:,ii) = mean(rf1,2);',varName{cii}); eval(cmdTxt)
                        end
                    else
                        cmdtxt = sprintf('rf = props1.%s;',prop_names{pni}); eval(cmdtxt);
                        cmdTxt = sprintf('%s(:,ii) = mean(exec_fun_on_cell_mat(rf,''mean'',gFR),2);',varName{cii});eval(cmdTxt);
                    end
                end
            end
        end
        all_exc_inh = [all_exc_inh all_exc(:,ii) all_inh(:,ii)];
    end

    %%
    break;
end


%% three circle Venn diagram, percentages of complementary and conjunctive cells compared among Lb/Lbs, Ar_L/Ars_L, and ArL_L
while 1
    all_resp = [];
    ntrials = 50;
    sic = {[Ab_On Abs_On Ab_Off Abs_Off Ab_Offc Abs_Offc Lb Lbs];[Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off Ar_Offc ArL_Offc Ars_Offc ArL_L Ar_L Ars_L];[M_On M_Off]}; % for heat map after light sitmulus
    clear all_gFR
    for ii = 1:length(sic)
        sit = sic{ii};
        tRs = o.Rs(:,sit);
        props1 = get_props_Rs(tRs,ntrials);
        gFR = props1.vals;
        gFR = cell_list_op(gFR,[],'or',1);
        all_gFR(:,ii) = gFR;
    end
    all_resp = [all_resp all_gFR];
    good_FR = all_gFR;
    per_unique =[];
    condMat = [];
    for rr = 1:size(good_FR,2)
        for cc = 1:size(good_FR,2)
            if rr == cc
                condMat(cc) = cc;
            else
                condMat(cc) = -cc;
            end
        end
        temp_unique = get_cell_list(good_FR,condMat,1);
        unique_gFR(:,rr) = get_cell_list(good_FR,condMat,0);
        per_unique(:,rr) = temp_unique(:,1);
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},size(good_FR,2));
    dataT = make_between_table({per_unique},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    
    any_gFR = cell_list_op(all_gFR,[],'or',1);
    pangfr = 100*exec_fun_on_cell_mat(any_gFR,'sum')./exec_fun_on_cell_mat(any_gFR,'length');
    [mpangfr,sempangfr] = findMeanAndStandardError(pangfr);
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = make_xdata([size(good_FR,2)],[1]);
%     h(h==1) = 0;
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
    set(gcf,'color','w');    set(gcf,'Position',[10 3 1.6 1.25]);
    MY = 39; ysp = 1.5; mY = 0; % responsive cells
    stp = 0.43; widths = [0.4 0.4 0.4 0.4 0.4 0.4]+0.1; gap = 0.09;
    adjust_axes(ff,[mY MY],stp,widths,gap,{'Estimated','Marginal Means'});
    % s = generate_shades(length(bins)-1);
    tcolors = [mData.colors(1:2);mData.dcolors(8)];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'B','M','NB'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    set_axes_top_text(ff.hf,ff.h_axes(1),'Complementary',[-0.04 0 0 0]);
    
    save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs_3_factors.pdf',ntrials),600);

    pos = get(gca,'Position');
%     text(1,ylims(2),sprintf('%s',event_type{ind}),'FontSize',6);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graph.pdf'),600);
    
    hf = get_figure(6,[10 7 1.5 1]);
    good_FRV = all_gFR;
    p_gFR_C = [];
    for rr = 1:size(good_FRV,2)
        for cc = 1:size(good_FRV,2)
            gFR_C = cell_list_op(good_FRV(:,rr),good_FRV(:,cc),'and',1);
            temp_gFR_C = 100*exec_fun_on_cell_mat(gFR_C,'sum')./exec_fun_on_cell_mat(gFR_C,'length');
            [m_p_gFR_C(rr,cc),sem_p_gFR_C(rr,cc)] = findMeanAndStandardError(temp_gFR_C);
            p_gFR_C(rr,cc) = mean(100*exec_fun_on_cell_mat(gFR_C,'sum')./exec_fun_on_cell_mat(gFR_C,'length'));
        end
    end
    
    gFR_C_all = cell_list_op(good_FRV(:,1),good_FRV(:,2),'and',1);    gFR_C_all = cell_list_op(gFR_C_all,good_FRV(:,3),'and',1);
    p_gFR_C_all = mean(100*exec_fun_on_cell_mat(gFR_C_all,'sum')./exec_fun_on_cell_mat(gFR_C_all,'length'));
    [mp_gFR_C_all,semp_gFR_C_all] = findMeanAndStandardError(100*exec_fun_on_cell_mat(gFR_C_all,'sum')./exec_fun_on_cell_mat(gFR_C_all,'length'));
%     AVenn = [p_gFR_C(1,1) p_gFR_C(2,2)]; IVenn = [p_gFR_C(1,2)];
%     [HVenn] = venn(AVenn,IVenn,'ErrMinMode','None');
    [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(good_FRV,0.5,0.05);
    AVenn = [p_gFR_C(1,1) p_gFR_C(2,2) p_gFR_C(3,3)]; IVenn = [p_gFR_C(1,2) p_gFR_C(1,3) p_gFR_C(2,3) p_gFR_C_all];
    [HVenn] = venn(AVenn,IVenn,'ErrMinMode','None');
    format_axes(gca);
    axis equal; axis off;
    changePosition(gca,[0.0 -0.0 -0.05 -0.05]); %put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -5 0]});
    pmchar=char(177);
%     text(-3.65,2.5,{'Brake',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(1,1),pmchar,sem_p_gFR_C(1,1))},'FontSize',6,'rotation',0);
%     text(4.5,-2.65,{'No-Brake',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(2,2),pmchar,sem_p_gFR_C(2,2))},'FontSize',6,'rotation',0);
%     text(0,4,{sprintf('Volun (%0.0f%c%0.0f%%)',m_p_gFR_C(3,3),pmchar,sem_p_gFR_C(3,3))},'FontSize',6,'rotation',0);
    set(HVenn(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); 
    set(HVenn(2),'FaceColor',tcolors{2},'FaceAlpha',0.75);
    set(HVenn(3),'FaceColor',tcolors{3},'FaceAlpha',0.75);
    ylims = ylim;
%     text(1,ylims(2)+0.51,sprintf('%s',event_type{ind}),'FontSize',6);
%     text(5,ylims(2)-0.51,sprintf('%.0f%%',IVenn(end)),'FontSize',6);
    btxt = (sprintf('%0.1f%c%0.1f%%',m_p_gFR_C(1,1),pmchar,sem_p_gFR_C(1,1)));
    nbtxt = (sprintf('%0.1f%c%0.1f%%',m_p_gFR_C(2,2),pmchar,sem_p_gFR_C(2,2)));
    inttxt = (sprintf('%0.1f%c%0.1f%%',m_p_gFR_C(1,2),pmchar,sem_p_gFR_C(1,2)));
    clc
    disp(btxt);
    disp(inttxt);
    disp(nbtxt);
    save_pdf(hf,mData.pdf_folder,sprintf('Venn_Diagram.pdf'),600);

    
    
    %%
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(all_gFR,0.5,0.05);
    
    per_comp1(:,1) = squeeze(uni(1,2,:));
    per_comp1(:,2) = squeeze(uni(1,3,:));
    per_comp1(:,3) = squeeze(uni(2,3,:));
    
    per_comp2(:,1) = squeeze(uni(2,1,:));
    per_comp2(:,2) = squeeze(uni(3,1,:));
    per_comp2(:,3) = squeeze(uni(3,2,:));
    
    [within,dvn,xlabels] = make_within_table({'PT','CondP'},[3,3]); inds = [1 2 3];
    dataTP = make_between_table({[per_comp1(:,[inds]) per_conj(:,[inds]) per_comp2(:,[inds])]},dvn);
    raP = RMA(dataTP,within,{0.05,{'hsd','bonferroni'}});
    raP.ranova
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,raP,{'CondP_by_PT','bonferroni'},[1 1 1]); [combs,p,h] = reduce_combs(combs,p,h,{[1:3],[4:6],[7:9]});
    xdata = make_xdata([3 3 3],[1 1.5]);
    
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.36],'widthHeightAdjustment',[10 -410]);
    set(gcf,'color','w');    set(gcf,'Position',[10 3 1.93 1.25]);
    MY = 70; ysp = 2; mY = 0; % responsive cells
    stp = 0.37; widths = [1.35 0.4 0.4 0.4 0.4 0.4]+0.1; gap = 0.09;
    adjust_axes(ff,[mY MY],stp,widths,gap,{'Estimated','Marginal Means'});
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.conj_comp_colors,3,1);
    axes(ff.h_axes(1,1));
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[mY MY]); format_axes(gca);xticks = xdata; xticklabels = {'Comp1','Conj','Comp2'};set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
%     put_sig_lines(gcf,gca,hbs(1:3),mVar(1:3)
%     set_axes_top_text(ff.hf,ff.h_axes(1),'Conjunctive',[0.0 0 0 0]);
    set_bar_graph_sub_xtick_text(gcf,gca,hbs,3,{'B-NB','B-M','NB-M'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('unique_cells_across_conditions_%s_volunM.pdf','L'),600);
    
   
    %%
    break;
end




% 
% 
% 
% while 1
%     titles = {'M-On','M-Off'};
%     si = [M_On M_Off];
%     Rs = o.Rs(:,si);mR = o.mR(:,si);
%     ntrials = 50;
%     props1 = get_props_Rs(Rs,ntrials);
%     %%
%     break;
% end
% %% See onsets and offsets
% while 1
%     an = 3; cn = 1;
%     motionOnsets = Rs{an,cn}.onsets;
%     motionOffsets = Rs{an,cn}.offsets;
%     hf = figure(1000);clf;
%     set(gcf,'color','w','units','inches'); set(gcf,'Position',[7 4 1.5 1]);
%     display_with_air_puff(ei{an}.b,motionOnsets,motionOffsets);
%     xlim([270 330]/60);ylim([0 2]);
%     plot([4.6 4.8],[1.45 1.45],'b');
%     plot([4.6 4.8],[1.9 1.9],'m');
%     changePosition(gca,[-0.07 0.15 0 -0.15]); box off;
%     ax = gca; ax.YAxis.Visible = 'off';
%     xlabel('Time (min)');
%     set(gca,'FontSize',6,'FontWeight','Normal','TickDir','out');
%     save_pdf(hf,mData.pdf_folder,sprintf('motionOnset_trials'),600);
%     break;
% end
% %% Show sample rasters
% an = 4; cn = 1;
% figure(2000);clf;imagesc(Rs{an,cn}.speed);colorbar;
% % plotRasters_simplest(Rs{an,cn})
% props1 = get_props_Rs(Rs,ntrials);
% plotRasters_simplest(Rs{an,cn},find(props1.inh{an,cn}))
% %%
% an = 4; cn = 1;
% ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
%         'spaceRowsCols',[0.15 0.04],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
%         [-60 -475]);
%     set(gcf,'color','w'); set(gcf,'Position',[10 4 3.4 1]);
% ff = sample_rasters(Rs{an,cn},[248 57 264 270 230],ff);
% for ii = 1:4
%         set(ff.h_axes(1,ii),'xtick',[1 13.5 27],'xticklabels',{'-1.5','0','1.5'});
%     end
% save_pdf(ff.hf,mData.pdf_folder,sprintf('motion_rastersD'),600);
% %% population vector and correlation
% while 1
%     selected_property = 'tuned';
%     an = 1;
%     resp = props1.vals;
% %     eval(cmdTxt);
%     ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 2],...
%         'spaceRowsCols',[0.01 0.05],'rightUpShifts',[0.16 0.11],'widthHeightAdjustment',...
%         [-130 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 1.7 1.5]);
%     [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
%     ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
%     for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
%     colormap_ig
%     for ii = 1:2
%         set(ff.h_axes(2,ii),'xtick',[1 13.5 27],'xticklabels',{'-1.5','0','1.5'});
%     end
%     save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_motion_%s.pdf',selected_property),600);
% 
%     % average correlation of all animals
%     ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 2],...
%         'spaceRowsCols',[0.01 0.05],'rightUpShifts',[0.16 0.25],'widthHeightAdjustment',...
%         [-130 -320]);    set(gcf,'color','w');    set(gcf,'Position',[10 8 1.7 0.8]);
%     ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
%     colormap_ig
%     for ii = 1:2
%         set(ff.h_axes(1,ii),'xtick',[1 13.5 27],'xticklabels',{'-1.5','0','1.5'},'ytick',[1 13.5 27],'yticklabels',{'-1.5','0','1.5'});
%     end
%     save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_motion_%s.pdf',selected_property),600);
%     %%
%     break;
% end
% 
% %% number of trials
% while 1
%     for rr = 1:size(Rs,1)
%         for cc = 1:2
%             tRs = Rs{rr,cc};
%             num_trials(rr,cc) = size(tRs.sp_rasters1,1);
%         end
%     end
%     mean(num_trials)
%     std(num_trials)/sqrt(5)
%     %%
%     break;
% end
% 
% %% average % of trials in which the cell responded
% while 1
%     mean_N_trials_resp = exec_fun_on_cell_mat(props1.N_Resp_Trials,'mean');
%     [within,dvn,xlabels] = make_within_table({'Cond'},[length(si)]);
%     dataT = make_between_table({mean_N_trials_resp},dvn);
%     ra = RMA(dataT,within);
%     ra.ranova
%     %%
%     [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
%     xdata = make_xdata([2],[1]);
%     hf = get_figure(5,[8 7 1.25 1]);
%     % s = generate_shades(length(bins)-1);
%     tcolors = mData.colors;
%     [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
%         'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
%         'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY + 5;
%     ylims = ylim;
%     format_axes(gca);
%     set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
%     xticks = xdata; xticklabels = {'MOn','MOff'};
%     set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
%     changePosition(gca,[0.065 0.0 -0.4 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Trials (%)',[0 0 0]});
%     save_pdf(hf,mData.pdf_folder,sprintf('motion_trials.pdf'),600);
%     %%
%     break;
% end
% %% Percentage of Responsive Cells
% while 1
%     resp = props1.vals;
%     perc_r = 100*exec_fun_on_cell_mat(resp,'sum')./exec_fun_on_cell_mat(resp,'length');
%     within = make_within_table({'Cond'},2);
%     dataT = make_between_table({perc_r},{'M_On','M_Off'});
%     ra = RMA(dataT,within);
%     ra.ranova
%     
%     [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
%      hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w'); hold on;
%     tcolors = mData.dcolors(17:20);
%     [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
%         'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
%         'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.05);
%     set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
% %     hatch(hbs(2),135,'k','-',2,0.1);
%     xticks = [xdata(1:end)]; xticklabels = {'Onset','Offset'};
%     set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30)
%     changePosition(gca,[0.2 0.03 -0.5 -0.11]);
%     put_axes_labels(gca,{[],[0 0 0]},{{'Motion Responsive','Cells (%)'},[0 0 0]});
%     save_pdf(hf,mData.pdf_folder,sprintf('percentage_motion_responsive'),600);
%     break;
% end
% %% Percentage of excitatory inhibitory responsive cells
% while 1
%     resp1 = props1.exc;     resp2 = props1.inh;
%     presp1 = 100 * exec_fun_on_cell_mat(resp1,'sum')./exec_fun_on_cell_mat(resp1,'length');
%     presp2 = 100 * exec_fun_on_cell_mat(resp2,'sum')./exec_fun_on_cell_mat(resp2,'length');
%     
%     [within,dvn,xlabels] = make_within_table({'MT','CT'},[2,2]);
%     dataT = make_between_table({presp1(:,1),presp2(:,1),presp1(:,2),presp2(:,2)},dvn);
%     ra = RMA(dataT,within);
%     ra.ranova
%     %%
%     [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'MT_by_CT','hsd'},[1 1 1]);
%     h([2 5]) = 0;
%     xdata = make_xdata([2,2],[1 2]);
%     hf = get_figure(5,[8 7 1.25 1]);
%     % s = generate_shades(length(bins)-1);
%     tcolors = mData.colors(3:end);
%     [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
%         'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
%         'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.001);
%     maxY = maxY + 3;
%     ylims = ylim;
%     format_axes(gca);
%     set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 22]); format_axes(gca);
%     xticks = xdata; xticklabels = {'Exc','Inh'};
%     set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
%     changePosition(gca,[0.1 0.0 -0.2 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
%     
%     save_pdf(hf,mData.pdf_folder,'motion_responsive_exc_inh',600);
%     
%     %%
%     [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
% 
%     xdata = make_xdata([2],[1 2]);
%     hf = get_figure(5,[8 7 1.25 1]);
%     % s = generate_shades(length(bins)-1);
%     tcolors = mData.colors(7:end);
%     [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
%         'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
%         'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY + 3;
%     ylims = ylim;
%     format_axes(gca);
%     set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 22]); format_axes(gca);
%     xticks = xdata; xticklabels = {'Exc','Inh'};
%     set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[]); xtickangle(45)
%     changePosition(gca,[-0.1 0.0 -0.5 -0.05]); %put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
%     
%     save_pdf(hf,mData.pdf_folder,'motion_responsive_exc_inh',600);
%     %%
%     break;
% end
% 
% 
% %% Percentage response fidelity excitatory inhibitory responsive cells
% while 1
%     gFR = props1.exc; rf = props1.N_Resp_Trials;
%     all_exc = exec_fun_on_cell_mat(rf,'mean',gFR);
%     gFR = props1.inh;rf = props1.N_Resp_Trials;
%     all_inh = exec_fun_on_cell_mat(rf,'mean',gFR);
%     
%     [within,dvn,xlabels] = make_within_table({'MT','CT'},[2,2]);
%     dataT = make_between_table({[all_exc(:,1),all_inh(:,1),all_exc(:,2),all_inh(:,2)]},dvn);
%     ra = RMA(dataT,within);
%     ra.ranova
%     
%     [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'MT_by_CT','hsd'},[1 1 1]);
%     h(h==1) = 0;
%     xdata = make_xdata([2,2],[1 1.5]);
%     hf = get_figure(5,[8 7 1.25 1]);
%     % s = generate_shades(length(bins)-1);
%     tcolors = mData.colors(3:end);
%     [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
%         'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
%         'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY + 5;
%     ylims = ylim;
%     format_axes(gca);
%     set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 50]); format_axes(gca);
%     xticks = xdata; xticklabels = {'Exc','Inh'};
%     set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 20 40]); xtickangle(45)
%     changePosition(gca,[0.1 0.0 -0.25 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Trials (%)'},[0 0 0]});
%     save_pdf(hf,mData.pdf_folder,'motion_responsive_exc_inh',600);
%     
%     %%
%     [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'MT','hsd'},[1 1 1]);
%     xdata = make_xdata([2],[1 1.5]);
%     hf = get_figure(5,[8 7 1.25 1]);
%     % s = generate_shades(length(bins)-1);
%     tcolors = mData.colors(7:end);
%     [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
%         'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
%         'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY + 5;
%     ylims = ylim;
%     format_axes(gca);
%     set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 50]); format_axes(gca);
%     xticks = xdata; xticklabels = {'MOn','MOff'};
%     set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 20 40]); xtickangle(45)
%     changePosition(gca,[0.1 0.0 -0.5 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Trials (%)'},[0 0 0]});
%     save_pdf(hf,mData.pdf_folder,'motion_responsive_exc_inh',600);
%     %%
%     break;
% end
% 
% %% Overlap Indices ImageSC 
% while 1
%     ntrials = 50;
%     si = [Lb_T ArL_L_T Lbs_T Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
%     props1A = get_props_Rs(o.Rs,ntrials);
%     respO = [props1A.good_FR(:,si)];% resp_speed];
%     resp = [props1.vals resp_speedAcc respO];% resp_speed];
%     [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
%     sz = size(mOI,1);
% %     mOI = OI_mat(:,:,4);
%     oM = ones(size(mOI));
%     mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
%     maxI = max([mOI(:);semOI(:)]);    
%     minI = min([mOI(:);semOI(:)]);
%     
%     mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
%     txl = [{'M-On','M-Off','Speed','Accel'} rasterNamesTxt(si)]; 
% %     mOI = mOI .* mask;
%     imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
%     imAlpha(mask1 == 1) = 0;
% %     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
%     hf = get_figure(5,[8 7 2 2]);
%     %
% %     axes(ff.h_axes(1));
%     im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
%     set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
%     format_axes(gca);
%     set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
%     set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
%     text(0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
%     box on
%     changePosition(gca,[0.03 0 -0.04 0]);
%     hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
%     save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean_motion.pdf',ntrials),600);
%     %%
%     break;
% end
% 
% %% agglomerative hierarchical clustering zMIT>zMID
% while 1
%     mOI1 = mOI;
%     mOI1(isnan(mOI1)) = 1;
%     Di = pdist(mOI1);
%     tree = linkage(Di);
%     figure(hf);clf
%     [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold','default');
%     hf = gcf;
%     set(hf,'Position',[7 3 2.5 1]);
%     set(H,'linewidth',1);
%     set(gca,'xticklabels',txl(TC));xtickangle(45);
%     format_axes(gca);
%     hx = ylabel({'Eucledian','Distance'});%changePosition(hx,[-0.051 0 0]);
%     changePosition(gca,[0 0.0 0.07 0.05]);
%     save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster_motion.pdf'),600);
%     %%
%     break;
% end
% 
% 
% 
% %% Overlap Indices ImageSC for presentation
% while 1
%     ntrials = 50;
%     si = [Lb_T ArL_L_T Lbs_T Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D];
%     props1A = get_props_Rs(o.Rs,ntrials);
%     respO = [props1A.good_FR(:,si)];% resp_speed];
%     resp = [props1.vals respO];% resp_speed];
%     [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
%     sz = size(mOI,1);
% %     mOI = OI_mat(:,:,4);
%     oM = ones(size(mOI));
%     mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
%     maxI = max([mOI(:);semOI(:)]);    
%     minI = min([mOI(:);semOI(:)]);
%     
%     mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
%     txl = [{'M-On','M-Off'} rasterNamesTxt(si)]; 
% %     mOI = mOI .* mask;
%     imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
%     imAlpha(mask1 == 1) = 0;
% %     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
%     hf = get_figure(5,[8 7 1.5 1.5]);
%     %
% %     axes(ff.h_axes(1));
%     im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
%     set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
%     format_axes(gca);
%     set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
%     set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
%     text(-1,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
%     box on
%     changePosition(gca,[0.03 0 -0.04 0]);
%     hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
%     save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean_motion.pdf',ntrials),600);
%     %%
%     break;
% end
% 
% %% agglomerative hierarchical clustering 
% while 1
%     mOI1 = mOI;
%     mOI1(isnan(mOI1)) = 1;
%     Di = pdist(mOI1);
%     tree = linkage(Di);
%     figure(hf);clf
%     [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold','default');
%     hf = gcf;
%     set(hf,'Position',[7 3 2.5 1]);
%     set(H,'linewidth',1);
%     set(gca,'xticklabels',txl(TC));xtickangle(45);
%     format_axes(gca);
%     hx = ylabel({'Eucledian','Distance'});%changePosition(hx,[-0.051 0 0]);
%     changePosition(gca,[0 0.0 0.07 0.05]);
%     save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster_motion.pdf'),600);
%     %%
%     break;
% end
% 
