%% compare percent responsive cells pooled light ON brake vs no-brake exc and inh ArL_L
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    allsic = {{[Lb Lbs];[ArL_L]};
    };
    
%     ind = 3;
%     sic = allsic{1:6};
    all_resp = []; all_resp_exc = []; all_resp_inh = []; all_resp_unt = [];
    for jj = 1:length(allsic)
        sic = allsic{jj};
        clear all_gFR all_exc all_inh all_unt
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
            gFR = props1.good_FR_and_untuned;
            gFR = cell_list_op(gFR,[],'or',1);
            all_unt(:,ii) = gFR;
        end
        all_resp = [all_resp all_gFR]; all_resp_exc = [all_resp_exc all_exc]; all_resp_inh = [all_resp_inh all_inh];
        all_resp_unt = [all_resp_unt all_unt];
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
    break;
end

%% compare percent response fidelity pooled light on brake vs no-brake ArL_L
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    allsic = {{[Lb Lbs];[ArL_L]};
    };
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
            gFR = props1.inh; rf = props1.N_Resp_Trials;
            all_inh(:,ii) = mean(exec_fun_on_cell_mat(rf,'mean',gFR),2);
        end
        all_resp = [all_resp all_gFR]; all_resp_exc = [all_resp_exc all_exc]; all_resp_inh = [all_resp_inh all_inh];
    end
    good_FR = [all_resp_exc all_resp_inh];

    [within,dvn,xlabels] = make_within_table({'Type','Cond'},[2,2]);
    dataT = make_between_table({good_FR},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    break;
end


%% compare percent responsive cells pooled light ON brake vs no-brake exc and inh Ar_L
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    allsic = {{[Lb Lbs];[Ars_L]};
    };
    
%     ind = 3;
%     sic = allsic{1:6};
    all_resp = []; all_resp_exc = []; all_resp_inh = []; all_resp_unt = [];
    for jj = 1:length(allsic)
        sic = allsic{jj};
        clear all_gFR all_exc all_inh all_unt
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
            gFR = props1.good_FR_and_untuned;
            gFR = cell_list_op(gFR,[],'or',1);
            all_unt(:,ii) = gFR;
        end
        all_resp = [all_resp all_gFR]; all_resp_exc = [all_resp_exc all_exc]; all_resp_inh = [all_resp_inh all_inh];
        all_resp_unt = [all_resp_unt all_unt];
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
    break;
end

%% compare percent response fidelity pooled light on brake vs no-brake Ar_L
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    allsic = {{[Lb Lbs];[ArL_L]};
    };
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
            gFR = props1.inh; rf = props1.N_Resp_Trials;
            all_inh(:,ii) = mean(exec_fun_on_cell_mat(rf,'mean',gFR),2);
        end
        all_resp = [all_resp all_gFR]; all_resp_exc = [all_resp_exc all_exc]; all_resp_inh = [all_resp_inh all_inh];
    end
    good_FR = [all_resp_exc all_resp_inh];

    [within,dvn,xlabels] = make_within_table({'Type','Cond'},[2,2]);
    dataT = make_between_table({good_FR},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    break;
end

%% for heatmap of all air onset and offset and arb and light (conditions and cell types) GRAND one
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    event_type = {'B-AOn','B-AOff','B-Arb','B-L','NB-AOn','NB-AOff','NB-Arb','NB-L','NB-ArbL'};
    sic = {[Ab_On Abs_On];[Ab_Off Abs_Off];[Ab_Offc Abs_Offc];[Lb Lbs];[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off];[Ar_Offc ArL_Offc Ars_Offc];[ArL_L];[Ar_L Ars_L];};
    
    event_type = {'2-B-AOn','7-B-AOn','2-B-AOff','7-B-AOff','2-B-Arb','7-B-Arb','1-B-L','6-B-L','3-NB-AOn','4-NB-AOn','5-NB-AOn',...
        '3-NB-AOff','4-NB-AOff','5-NB-AOff','3-NB-Arb','4-NB-Arb','5-NB-Arb','3-NB-AL','4-NB-L','5-NB-AL'};
    sic = {[Ab_On];[Abs_On];[Ab_Off];[Abs_Off];[Ab_Offc];[Abs_Offc];[Lb];[Lbs];[Ar_On];[ArL_On];[Ars_On];[Ar_Off];[ArL_Off];[Ars_Off];[Ar_Offc];[ArL_Offc];[Ars_Offc];[ArL_L];[Ar_L];[Ars_L]};
    
    event_type = {'2-B-AOn','7-B-AOn','2-B-AOff','7-B-AOff','2-B-Arb','7-B-Arb','1-B-L','6-B-L','3-NB-AOn','4-NB-AOn','5-NB-AOn',...
        '3-NB-OArb','4-NB-OArb','5-NB-OArb','3-NB-AOff','4-NB-AOff','5-NB-AOff','3-NB-Arb','4-NB-Arb','5-NB-Arb','3-NB-AL','4-NB-L','5-NB-AL'};
    sic = {[Ab_On];[Abs_On];[Ab_Off];[Abs_Off];[Ab_Offc];[Abs_Offc];[Lb];[Lbs];[Ar_On];[ArL_On];[Ars_On];[Ar_Onc];[ArL_Onc];[Ars_Onc];[Ar_Off];[ArL_Off];[Ars_Off];[Ar_Offc];[ArL_Offc];[Ars_Offc];[ArL_L];[Ar_L];[Ars_L]};
    
    event_type = {'2-B-AOn','7-B-AOn','2-B-AOff','7-B-AOff','2-B-Arb','7-B-Arb','1-B-L','6-B-L','3-NB-AOn','4-NB-AOn','5-NB-AOn',...
        '3-NB-AOff','4-NB-AOff','5-NB-AOff','3-NB-Arb','4-NB-Arb','5-NB-Arb','4-NB-L'};
    sic = {[Ab_On];[Abs_On];[Ab_Off];[Abs_Off];[Ab_Offc];[Abs_Offc];[Lb];[Lbs];[Ar_On];[ArL_On];[Ars_On];[Ar_Off];[ArL_Off];[Ars_Off];[Ar_Offc];[ArL_Offc];[Ars_Offc];[ArL_L]};
    
    
    clear all_gV all_gFR all_exc all_inh;
    all_exc_inh = []; xtl = []; 
    
    for ii = 1:length(sic)
        sit = sic{ii};
        tRs = o.Rs(:,sit);
        props1 = get_props_Rs(tRs,ntrials);
        gFR = props1.vals;
        gFR = cell_list_op(gFR,[],'or',1);
        all_gV(:,ii) = gFR;
        gFR = props1.good_FR_and_untuned;
        gFR = cell_list_op(gFR,[],'or',1);
        all_gFR(:,ii) = gFR;
        gFR = props1.exc; %gFR = props1.good_FR_and_exc;
        gFR = cell_list_op(gFR,[],'or',1);
        all_exc(:,ii) = gFR;
        gFR = props1.inh; %gFR = props1.good_FR_and_inh;
        gFR = cell_list_op(gFR,[],'or',1);
        all_inh(:,ii) = gFR;
        all_exc_inh = [all_exc_inh all_exc(:,ii) all_inh(:,ii)];
        xtl = [xtl {sprintf('%s-Exc',event_type{ii})} {sprintf('%s-Inh',event_type{ii})}];
    end
    %%
    break;
end

%% compare percent responsive cells pooled light ON brake vs no-brake exc and inh and Control Light (from Condition 3 and 5 where there is no light stimulus)
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    allsic = {{[Lb Lbs];[ArL_L];[Ar_L];[Ars_L]};
    };
    allsic = {{[ArL_L];[Ar_L];[ Ars_L]};
        };

    all_resp = []; all_resp_exc = []; all_resp_inh = []; all_resp_unt = [];
    for jj = 1:length(allsic)
        sic = allsic{jj};
        clear all_gFR all_exc all_inh all_unt
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
            gFR = props1.good_FR_and_untuned;
            gFR = cell_list_op(gFR,[],'or',1);
            all_unt(:,ii) = gFR;
        end
        all_resp = [all_resp all_gFR]; all_resp_exc = [all_resp_exc all_exc]; all_resp_inh = [all_resp_inh all_inh];
        all_resp_unt = [all_resp_unt all_unt];
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
    [within,dvn,xlabels] = make_within_table({'Type','Cond'},[2,length(sic)]);
    dataT = make_between_table({per_active},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Type_by_Cond','hsd'},[1.5 1 1]);
	xdata = make_xdata([length(sic) length(sic)],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;    
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'NB-Exc','NB-Exc-CP','NB-Exc-CA','NB-Inh','NB-Inh-CP','NB-Inh-CP'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 5 10]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.1 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;    
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'NB','NB-C'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 5 10]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.3 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end

%% compare percent response fidelity pooled light on brake vs no-brake and Control Light (from Condition 3 and 5 where there is no light stimulus)
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    allsic = {{[Lb Lbs];[ArL_L];[Ar_L Ars_L]};
    };
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
            gFR = props1.inh; rf = props1.N_Resp_Trials;
            all_inh(:,ii) = mean(exec_fun_on_cell_mat(rf,'mean',gFR),2);
        end
        all_resp = [all_resp all_gFR]; all_resp_exc = [all_resp_exc all_exc]; all_resp_inh = [all_resp_inh all_inh];
    end
    good_FR = [all_resp_exc all_resp_inh];

    [within,dvn,xlabels] = make_within_table({'Type','Cond'},[2,3]);
    dataT = make_between_table({good_FR},dvn);
    ra = RMA(dataT,within);
    ra.ranova
  
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Type_by_Cond','hsd'},[1.5 1 1]);
    h(2)= 0;
	xdata = make_xdata([3 3],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-Exc','NB-Exc','B-Inh','NB-Inh'};
    xticklabels = {'B-Exc','NB-Exc','NB-Exc-C','B-Inh','NB-Inh','NB-Inh-C'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 35 70]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.2 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{'Trials (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end



%% compare zMI pooled light on brake vs no-brake
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    allsic = {{[Lb Lbs];[ArL_L]};
    };
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
            gFR = props1.vals; rf = props1.zMI;
            all_gFR(:,ii) = mean(exec_fun_on_cell_mat(rf,'mean',gFR),2);
            gFR = props1.exc; rf = props1.zMI;
            all_exc(:,ii) = mean(exec_fun_on_cell_mat(rf,'mean',gFR),2);
            gFR = props1.inh; rf = props1.zMI;
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
    h(2)= 0;
	xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-Exc','NB-Exc','B-Inh','NB-Inh'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 35 70]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.5 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{'Trials (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end



%% compare percent responsive cells pooled air on or off brake vs no-brake Unt
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    event_type = {'Air ON','Air OFF','Light ON'};

    allsic = {{[Ab_On Abs_On Ab_Off Abs_Off];[Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off];[Lb Lbs];[ArL_L]};};

    all_resp = []; all_resp_exc = []; all_resp_inh = [];
    all_respV = [];
    for jj = 1:length(allsic)
        sic = allsic{jj};
        clear all_gFR all_exc all_inh all_gV
        for ii = 1:length(sic)
            sit = sic{ii};
            tRs = o.Rs(:,sit);
            props1 = get_props_Rs(tRs,ntrials);
            gFR = props1.vals;
            gFR = cell_list_op(gFR,[],'or',1);
            all_gV(:,ii) = gFR;
            gFR = props1.good_FR_and_untuned;
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
        all_respV = [all_respV all_gV]; 
    end
    good_FR = all_gFR;%[all_resp_exc all_resp_inh];
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
    break;
end

%% compare percent response fidelity pooled air with light brake vs no-brake
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    event_type = {'Air ON','Air OFF','Light ON'};
    allsic = {{[Ab_On Abs_On Ab_Off Abs_Off];[Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off];[Lb Lbs];[ArL_L]};};
    all_resp = []; all_resp_exc = []; all_resp_inh = [];
    for jj = 1:length(allsic)
        sic = allsic{jj};
        clear all_gFR all_exc all_inh
        for ii = 1:length(sic)
            sit = sic{ii};
            tRs = o.Rs(:,sit);
            props1 = get_props_Rs(tRs,ntrials);
            gFR = props1.good_FR_and_untuned; rf = props1.N_Resp_Trials;
            all_gFR(:,ii) = mean(exec_fun_on_cell_mat(rf,'mean',gFR),2);
            gFR = props1.exc; rf = props1.N_Resp_Trials;
            all_exc(:,ii) = mean(exec_fun_on_cell_mat(rf,'mean',gFR),2);
            gFR = props1.inh;rf = props1.N_Resp_Trials;
            all_inh(:,ii) = mean(exec_fun_on_cell_mat(rf,'mean',gFR),2);
        end
        all_resp = [all_resp all_gFR]; all_resp_exc = [all_resp_exc all_exc]; all_resp_inh = [all_resp_inh all_inh];
    end
    good_FR = all_resp;%[all_resp_exc all_resp_inh];

    [within,dvn,xlabels] = make_within_table({'Type','Cond'},[2,2]);
    dataT = make_between_table({good_FR},dvn);
    ra = RMA(dataT,within);
    ra.ranova
  
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Type_by_Cond','hsd'},[1.5 1 1]);
    h(h==1) = 0;
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
    set_axes_limits(gca,[0.35 xdata(end)+.65],[50 75]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-A','NB-A','B-L','NB-L'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[50 70]); xtickangle(45)
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
    xticks = xdata; xticklabels = {'B-Act','NB-Act','B-Sup','NB-Sup'};
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
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[50 75]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-A','NB-A','B-L','NB-L'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[50 70]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.5 -0.04]); %put_axes_labels(gca,{[],[0 0 0]},{{'Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_allf3.pdf',ntrials),600);
    %%
    break;
end


%% compare percent responsive cells pooled
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    event_type = {'Air ON','Air OFF','Light ON'};
%     sic = {[Lb Lbs Ab_On Abs_On];[Ab_Off Abs_Off];ArL_L;[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off]};%;M_On;M_Off};
    allsic = {{[Ab_On Abs_On];[Ar_On ArL_On Ars_On]};
    {[Ab_Off Abs_Off];[Ar_Off ArL_Off Ars_Off]};
    {[Lb Lbs];ArL_L};{M_On};{M_Off}};
%     ind = 3;
%     sic = allsic{1:6};
    all_resp = [];
    for jj = 1:length(allsic)
        sic = allsic{jj};
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
    end
    good_FR = all_resp;
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
    [within,dvn,xlabels] = make_within_table({'Cond'},[size(good_FR,2)]);
    dataT = make_between_table({per_active},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
  
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
	xdata = make_xdata([2 2 2 2],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
    htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-A-On','NB-A-On','B-A-Off','NB-A-Off','B-L-On','NB-L-On','M-On','M-Off'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 15 30]); xtickangle(45)
    changePosition(gca,[0.05 0.01 0.01 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{'Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end


%% compare percent responsive cells pooled air light on vs off with brake no brake
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
%     sic = {[Lb Lbs Ab_On Abs_On];[Ab_Off Abs_Off];ArL_L;[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off]};%;M_On;M_Off};
    allsic = {{[Ab_On Abs_On Ab_Off Abs_Off];[Lb Lbs]};
    {[Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off];[ArL_L]};
    };
%     ind = 3;
%     sic = allsic{1:6};
    all_resp = [];
    for jj = 1:length(allsic)
        sic = allsic{jj};
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
    end
    good_FR = all_resp;
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
    pall_resp = 100*exec_fun_on_cell_mat(all_resp,'sum')./exec_fun_on_cell_mat(all_resp,'length');
    [within,dvn,xlabels] = make_within_table({'Cond','Event'},[2 2]);
    dataT = make_between_table({per_active},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
  
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond_by_Event','hsd'},[1.5 1 1]);
    h(2) = 0;
	xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 1 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.05);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 70]); format_axes(gca);
    xticks = xdata; xticklabels = {'Air','Light','Air','Light'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 25 50]); xtickangle(45)
    changePosition(gca,[0.09 0.01 -0.2 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{'Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);
    hf = get_figure(5,[8 7 1 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(5:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 70]); format_axes(gca);
    xticks = xdata; xticklabels = {'Air','Light','AOn','AOff'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 25 50]); xtickangle(45)
    changePosition(gca,[0.09 0.01 -0.45 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Event','hsd'},[1.5 1 1]);
	xdata = make_xdata([2],[1 1.5]);
    hf = get_figure(5,[8 7 1 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(7:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 70]); format_axes(gca);
    xticks = xdata; xticklabels = {'Air','Light','AOn','AOff'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 25 50]); xtickangle(45)
    changePosition(gca,[0.09 0.01 -0.45 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end


%% Brake vs Brake or No-Brake versus No-Brake for ON OFF comparison
while 1
    all_resp = [];
    for ind = 1
    ntrials = 50;
    event_type = {'Brake','No-Brake'};
%     sic = {[Lb Lbs Ab_On Abs_On];[Ab_Off Abs_Off];ArL_L;[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off]};%;M_On;M_Off};
    allsic = {{[Ab_On Abs_On];[Ab_Off Abs_Off]};
        {[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off]}};%;M_On;M_Off};
%     ind = 3;
    sic = allsic{ind};
    clear all_gFR
    for ii = 1:length(sic)
        sit = sic{ii};
        tRs = o.Rs(:,sit);
        props1 = get_props_Rs(tRs,ntrials);
        gFR = props1.good_FR;
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
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = make_xdata([size(good_FR,2)],[1]);
%     h(h==1) = 0;
    hf = get_figure(5,[8 7 1.5 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.17);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Air ON','Air OFF'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 5 10]); xtickangle(45)
    changePosition(gca,[0.035 0.01 -0.05 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -1 0]});
    pos = get(gca,'Position');
    text(1,ylims(2),sprintf('%s',event_type{ind}),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('unique_cells_across_conditionsSimpBB_%s.pdf',event_type{ind}),600);
    
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
    text(-3.65,2.5,{'Air ON',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(1,1),pmchar,sem_p_gFR_C(1,1))},'FontSize',6,'rotation',0);
    text(3.5,-2.65,{'Air OFF',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(2,2),pmchar,sem_p_gFR_C(2,2))},'FontSize',6,'rotation',0);
    set(HVenn(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); set(HVenn(2),'FaceColor',tcolors{2},'FaceAlpha',0.75);
%     text(-6,6,{'Volunteer',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(3,3),pmchar,sem_p_gFR_C(3,3))},'FontSize',6,'rotation',0);
    ylims = ylim;
    text(1,ylims(2)+0.51,sprintf('%s',event_type{ind}),'FontSize',6);
    text(5,ylims(2)-0.51,sprintf('%.0f%%',IVenn),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('conjunctive_cells_venn_diagramSimpBB_%s.pdf',event_type{ind}),600);
    end
    %%
    break;
end

%% Brake vs No-Brake Air on and Air off Exc and Inh (for VENN Diagrams) and For figure air offset heat map
while 1
    %%
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    event_type = {'Air ON','Air OFF','Light ON'};
    sic = {[Ab_On Abs_On];[Ab_Off Abs_Off];[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off];[Ab_Offc Abs_Offc];[Ar_Offc ArL_Offc Ars_Offc]};
    clear all_gV all_gFR all_exc all_inh;
    for ii = 1:length(sic)
        sit = sic{ii};
        tRs = o.Rs(:,sit);
        props1 = get_props_Rs(tRs,ntrials);
        gFR = props1.vals;
        gFR = cell_list_op(gFR,[],'or',1);
        all_gV(:,ii) = gFR;
        gFR = props1.good_FR_and_untuned;
        gFR = cell_list_op(gFR,[],'or',1);
        all_gFR(:,ii) = gFR;
        gFR = props1.exc; %gFR = props1.good_FR_and_exc;
        gFR = cell_list_op(gFR,[],'or',1);
        all_exc(:,ii) = gFR;
        gFR = props1.inh; %gFR = props1.good_FR_and_inh;
        gFR = cell_list_op(gFR,[],'or',1);
        all_inh(:,ii) = gFR;
    end
    good_FR = [all_exc(:,1),all_inh(:,1),all_exc(:,2),all_inh(:,2),all_exc(:,3),all_inh(:,3),all_exc(:,4),all_inh(:,4)];
%     good_FR = [all_exc(:,1),all_inh(:,1),all_gFR(:,1),all_exc(:,2),all_inh(:,2),all_gFR(:,2),all_exc(:,3),all_inh(:,3),all_gFR(:,3),all_exc(:,4),all_inh(:,4),all_gFR(:,4)];
    
    good_FR_p = cell_list_op(all_exc,all_inh,'or');
%     good_FR_p = cell_list_op(good_FR_p,all_gFR,'or');
    
    inds = 1:2;
    good_FR_pocondt = [cell_list_op(all_exc(:,inds),[],'or',1) cell_list_op(all_inh(:,inds),[],'or',1) cell_list_op(all_gFR(:,inds),[],'or',1)];
    inds = 3:4;
    good_FR_pocondt = [good_FR_pocondt cell_list_op(all_exc(:,inds),[],'or',1) cell_list_op(all_inh(:,inds),[],'or',1) cell_list_op(all_gFR(:,inds),[],'or',1)];
    
    inds = [1:2 5];
    good_FR_pocondt_arb = [cell_list_op(all_exc(:,inds),[],'or',1) cell_list_op(all_inh(:,inds),[],'or',1) cell_list_op(all_gFR(:,inds),[],'or',1)];
    inds = [3:4 6];
    good_FR_pocondt_arb = [good_FR_pocondt_arb cell_list_op(all_exc(:,inds),[],'or',1) cell_list_op(all_inh(:,inds),[],'or',1) cell_list_op(all_gFR(:,inds),[],'or',1)];
    
    good_FR_bnb(:,1) = cell_list_op(good_FR_p(:,1:2),[],'and',1);
    good_FR_bnb(:,2) = cell_list_op(good_FR_p(:,3:4),[],'and',1);
    all_good_FR = cell_list_op(good_FR_bnb,[],'or'); 
    p_agfr = 100*exec_fun_on_cell_mat(all_good_FR,'sum')./exec_fun_on_cell_mat(all_good_FR,'length');
    [mpagfr,sempagfr] = findMeanAndStandardError(p_agfr(:,1));
    
    all_good_bnb1(:,1) = cell_list_op(good_FR_pocondt(:,1:3),[],'or',1); all_good_bnb1(:,2) = cell_list_op(good_FR_pocondt(:,4:6),[],'or',1);
    all_good_FR = cell_list_op(all_good_bnb1,[],'or');
    p_agfr = 100*exec_fun_on_cell_mat(all_good_FR,'sum')./exec_fun_on_cell_mat(all_good_FR,'length');
    [mpagfr1,sempagfr1] = findMeanAndStandardError(p_agfr(:,1));
    disp('Done');
    %%
    avar = [all_gV(:,[1 3]) all_gV(:,[2 4]) all_gV(:,[5 6])];
    pavar = find_percent(avar);
    [within,dvn,xlabels] = make_within_table({'Type','Cond'},[3,2]);
    dataT = make_between_table({pavar},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    print_for_manuscript(ra)
    %%
    % for different VENN Diagrams change good_FRV = all_exc(:,[1 3]) or
    % all_exc(:,[2 4]) or all_exc(:,[5 6]) same thing for all_inh
    good_FRV = [cell_list_op(all_gV(:,[1:2 5]),[],'or',1) cell_list_op(all_gV(:,[3:4 6]),[],'or',1)]; 
    cell_any = descriptiveStatistics(find_percent(cell_list_op(good_FRV,[],'or',1)));
    tcolors = mData.colors;
    hf = get_figure(6,[10 7 1.5 1]);
%     good_FRV = good_FR_bnb; %good_FRV = all_inh(:,[5 6]);
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
    xticks = xdata; xticklabels = {'Brake','Conjunc','No-Brake'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.09 0.01 -0.3 0]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all3.pdf',ntrials),600);
    %% for heat maps
    good_FR = [all_exc(:,1),all_inh(:,1),all_exc(:,2),all_inh(:,2),all_exc(:,3),all_inh(:,3),all_exc(:,4),all_inh(:,4)]; % for comparison of air onset with air offset across brake and no-brake
    good_FR = [all_exc(:,1),all_inh(:,1),all_exc(:,2),all_inh(:,2),all_exc(:,5),all_inh(:,5),all_exc(:,3),all_inh(:,3),all_exc(:,4),all_inh(:,4),all_exc(:,6),all_inh(:,6)]; % for comparison of air onset offset and arb
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(good_FR,0.5,0.05);
    mOI = mCI; semOI = semCI;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
%     txl = {'B-A-On-Exc','B-A-On-Inh','B-A-On-Unt','B-A-Off-Exc','B-A-Off-Inh','B-A-Off-Unt','NB-A-On-Exc','NB-A-On-Inh','NB-A-On-Unt','NB-A-Off-Exc','NB-A-Off-Inh','NB-A-Off-Unt'};%,'M-On','M-Off'};
    txl = {'Exc','Inh'}; txl = repmat(txl,1,4);
    txl = {'Exc','Inh'}; txl = repmat(txl,1,6);
%     txl = {'Air-On','Air-Off'};txl = repmat(txl,1,2);
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.25 1.25]);
    hf = get_figure(5,[8 7 1.75 1.75]);
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
    set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[-0.02 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',ntrials),600);
    
    %% with NaN
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(good_FR,0.5,0.05);
    mOI = mCI; semOI = semCI;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
%     txl = {'B-A-On-Exc','B-A-On-Inh','B-A-On-Unt','B-A-Off-Exc','B-A-Off-Inh','B-A-Off-Unt','NB-A-On-Exc','NB-A-On-Inh','NB-A-On-Unt','NB-A-Off-Exc','NB-A-Off-Inh','NB-A-Off-Unt'};%,'M-On','M-Off'};
    txl = {'Exc','Inh','Unt'}; txl = repmat(txl,1,4);
    txl = {'Exc','Inh'}; txl = repmat(txl,1,4);
%     txl = {'Air-On','Air-Off'};txl = repmat(txl,1,2);
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.25 1.25]);
    hf = get_figure(5,[8 7 1.75 1.75]);
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
    changePosition(gca,[-0.02 0 -0.04 0]);
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
    set(hf,'Position',[7 3 2 1]);
%     set(hf,'Position',[7 3 2 1]);
    set(H,'linewidth',0.5);
    txl = {' B-AOn-Exc',' B-AOn-Inh',' B-AOff-Exc',' B-AOff-Inh','NB-AOn-Exc','NB-AOn-Inh','NB-AOff-Exc','NB-AOff-Inh'};
    txl = {' B-AOn-Exc',' B-AOn-Inh',' B-AOff-Exc',' B-AOff-Inh',' B-Arb-Exc',' B-Arb-Inh','NB-AOn-Exc','NB-AOn-Inh','NB-AOff-Exc','NB-AOff-Inh','NB-Arb-Exc','NB-Arb-Inh'};
%     txl = {' Tun-B-A',' Tun-B-L','Tun-NB-A','Tun-NB-L',' Unt-B-A',' Unt-B-L','Unt-NB-A','Unt-NB-L'};
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel({'Eucledian','Distance'});changePosition(hx,[0 -0.3 0]);
    changePosition(gca,[0.07 0.0 0.0 -0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);
    
    %%
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(good_FR_p,0.5,0.05);
    for an = 1:5
        selA1(an) = all_CI{an}(1,2);
        selA2(an) = all_CI{an}(4,3);
    end
%     [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(all_gFR,0.5,0.05);
    mOI = mCI; semOI = semCI;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
%     mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = 0;%min([mOI(:);semOI(:)]);
    
%     mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
%     txl = {'B-A-On-Exc','B-A-On-Inh','B-A-On-Unt','B-A-Off-Exc','B-A-Off-Inh','B-A-Off-Unt','NB-A-On-Exc','NB-A-On-Inh','NB-A-On-Unt','NB-A-Off-Exc','NB-A-Off-Inh','NB-A-Off-Unt'};%,'M-On','M-Off'};
    txl = {'Air-On','Air-Off'};txl = repmat(txl,1,2);
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
%     imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.25 1.25]);
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
    changePosition(gca,[0.0 0 -0.05 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.07]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',ntrials),600);
    %%
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(good_FR_pocondt,0.5,0.05);
    mOI = mCI; semOI = semCI;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
%     mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = 0;%min([mOI(:);semOI(:)]);
    
%     mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
%     txl = {'B-A-On-Exc','B-A-On-Inh','B-A-On-Unt','B-A-Off-Exc','B-A-Off-Inh','B-A-Off-Unt','NB-A-On-Exc','NB-A-On-Inh','NB-A-On-Unt','NB-A-Off-Exc','NB-A-Off-Inh','NB-A-Off-Unt'};%,'M-On','M-Off'};
    txl = {'Exc','Inh','Unt'};txl = repmat(txl,1,2);
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
%     imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.5 1.5]);
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
    changePosition(gca,[0.01 0 -0.05 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.07]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',ntrials),600);
    %%
    break;
end

%% Brake vs No-Brake Air on Exc and Inh (for drawing VENN diagrams brake vs no-brake)
while 1
    tcolors = mData.colors;
    hf = get_figure(6,[10 7 1.5 1]);
    good_FRV = good_FR(:,1:2);% good_FRV = good_FR(:,3:4);
    good_FRV = all_resp_exc;%all_resp_unt;
    [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(good_FRV,0.5,0.05);
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
    save_pdf(hf,mData.pdf_folder,sprintf('conjunctive_cells_venn_diagram_%s.pdf',''),600);
    %%
    resp = all_resp;
   [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = {'L-b','L-nb','AOn-b','AOff-b','AOn-nb','AOff-nb'};%,'M-On','M-Off'};
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.5 1.5]);
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
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',ntrials),600);
    %%
    break;
end

%% Brake vs No-Brake and Voluntary Motion On Off
while 1
    all_resp = [];
    for ind = 3
    ntrials = 50;
    event_type = {'Air','Light','All'};
    allsic = {{[Ab_On Abs_On Ab_Off Abs_Off];[Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off];[M_On M_Off]};
    {[Lb Lbs];ArL_L;[M_On M_Off]};
    {[Ab_On Abs_On Ab_Off Abs_Off Lb Lbs];[Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off ArL_L];[M_On M_Off]}};
%     ind = 3;
    sic = allsic{ind};
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
    hf = get_figure(5,[8 7 1.5 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors(1:6);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.17);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Brake','No-Brake','Motion'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.035 0.01 -0.25 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
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
    save_pdf(hf,mData.pdf_folder,sprintf('conjunctive_cells_venn_diagram_%s_volunM.pdf',event_type{ind}),600);
    end
    %%
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(all_gFR,0.5,0.05);
    
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
%     h(h==1) = 0;
    hf = get_figure(5,[8 7 1.5 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(4:6);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.17);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-NB','B-M','NB-M'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.035 0.01 -0.25 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    pos = get(gca,'Position');
%     text(1,ylims(2),sprintf('%s',event_type{ind}),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('unique_cells_across_conditions_%s_volunM.pdf',event_type{ind}),600);
    
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


%% Brake vs No-Brake Air on and Air off and light ON Exc and Inh
while 1
    %%
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    sic = {[Ab_On Abs_On];[Ab_Off Abs_Off];[Lb Lbs];[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off];[ArL_L]};
    clear all_gV all_gFR all_exc all_inh all_good_FR;
    for ii = 1:length(sic)
        sit = sic{ii};
        tRs = o.Rs(:,sit);
        props1 = get_props_Rs(tRs,ntrials);
        gFR = props1.vals;
        gFR = cell_list_op(gFR,[],'or',1);
        all_gV(:,ii) = gFR;
        gFR = props1.good_FR_and_untuned;
        gFR = cell_list_op(gFR,[],'or',1);
        all_gFR(:,ii) = gFR;
        gFR = props1.exc; %gFR = props1.good_FR_and_exc;
        gFR = cell_list_op(gFR,[],'or',1);
        all_exc(:,ii) = gFR;
        gFR = props1.inh; %gFR = props1.good_FR_and_inh;
        gFR = cell_list_op(gFR,[],'or',1);
        all_inh(:,ii) = gFR;
%         gFR = cell_list_op(props1.good_FR,props1.vals,'or'); %gFR = props1.good_FR_and_inh;
%         gFR = cell_list_op(gFR,[],'or',1);
%         all_good_FR(:,ii) = gFR;
        
    end
    good_FR = [all_exc(:,1),all_inh(:,1),all_exc(:,2),all_inh(:,2),all_exc(:,3),all_inh(:,3),all_exc(:,4),all_inh(:,4)];
    good_FR = [all_exc(:,1),all_inh(:,1),all_gFR(:,1),all_exc(:,2),all_inh(:,2),all_gFR(:,2),all_exc(:,3),all_inh(:,3),all_gFR(:,3),all_exc(:,4),all_inh(:,4),all_gFR(:,4)];
    
    good_FR_AL = [cell_list_op(all_gV(:,1:2),[],'or',1) all_gV(:,3) cell_list_op(all_gV(:,4:5),[],'or',1) all_gV(:,6)];
    
    good_FR_p = cell_list_op(all_exc,all_inh,'or');
%     good_FR_p = cell_list_op(good_FR_p,all_gFR,'or');
    
    inds = 1:2;
    good_FR_pocondt = [cell_list_op(all_exc(:,inds),[],'or',1) cell_list_op(all_inh(:,inds),[],'or',1) cell_list_op(all_gFR(:,inds),[],'or',1)];
    inds = 3:4;
    good_FR_pocondt = [good_FR_pocondt cell_list_op(all_exc(:,inds),[],'or',1) cell_list_op(all_inh(:,inds),[],'or',1) cell_list_op(all_gFR(:,inds),[],'or',1)];
    
   clear all_good_FR
    all_good_BNB = [cell_list_op(good_FR_AL(:,1:2),[],'or',1) cell_list_op(good_FR_AL(:,3:4),[],'or',1)];
    all_good_FR = cell_list_op(all_good_BNB,[],'or',1);
    p_agfr = 100*exec_fun_on_cell_mat(all_good_FR,'sum')./exec_fun_on_cell_mat(all_good_FR,'length');
    [mpagfr1,sempagfr1] = findMeanAndStandardError(p_agfr(:,1));
    disp('Done');
    
    %%
    tcolors = mData.colors;
    hf = get_figure(6,[10 7 1.5 1]);
    good_FRV = all_good_BNB;
    [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(good_FRV,0.5,0.05);
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
    save_pdf(hf,mData.pdf_folder,sprintf('conjunctive_cells_venn_diagram_%s.pdf',''),600);
  
    %%
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(good_FR,0.5,0.05);
    mOI = mCI; semOI = semCI;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
%     txl = {'B-A-On-Exc','B-A-On-Inh','B-A-On-Unt','B-A-Off-Exc','B-A-Off-Inh','B-A-Off-Unt','NB-A-On-Exc','NB-A-On-Inh','NB-A-On-Unt','NB-A-Off-Exc','NB-A-Off-Inh','NB-A-Off-Unt'};%,'M-On','M-Off'};
    txl = {'Exc','Inh','Unt'}; txl = repmat(txl,1,4);
    txl = {'Exc','Inh'}; txl = repmat(txl,1,6);
%     txl = {'Air-On','Air-Off'};txl = repmat(txl,1,2);
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.25 1.25]);
    hf = get_figure(5,[8 7 1.75 1.75]);
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
    changePosition(gca,[-0.02 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',ntrials),600);
    %%
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(good_FR_AL,0.5,0.05);
    for an = 1:5
        selA1(an) = all_CI{an}(1,2);
        selA2(an) = all_CI{an}(4,3);
    end
%     [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(all_gFR,0.5,0.05);
    mOI = mCI; semOI = semCI;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = 0;%min([mOI(:);semOI(:)]);
    
%     mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
%     txl = {'B-A-On-Exc','B-A-On-Inh','B-A-On-Unt','B-A-Off-Exc','B-A-Off-Inh','B-A-Off-Unt','NB-A-On-Exc','NB-A-On-Inh','NB-A-On-Unt','NB-A-Off-Exc','NB-A-Off-Inh','NB-A-Off-Unt'};%,'M-On','M-Off'};
    txl = {'Air','Light'};txl = repmat(txl,1,2);
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1 1]);
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
    changePosition(gca,[0.0 0 -0.08 -0.02]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.07]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',ntrials),600);

    %%
    break;
end

