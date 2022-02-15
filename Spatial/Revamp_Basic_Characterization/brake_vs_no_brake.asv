function brake_vs_no_brake

%% Load Data
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei = evalin('base','ei'); 
    
    selContexts = [1 4 6 2 7 2 7 3 4 5 3 4 5 0 0 3 4 5 3 4 5];
    rasterNames = {'light22T','light22T','light22T','airOnsets22T','airOnsets22T','airOffsets22T',...
                    'airOffsets22T','airOnsets22T','airOnsets22T','airOnsets22T',...
                    'airOffsets22T','airOffsets22T','airOffsets22T','motionOnsets','motionOffsets',...
                    'airD','airD','airD','airIT','airIT','airIT'};
%     rasterNamesTxt = {'Lb','ArL-L','Lb*','Ab-On','Ab*-On','Ab-Off','Ab*-Off',...
%         'Ar-On','ArL-On','Ar*-On','Ar-Off','ArL-Off','Ar*-Off','M-On','M-Off',...
%         'Ar-D','ArL-D','Ar*-D','Ar-T','ArL-T','Ar*-T'};
    rasterNamesTxt = {'B-L','NB-L','B-L*','B-A-On','B-A*-On','B-A-Off','B-A*-Off',...
        'NB-A-On','NB-AL-On','NB-A*-On','NB-A-Off','NB-AL-Off','NB-A*-Off','M-On','M-Off',...
        'NB-A-D','NB-AL-D','NB-A*-D','NB-A-T','NB-AL-T','NB-A*-T'};
    o = get_data1(ei,selContexts,rasterNames);
    
    Lb = 1; ArL_L = 2; Lbs = 3; Ab_On = 4; Abs_On = 5; Ab_Off = 6; Abs_Off = 7; 
    Ar_On = 8; ArL_On = 9; Ars_On = 10; Ar_Off = 11; ArL_Off = 12; Ars_Off = 13;
    M_On = 14; M_Off = 15; Ar_D = 16; ArL_D = 17; Ars_D = 18; Ar_T = 19; ArL_T = 20; Ars_T = 21;
    break
end
n = 0;

%% combine distance time rasters 
siD = [Ar_D ArL_D Ars_D]; siT = [Ar_T ArL_T Ars_T];
RsD = o.Rs(:,siD);mRD = o.mR(:,siD); RsT = o.Rs(:,siT);mRT = o.mR(:,siT);
respDT = combine_distance_time_rasters(RsD,RsT);

%% Show sample rasters
while 1
    rtype = {'B_L','NB_L','B_AOn','B_AOff','NB_AOn','NB_AOff','NB-A'};
    cellNums = {[92 168 387 436];
        [70 59 328 558];
        [140 39 17 66 193 140];
        [575 581 373 532];
        [142 66 54 567 429 235 164];
        [439 42 215 353];
        [552 567 165 103]};
    an = 1; cn = 8;
    ntrials = 50;
    asi = [Lb ArL_L Ab_On Abs_Off Ar_On Ar_Off Ar_D Ar_T];
    si = asi(cn);
%     si = [Lb];
    Rs = o.Rs(:,si);
    props1 = get_props_Rs(Rs,ntrials);
%     plotRasters_simplest(Rs{an},find(props1.vals{an}))
    plotRasters_simplest(Rs{an},find(props1.good_zMI{an}))
    break;
    % find(resp_valsC{an}(:,cn));
    R = Rs{an};
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.15 0.04],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
        [-60 -475]);
    set(gcf,'color','w'); set(gcf,'Position',[10 4 3.4 1]);
    plot_time_rasters(R,cellNums{cn},ff);
    for ii = 1:4
        set(ff.h_axes(1,ii),'xtick',[1 18 36],'xticklabels',{'-2','0','2'});
    end
    colormap_ig
    save_pdf(ff.hf,mData.pdf_folder,sprintf('rasters_%s',rtype{cn}),600);
    break;
end

%% Show sample rasters distance
while 1
    rtype = {'B_L','NB_L','B_AOn','B_AOff','NB_AOn','NB_AOff','NB-A'};
    cellNums = {[92 168 387 436];
        [70 59 328 558];
        [140 39 17 66 193 140];
        [575 581 373 532];
        [142 66 54 567 429 235 164];
        [439 42 215 353];
        [11 54 116 165]};
    an = 1; cn = 7;
    ntrials = 50;
    asi = [Lb ArL_L Ab_On Abs_Off Ar_On Ar_Off Ar_D Ar_T];
    si = asi(cn);
%     si = [Lb];
    Rs = o.Rs(:,si);
    props1 = get_props_Rs(Rs,ntrials);
%     plotRasters_simplest(Rs{an},find(props1.vals{an}))
%     plotRasters_simplest(Rs{an},find(props1.good_zMI{an}))
%     break;
    % find(resp_valsC{an}(:,cn));
    R = Rs{an};
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.15 0.04],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
        [-60 -475]);
    set(gcf,'color','w'); set(gcf,'Position',[10 4 3.4 1]);
    plot_dist_rasters(R,cellNums{cn},ff);
    for ii = 1:4
        set(ff.h_axes(1,ii),'xtick',[1 24.5 49],'xticklabels',{'0','75','150'});
    end
    colormap_ig
    save_pdf(ff.hf,mData.pdf_folder,sprintf('rasters_%s',rtype{cn}),600);
    break;
end

%% Show sample rasters time
while 1
    rtype = {'B_L','NB_L','B_AOn','B_AOff','NB_AOn','NB_AOff','NB-A','NB-A'};
    cellNums = {[92 168 387 436];
        [70 59 328 558];
        [140 39 17 66 193 140];
        [575 581 373 532];
        [142 66 54 567 429 235 164];
        [439 42 215 353];
        [11 54 116 165];
        [328 565 390 406]};
    an = 1; cn = 8;
    ntrials = 50;
    asi = [Lb ArL_L Ab_On Abs_Off Ar_On Ar_Off Ar_D Ar_T];
    si = asi(cn);
%     si = [Lb];
    Rs = o.Rs(:,si);
    props1 = get_props_Rs(Rs,ntrials);
%     plotRasters_simplest(Rs{an},find(props1.vals{an}))
%     plotRasters_simplest(Rs{an},find(props1.good_zMI{an}))
%     break;
    % find(resp_valsC{an}(:,cn));
    R = Rs{an};
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.15 0.04],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
        [-60 -475]);
    set(gcf,'color','w'); set(gcf,'Position',[10 4 3.4 1]);
    plot_time_rasters(R,cellNums{cn},ff);
    for ii = 1:4
        set(ff.h_axes(1,ii),'xtick',[1 68 136],'xticklabels',{'0','7.5','15'});
    end
    colormap_ig
    save_pdf(ff.hf,mData.pdf_folder,sprintf('rasters_%s',rtype{cn}),600);
    break;
end


%% compare the zMIs
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    si = [Lb Lbs Ab_On Abs_On Ab_Off Abs_Off ArL_L Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off M_On M_Off];
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    good_FR = props1.good_FR;
    all_zMIs = props1.zMI;
    zMIs = [];
    for rr = 1:size(all_zMIs,1)
        for cc = 1:size(all_zMIs,2)
            resp = good_FR{rr,cc};
            tzmis = all_zMIs{rr,cc};
            zMIs(rr,cc) = nanmean(tzmis(resp));
        end
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},[15]);
    dataT = make_between_table({zMIs},dvn);
    ra = RMA(dataT,within);
    ra.ranova
%     ra.mauchly
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = make_xdata([2 2 2 1 3 3 2],[1 1.5]);
    ptab = 0;
    if ptab h(h==1) = 0; end
    hf = get_figure(5,[8 7 1.75 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    

    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'yspacing',0.5);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rasterNamesTxt(si);
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 round(maxY/2) round(maxY)]); xtickangle(45);
    if ptab
    changePosition(gca,[0.08 0.01 0.0 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. zMI',[-1 3 0]});
    ha = gca; ptable = extras.pvalsTable;
    display_p_table(ha,hbs,[0 -0.07 0 0.9],ptable);
    else
    changePosition(gca,[0.01 -0.01 0.06 0.01]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. zMI',[0 0 0]});
    end
    save_pdf(hf,mData.pdf_folder,sprintf('zMIs_good_FR_all_Conditions.pdf'),600);
    %%
    break;
end

%% compare percent responsive cells
while 1
    ntrials = 40; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    si = [Lb Lbs Ab_On Abs_On Ab_Off Abs_Off ArL_L Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off M_On M_Off];
    props1 = get_props_Rs(o.Rs,ntrials);
    good_FR = props1.good_FR_and_tuned(:,si);
    good_FR_any = cell_list_op(good_FR,[],'or');
    good_FR_all = cell_list_op(good_FR,[],'and');
    per_active =[]; 
    for rr = 1:size(good_FR,1)
        for cc = 1:size(good_FR,2)
            tts = good_FR{rr,cc};
            per_active(rr,cc) = 100*sum(tts)/length(tts);
        end
        tts = good_FR_any{rr,1};        per_active_any(rr) = 100*sum(tts)/length(tts);
        tts = good_FR_all{rr,1};        per_active_all(rr) = 100*sum(tts)/length(tts);
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},[length(si)]);
    dataT = make_between_table({per_active},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = make_xdata([2 2 2 1 3 3 2],[1 1.5]);
    h(h==1) = 0;
    hf = get_figure(6,[8 7 2.25 2]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors;
    [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.5);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rasterNamesTxt(si);
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 5 10]); xtickangle(45);
%     changePosition(gca,[0.08 0.01 0.0 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{'Responsive Cells (%)',[-0.7 +15 0]});
    changePosition(gca,[0.06 0.01 0.02 -0.55]); put_axes_labels(gca,{[],[0 0 0]},{'Responsive Cells (%)',[-1.25 -5 0]});
    ha = gca; ptable = extras.pvalsTable;
    ha = display_p_table_img(ha,hbs,[0 0.23 0 0.37],ptable); %ytickangle(20)
%     htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    
%     [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
% %     h(h==1) = 0;
%     hf = get_figure(5,[8 7 2 1.5]);
%     % s = generate_shades(length(bins)-1);
%     tcolors = colors;
%     [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
%         'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
%         'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
%     maxY = maxY + 5;
%     ylims = ylim;
%     format_axes(gca);
%     htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
% %     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
%     set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
%     xticks = xdata; xticklabels = rasterNamesTxt(si);
%     set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 50 100]); xtickangle(45)
%     changePosition(gca,[0.01 0.01 0.05 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Responsive Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end

%% compare percent responsive cells clustering based
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    si = [Lb Lbs Ab_On Abs_On Ab_Off Abs_Off ArL_L Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off M_On M_Off];
    props1 = get_props_Rs(o.Rs,ntrials);
    good_FR = props1.clus_based(:,si);
    good_FR_any = cell_list_op(good_FR,[],'or');
    good_FR_all = cell_list_op(good_FR,[],'and');
    per_active =[]; 
    for rr = 1:size(good_FR,1)
        for cc = 1:size(good_FR,2)
            tts = good_FR{rr,cc};
            per_active(rr,cc) = 100*sum(tts)/length(tts);
        end
        tts = good_FR_any{rr,1};        per_active_any(rr) = 100*sum(tts)/length(tts);
        tts = good_FR_all{rr,1};        per_active_all(rr) = 100*sum(tts)/length(tts);
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},[length(si)]);
    dataT = make_between_table({per_active},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = make_xdata([2 2 2 1 3 3 2],[1 1.5]);
    h(h==1) = 0;
    hf = get_figure(5,[8 7 2.25 2]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors;
    [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.5);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rasterNamesTxt(si);
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 20 40]); xtickangle(45);
%     changePosition(gca,[0.08 0.01 0.0 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{'Responsive Cells (%)',[-0.7 +15 0]});
    changePosition(gca,[0.06 0.01 0.02 -0.55]); put_axes_labels(gca,{[],[0 0 0]},{'Responsive Cells (%)',[-1.25 -5 0]});
    ha = gca; ptable = extras.pvalsTable;
    ha = display_p_table_img(ha,hbs,[0 0.23 0 0.37],ptable); %ytickangle(20)
%     htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    
%     [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
% %     h(h==1) = 0;
%     hf = get_figure(5,[8 7 2 1.5]);
%     % s = generate_shades(length(bins)-1);
%     tcolors = colors;
%     [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
%         'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
%         'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
%     maxY = maxY + 5;
%     ylims = ylim;
%     format_axes(gca);
%     htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
% %     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
%     set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
%     xticks = xdata; xticklabels = rasterNamesTxt(si);
%     set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 50 100]); xtickangle(45)
%     changePosition(gca,[0.01 0.01 0.05 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Responsive Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end

%% compare percent responsive cells pooled light ON brake vs no-brake exc and inh
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
  
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Type_by_Cond','hsd'},[1.5 1 1]);
	xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-Act','NB-Act','B-Sup','NB-Sup'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 5 10]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.4 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{'Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end

%% compare percent response fidelity pooled light on brake vs no-brake
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
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',30,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-Act','NB-Act','B-Sup','NB-Sup'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 35 70]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.3 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{'Response','Fidelity (% trials)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
    %%
    break;
end

%% compare percent responsive cells pooled air on or off brake vs no-brake exc and inh
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    event_type = {'Air ON','Air OFF','Light ON'};
%     sic = {[Lb Lbs Ab_On Abs_On];[Ab_Off Abs_Off];ArL_L;[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off]};%;M_On;M_Off};
    allsic = {{[Ab_On Abs_On];[Ar_On ArL_On Ars_On]};
    };
    allsic = {{[Ab_Off Abs_Off];[Ar_Off ArL_Off Ars_Off]};
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
  
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Type_by_Cond','hsd'},[1.5 1 1]);
	xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-Act','NB-Act','B-Sup','NB-Sup'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 15 30]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.4 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{'Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
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
    allsic = {{[Ab_Off Abs_Off];[Ar_Off ArL_Off Ars_Off]};
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
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',30,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-Act','NB-Act','B-Sup','NB-Sup'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 35 70]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.3 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{'Response','Fidelity (% trials)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all.pdf',ntrials),600);
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


%% compare percent silent cells or active cells
while 1
    si = [Lb Lbs Ab_On Abs_On Ab_Off Abs_Off ArL_L Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off M_On M_Off];
    props1 = get_props_Rs(o.Rs,50); 
    good_FR = props1.N_Resp_Trials(:,si);
    per_silent =[]; per_active = [];
    for rr = 1:size(good_FR,1)
        for cc = 1:size(good_FR,2)
            tts = good_FR{rr,cc};
            per_silent(rr,cc) = 100*sum(tts == 0)/length(tts);
            per_active(rr,cc) = 100*sum(tts >= 50)/length(tts);
        end
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},[size(good_FR,2)]);
    dataT = make_between_table({per_silent},dvn);
%     dataT = make_between_table({per_active},dvn);
    ra = RMA(dataT,within,{'lsd','hsd'});
    ra.ranova
%     ra.mauchly

    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = make_xdata([2 2 2 1 3 3 2],[1 1.5]);
    h(h==1) = 0;
    hf = get_figure(5,[8 7 2.2 1.5]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.5);
    
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rasterNamesTxt(si);
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 50 100]); xtickangle(30);
    changePosition(gca,[0.06 0.01 0.03 -0.55]); put_axes_labels(gca,{[],[0 0 0]},{'Silent Cells (%)',[-1.25 -10 0]});
    ha = gca; ptable = extras.pvalsTable;
    display_p_table_img(ha,hbs,[0 0.23 0 0.37],ptable);%ytickangle(10);
    save_pdf(hf,mData.pdf_folder,sprintf('silent_cells_across_conditions.pdf'),600);
    %%
    break;
end


%% percent unique cells
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    si = [Lb Lbs Ab_On Abs_On Ab_Off Abs_Off ArL_L Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off M_On M_Off];
%     si = [Ar_Off ArL_Off Ars_Off];
    props1 = get_props_Rs(o.Rs,ntrials);
    good_FR = props1.good_FR(:,si);
    good_FR = props1.clus_based(:,si);
    for ii = 1:5
        good_FRmat = cell2mat(good_FR(ii,:)); cgfr = corr(good_FRmat); cgfr(cgfr==1) = NaN;
        acgfr(:,:,ii) = cgfr;
    end
    mcgfr = mean(acgfr,3);
    figure(100);clf;imagesc(mcgfr);colorbar;set(gca,'Ydir','normal')
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
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1.5 1 1]);
    h(h==1)=0;
    xdata = make_xdata([2 2 2 1 3 3 2],[1 1.5]);
    hf = get_figure(5,[8 7 2 1.5]);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
    htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
%     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rasterNamesTxt(si);
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 5 10]); xtickangle(45)
    changePosition(gca,[0.01 0.01 0.05 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Responsive Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('unique_cells_across_conditions_%d_all.pdf',ntrials),600);
    
    break;
end

%% compare the firing rate
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    si = [Lb Lbs Ab_On Abs_On Ab_Off Abs_Off ArL_L Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off M_On M_Off];
%     si = [Ar_t_D Ar_i_T ArL_i_T ArL_t_D Ars_t_D Ars_i_T];
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    good_FR = props1.good_FR;
    all_zMIs = props1.mean_FR;% all_zMIs = props1.max_FR;
    zMIs = [];
    for rr = 1:size(all_zMIs,1)
        for cc = 1:size(all_zMIs,2)
            resp = good_FR{rr,cc};
            tzmis = all_zMIs{rr,cc};
            zMIs(rr,cc) = nanmean(tzmis(resp));
        end
    end
    [within,dvn,xlabels] = make_within_table({'Cond'},[15]);
    dataT = make_between_table({zMIs},dvn);
    ra = RMA(dataT,within);
    ra.ranova
%     ra.mauchly
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = make_xdata([2 2 2 1 3 3 2],[1 1.5]);
    ptab = 0;
    if ptab h(h==1) = 0; end
    hf = get_figure(5,[8 7 1.75 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    

    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'yspacing',0.5);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rasterNamesTxt(si);
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    if ptab
    changePosition(gca,[0.08 0.01 0.0 -0.5]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. FR (AU)',[-1 3 0]});
    ha = gca; ptable = extras.pvalsTable;
    display_p_table(ha,hbs,[0 -0.07 0 0.9],ptable);
    else
    changePosition(gca,[0.01 -0.01 0.06 0.01]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. FR (AU)',[0 0 0]});
    end
    save_pdf(hf,mData.pdf_folder,sprintf('FR_all_Conditions.pdf'),600);
    %%
    break;
end

%% Lb versus AbOn versus AbOff
while 1
    all_resp = [];
    
    ntrials = 50;
%     sic = {[Lb Lbs Ab_On Abs_On];[Ab_Off Abs_Off];ArL_L;[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off]};%;M_On;M_Off};
    sic = {Lb;Ab_On};
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
    xticks = xdata; xticklabels = {'Brake','No-Brake','IT-Interval'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.035 0.01 -0.05 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -5 0]});
    pos = get(gca,'Position');
    text(1,ylims(2),sprintf('%s',event_type{ind}),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('unique_cells_across_conditions_%s.pdf',event_type{ind}),600);
    
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
    text(-3.65,2.5,{'Brake',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(1,1),pmchar,sem_p_gFR_C(1,1))},'FontSize',6,'rotation',0);
    text(4.5,-2.65,{'No-Brake',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(2,2),pmchar,sem_p_gFR_C(2,2))},'FontSize',6,'rotation',0);
    set(HVenn(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); set(HVenn(2),'FaceColor',tcolors{2},'FaceAlpha',0.75);
%     text(-6,6,{'Volunteer',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(3,3),pmchar,sem_p_gFR_C(3,3))},'FontSize',6,'rotation',0);
    ylims = ylim;
    text(1,ylims(2)+0.51,sprintf('%s',event_type{ind}),'FontSize',6);
    text(5,ylims(2)-0.51,sprintf('%.0f%%',IVenn),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('conjunctive_cells_venn_diagram_%s.pdf',event_type{ind}),600);
    
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

%% Brake (before-after vs No-Brake Air ON
while 1
    all_resp = [];
    for ind = 1
    ntrials = 50;
    event_type = {'Air ON'};
%     sic = {[Lb Lbs Ab_On Abs_On];[Ab_Off Abs_Off];ArL_L;[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off]};%;M_On;M_Off};
    allsic = {{[Ab_On Abs_On];[Ar_On ArL_On Ars_On];}};
    
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
    xticks = xdata; xticklabels = {'Brake-B','No-Brake','Brake-A'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.035 0.01 -0.05 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -5 0]});
    pos = get(gca,'Position');
    text(1,ylims(2),sprintf('%s',event_type{ind}),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('unique_cells_across_conditions_%s.pdf',event_type{ind}),600);
    
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
    
    gFR_C_all = cell_list_op(good_FRV(:,1),good_FRV(:,2),'and',1);    %gFR_C_all = cell_list_op(gFR_C_all,good_FRV(:,3),'and',1);
    p_gFR_C_all = mean(100*exec_fun_on_cell_mat(gFR_C_all,'sum')./exec_fun_on_cell_mat(gFR_C_all,'length'));
    [mp_gFR_C_all,semp_gFR_C_all] = findMeanAndStandardError(100*exec_fun_on_cell_mat(gFR_C_all,'sum')./exec_fun_on_cell_mat(gFR_C_all,'length'));
    AVenn = [p_gFR_C(1,1) p_gFR_C(2,2)]; IVenn = [p_gFR_C(1,2)];
%     [HVenn] = venn(AVenn,IVenn,'ErrMinMode','None');
%     AVenn = [p_gFR_C(1,1) p_gFR_C(2,2) p_gFR_C(3,3)]; IVenn = [p_gFR_C(1,2) p_gFR_C(1,3) p_gFR_C(2,3) p_gFR_C_all];
    [HVenn] = venn(AVenn,IVenn,'ErrMinMode','None');
    format_axes(gca);
    axis equal; axis off;
    changePosition(gca,[0.0 -0.0 -0.05 -0.05]); %put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -5 0]});
    pmchar=char(177);
    text(-3.65,2.5,{'Brake-B',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(1,1),pmchar,sem_p_gFR_C(1,1))},'FontSize',6,'rotation',0);
    text(4.5,-2.65,{'No-Brake',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(2,2),pmchar,sem_p_gFR_C(2,2))},'FontSize',6,'rotation',0);
%     text(0,4,{sprintf('Brake-A (%0.0f%c%0.0f%%)',m_p_gFR_C(3,3),pmchar,sem_p_gFR_C(3,3))},'FontSize',6,'rotation',0);
    set(HVenn(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); 
    set(HVenn(2),'FaceColor',tcolors{2},'FaceAlpha',0.75);
%     set(HVenn(3),'FaceColor',tcolors{3},'FaceAlpha',0.75);
    ylims = ylim;
%     text(1,ylims(2)+0.51,sprintf('%s',event_type{ind}),'FontSize',6);
    text(5,ylims(2)-0.51,sprintf('%.0f%%',IVenn(end)),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('conjunctive_cells_venn_diagram_%s.pdf',event_type{ind}),600);
    end
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

%% Brake (before-after vs No-Brake
while 1
    all_resp = [];
    for ind = 1:3
    ntrials = 50;
    event_type = {'Air ON','Air OFF','Light ON'};
%     sic = {[Lb Lbs Ab_On Abs_On];[Ab_Off Abs_Off];ArL_L;[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off]};%;M_On;M_Off};
    allsic = {{[Ab_On];[Ar_On ArL_On Ars_On]; Abs_On};
    {[Ab_Off];[Ar_Off ArL_Off Ars_Off];Abs_Off};
    {[Lb];ArL_L;Lbs}};%;M_On;M_Off};
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
    xticks = xdata; xticklabels = {'Brake-B','No-Brake','Brake-A'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.035 0.01 -0.05 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -5 0]});
    pos = get(gca,'Position');
    text(1,ylims(2),sprintf('%s',event_type{ind}),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('unique_cells_across_conditions_%s.pdf',event_type{ind}),600);
    
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
    AVenn = [p_gFR_C(1,1) p_gFR_C(2,2) p_gFR_C(3,3)]; IVenn = [p_gFR_C(1,2) p_gFR_C(1,3) p_gFR_C(2,3) p_gFR_C_all];
    [HVenn] = venn(AVenn,IVenn,'ErrMinMode','None');
    format_axes(gca);
    axis equal; axis off;
    changePosition(gca,[0.0 -0.0 -0.05 -0.05]); %put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -5 0]});
    pmchar=char(177);
    text(-3.65,2.5,{'Brake-B',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(1,1),pmchar,sem_p_gFR_C(1,1))},'FontSize',6,'rotation',0);
    text(4.5,-2.65,{'No-Brake',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(2,2),pmchar,sem_p_gFR_C(2,2))},'FontSize',6,'rotation',0);
    text(0,4,{sprintf('Brake-A (%0.0f%c%0.0f%%)',m_p_gFR_C(3,3),pmchar,sem_p_gFR_C(3,3))},'FontSize',6,'rotation',0);
    set(HVenn(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); 
    set(HVenn(2),'FaceColor',tcolors{2},'FaceAlpha',0.75);
    set(HVenn(3),'FaceColor',tcolors{3},'FaceAlpha',0.75);
    ylims = ylim;
%     text(1,ylims(2)+0.51,sprintf('%s',event_type{ind}),'FontSize',6);
    text(5,ylims(2)-0.51,sprintf('%.0f%%',IVenn(end)),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('conjunctive_cells_venn_diagram_%s.pdf',event_type{ind}),600);
    end
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

%% Brake vs No-Brake simplified air and light
while 1
    all_resp = [];
    for ind = 2
    ntrials = 50;
    event_type = {'Air','Light'};
    allsic = {{[Ab_On Abs_On Ab_Off Abs_Off];[Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off]};
    {[Lb Lbs];ArL_L}};%;M_On;M_Off};
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
    xticks = xdata; xticklabels = {'Brake','No-Brake','IT-Interval'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.035 0.01 -0.05 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -5 0]});
    pos = get(gca,'Position');
    text(1,ylims(2),sprintf('%s',event_type{ind}),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('unique_cells_across_conditions_%s.pdf',event_type{ind}),600);
    
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
    text(-3.65,2.5,{'Brake',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(1,1),pmchar,sem_p_gFR_C(1,1))},'FontSize',6,'rotation',0);
    text(4.5,-2.65,{'No-Brake',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(2,2),pmchar,sem_p_gFR_C(2,2))},'FontSize',6,'rotation',0);
    set(HVenn(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); set(HVenn(2),'FaceColor',tcolors{2},'FaceAlpha',0.75);
%     text(-6,6,{'Volunteer',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(3,3),pmchar,sem_p_gFR_C(3,3))},'FontSize',6,'rotation',0);
    ylims = ylim;
    text(1,ylims(2)+0.51,sprintf('%s',event_type{ind}),'FontSize',6);
    text(5,ylims(2)-0.51,sprintf('%.0f%%',IVenn),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('conjunctive_cells_venn_diagram_%s.pdf',event_type{ind}),600);
    end
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

%% Brake vs No-Brake OVERLAP HEATMAP
while 1
    all_resp = [];
    ind = 1;
    ntrials = 50;
    event_type = {'OverlapInd'};
    allsic = {{[Ab_On Abs_On Ab_Off Abs_Off];[Lb Lbs];
        [Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off];ArL_L;
        [M_On M_Off]};[Ar_D ArL_D Ars_D];[Ar_T ArL_T Ars_T]};
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
    xticks = xdata; xticklabels = {'Brake','No-Brake','IT-Interval'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.035 0.01 -0.05 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -5 0]});
    pos = get(gca,'Position');
    text(1,ylims(2),sprintf('%s',event_type{ind}),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('unique_cells_across_conditions_%s.pdf',event_type{ind}),600);
    
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

%% Brake vs No-Brake Air on Exc and Inh
while 1
    
    hf = get_figure(6,[10 7 1.5 1]);
    good_FRV = good_FR(:,1:2); good_FRV = good_FR(:,3:4);
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

%% Brake vs No-Brake
while 1
    all_resp = [];
    for ind = 1
    ntrials = 50;
    event_type = {'Air ON','Air OFF','Light ON'};
%     sic = {[Lb Lbs Ab_On Abs_On];[Ab_Off Abs_Off];ArL_L;[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off]};%;M_On;M_Off};
    allsic = {{[Ab_On Abs_On];[Ar_On ArL_On Ars_On]};
    {[Ab_Off Abs_Off];[Ar_Off ArL_Off Ars_Off]};
    {[Lb Lbs];ArL_L}};%;M_On;M_Off};
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
    [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(good_FR,0.5,0.05);

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
    xticks = xdata; xticklabels = {'Brake','No-Brake','IT-Interval'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.035 0.01 -0.05 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -5 0]});
    pos = get(gca,'Position');
    text(1,ylims(2),sprintf('%s',event_type{ind}),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('unique_cells_across_conditions_%s.pdf',event_type{ind}),600);
    
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
    text(-3.65,2.5,{'Brake',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(1,1),pmchar,sem_p_gFR_C(1,1))},'FontSize',6,'rotation',0);
    text(4.5,-2.65,{'No-Brake',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(2,2),pmchar,sem_p_gFR_C(2,2))},'FontSize',6,'rotation',0);
    set(HVenn(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); set(HVenn(2),'FaceColor',tcolors{2},'FaceAlpha',0.75);
%     text(-6,6,{'Volunteer',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(3,3),pmchar,sem_p_gFR_C(3,3))},'FontSize',6,'rotation',0);
    ylims = ylim;
    text(1,ylims(2)+0.51,sprintf('%s',event_type{ind}),'FontSize',6);
    text(5,ylims(2)-0.51,sprintf('%.0f%%',IVenn),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('conjunctive_cells_venn_diagram_%s.pdf',event_type{ind}),600);
    end
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
    tcolors = mData.dcolors(4:6);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',6,'barWidth',0.5,'sigLinesStartYFactor',0.17);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Brake','No-Brake','Motion'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.035 0.01 -0.25 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 0 0]});
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
    resp = all_resp;
   [OI,mOI,semOI,OI_mat,p_vals,h_vals,CI,mCI,semCI,CI_mat] = get_overlap_index(resp,0.5,0.05);
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

%% Brake vs No-Brake and Spatial
while 1
    all_resp = [];
    for ind = 3
    ntrials = 50;
    event_type = {'Air','Light','All'};
    allsic = {{[Ab_On Abs_On Ab_Off Abs_Off];[Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off];[Ar_D ArL_D Ars_D]};
    {[Lb Lbs];ArL_L;[Ar_D ArL_D Ars_D]};
    {[Ab_On Abs_On Ab_Off Abs_Off Lb Lbs];[Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off];[M_On M_Off];[Ar_D ArL_D Ars_D]}};
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
    xticks = xdata; xticklabels = {'Brake','No-Brake','Motion','Dist'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.035 0.01 -0.05 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -5 0]});
    pos = get(gca,'Position');
    text(1,ylims(2),sprintf('%s',event_type{ind}),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('unique_cells_across_conditions_%s_dist.pdf',event_type{ind}),600);
    
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
    
    [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(good_FRV,0.5,0.05);
    
    gFR_C_all = cell_list_op(good_FRV(:,1),good_FRV(:,2),'and',1);    gFR_C_all = cell_list_op(gFR_C_all,good_FRV(:,3),'and',1);
    p_gFR_C_all = mean(100*exec_fun_on_cell_mat(gFR_C_all,'sum')./exec_fun_on_cell_mat(gFR_C_all,'length'));
    [mp_gFR_C_all,semp_gFR_C_all] = findMeanAndStandardError(100*exec_fun_on_cell_mat(gFR_C_all,'sum')./exec_fun_on_cell_mat(gFR_C_all,'length'));
%     AVenn = [p_gFR_C(1,1) p_gFR_C(2,2)]; IVenn = [p_gFR_C(1,2)];
%     [HVenn] = venn(AVenn,IVenn,'ErrMinMode','None');
    AVenn = [p_gFR_C(1,1) p_gFR_C(2,2) p_gFR_C(3,3)]; IVenn = [p_gFR_C(1,2) p_gFR_C(1,3) p_gFR_C(2,3) p_gFR_C_all];
    [HVenn] = venn(AVenn,IVenn,'ErrMinMode','None');
    format_axes(gca);
    axis equal; axis off;
    changePosition(gca,[0.0 -0.0 -0.05 -0.05]); %put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -5 0]});
    pmchar=char(177);
    text(-3.65,2.5,{'Brake',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(1,1),pmchar,sem_p_gFR_C(1,1))},'FontSize',6,'rotation',0);
    text(4.5,-2.65,{'No-Brake',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(2,2),pmchar,sem_p_gFR_C(2,2))},'FontSize',6,'rotation',0);
    text(0,4,{sprintf('Dist (%0.0f%c%0.0f%%)',m_p_gFR_C(3,3),pmchar,sem_p_gFR_C(3,3))},'FontSize',6,'rotation',0);
    set(HVenn(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); 
    set(HVenn(2),'FaceColor',tcolors{2},'FaceAlpha',0.75);
    set(HVenn(3),'FaceColor',tcolors{3},'FaceAlpha',0.75);
    ylims = ylim;
%     text(1,ylims(2)+0.51,sprintf('%s',event_type{ind}),'FontSize',6);
%     text(5,ylims(2)-0.51,sprintf('%.0f%%',IVenn(end)),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('conjunctive_cells_venn_diagram_%s_dist.pdf',event_type{ind}),600);
    end
    %%
    any_cond = cell_list_op(good_FRV,[],'or',1);
    all_cond = cell_list_op(good_FRV,[],'and',1);
    p_any_cond = 100*exec_fun_on_cell_mat(any_cond,'sum')./exec_fun_on_cell_mat(any_cond,'length');
    p_all_cond = 100*exec_fun_on_cell_mat(all_cond,'sum')./exec_fun_on_cell_mat(all_cond,'length');
    [mpac,sempac] = findMeanAndStandardError(p_any_cond);
    [mpall,sempall] = findMeanAndStandardError(p_all_cond);
    
   [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(good_FRV,0.5,0.05);
   mOI = mCI; semOI = semCI;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
%     mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = {'Brake','No-Brake','Motion','Dist'};%,'M-On','M-Off'};
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
%     imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 1.65 1.5]);
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
    changePosition(gca,[0.03 0 -0.06 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',ntrials),600);
    %%
    break;
end

%% Brake vs No-Brake and temporal
while 1
    all_resp = [];
    for ind = 2
    ntrials = 50;
    event_type = {'Air','Light'};
    allsic = {{[Ab_On Abs_On Ab_Off Abs_Off];[Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off];[Ar_T ArL_T Ars_T]};
    {[Lb Lbs];ArL_L;[Ar_D ArL_D Ars_D]}};
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
    xticks = xdata; xticklabels = {'Brake','No-Brake','Time'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.035 0.01 -0.05 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -5 0]});
    pos = get(gca,'Position');
    text(1,ylims(2),sprintf('%s',event_type{ind}),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('unique_cells_across_conditions_%s_time.pdf',event_type{ind}),600);
    
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
    AVenn = [p_gFR_C(1,1) p_gFR_C(2,2) p_gFR_C(3,3)]; IVenn = [p_gFR_C(1,2) p_gFR_C(1,3) p_gFR_C(2,3) p_gFR_C_all];
    [HVenn] = venn(AVenn,IVenn,'ErrMinMode','None');
    format_axes(gca);
    axis equal; axis off;
    changePosition(gca,[0.0 -0.0 -0.05 -0.05]); %put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -5 0]});
    pmchar=char(177);
    text(-3.65,2.5,{'Brake',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(1,1),pmchar,sem_p_gFR_C(1,1))},'FontSize',6,'rotation',0);
    text(4.5,-2.65,{'No-Brake',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(2,2),pmchar,sem_p_gFR_C(2,2))},'FontSize',6,'rotation',0);
    text(0,4,{sprintf('Time (%0.0f%c%0.0f%%)',m_p_gFR_C(3,3),pmchar,sem_p_gFR_C(3,3))},'FontSize',6,'rotation',0);
    set(HVenn(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); 
    set(HVenn(2),'FaceColor',tcolors{2},'FaceAlpha',0.75);
    set(HVenn(3),'FaceColor',tcolors{3},'FaceAlpha',0.75);
    ylims = ylim;
%     text(1,ylims(2)+0.51,sprintf('%s',event_type{ind}),'FontSize',6);
%     text(5,ylims(2)-0.51,sprintf('%.0f%%',IVenn(end)),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('conjunctive_cells_venn_diagram_%s_time.pdf',event_type{ind}),600);
    end
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

%% Vol Motion ON OFF, Spatial, and Temporal
while 1
    all_resp = [];
    for ind = 1
    ntrials = 50;
    event_type = {'MoSpTe'};
    allsic = {{[M_On M_Off];[Ar_D ArL_D Ars_D];[Ar_T ArL_T Ars_T]}};
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
    xticks = xdata; xticklabels = {'Volun','Dist','Time'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.035 0.01 -0.05 0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -5 0]});
    pos = get(gca,'Position');
    text(1,ylims(2),sprintf('%s',event_type{ind}),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('unique_cells_across_conditions_%s_time.pdf',event_type{ind}),600);
    
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
    AVenn = [p_gFR_C(1,1) p_gFR_C(2,2) p_gFR_C(3,3)]; IVenn = [p_gFR_C(1,2) p_gFR_C(1,3) p_gFR_C(2,3) p_gFR_C_all];
    [HVenn] = venn(AVenn,IVenn,'ErrMinMode','None');
    format_axes(gca);
    axis equal; axis off;
    changePosition(gca,[0.0 -0.0 -0.05 -0.05]); %put_axes_labels(gca,{[],[0 0 0]},{'Unique Cells (%)',[0 -5 0]});
    pmchar=char(177);
    text(-3.65,2.5,{'Volun',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(1,1),pmchar,sem_p_gFR_C(1,1))},'FontSize',6,'rotation',0);
    text(4.5,-2.65,{'Dist',sprintf('%0.0f%c%0.0f%%',m_p_gFR_C(2,2),pmchar,sem_p_gFR_C(2,2))},'FontSize',6,'rotation',0);
    text(0,4,{sprintf('Time (%0.0f%c%0.0f%%)',m_p_gFR_C(3,3),pmchar,sem_p_gFR_C(3,3))},'FontSize',6,'rotation',0);
    set(HVenn(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); 
    set(HVenn(2),'FaceColor',tcolors{2},'FaceAlpha',0.75);
    set(HVenn(3),'FaceColor',tcolors{3},'FaceAlpha',0.75);
    ylims = ylim;
%     text(1,ylims(2)+0.51,sprintf('%s',event_type{ind}),'FontSize',6);
%     text(5,ylims(2)-0.51,sprintf('%.0f%%',IVenn(end)),'FontSize',6);
    save_pdf(hf,mData.pdf_folder,sprintf('conjunctive_cells_venn_diagram_%s_time.pdf',event_type{ind}),600);
    end
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
    set(hf,'Position',[7 3 1.75 1]);
    set(H,'linewidth',1);
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel('Eucledian Distance');changePosition(hx,[0 -0.1 0]);
    changePosition(gca,[0.03 0.0 0.05 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);
    %%
    break;
end

%% Speed Figure
while 1
    si = [Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off];
    Rs = o.Rs(:,si);
    for ii = 1:length(ei)
        b1 = ei{ii}.b;
        for jj = 1:10
            alds(ii,jj) = b1.dist(b1.stim_r(jj+10)) - b1.dist(b1.air_puff_r(jj+20));
        end
    end
    ald = round(mean(alds(:)));
    %%
     ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[2 3],'spaceRowsCols',[0.1 0.08],'rightUpShifts',[0.13 0.2],'widthHeightAdjustment',[-125 -175]);
    set(gcf,'color','w'); set(gcf,'Position',[5 5 2 1]);
%     [Y,E] = discretize(1:49,3);
    all_speeds = []; cTxt = {'NB-A','NB-AL','NB-A*','Ar-Off','ArL-Off','Ar*-Off'}; 
    for cn = 1:6
        mean_speed_over_trials = [];
        aThisSpeed = [];
        for an = 1:size(Rs,1)
            thisSpeed = nanmean(Rs{an,cn}.speed);
            mean_speed_over_trials(an,:) = thisSpeed;
        end
        if ismember(cn,[1 2 3])
            axes(ff.h_axes(1,cn));
            delete(gca);
            continue;
        else
            axes(ff.h_axes(2,cn-3));
        end
        hold on;
        cis = Rs{1,cn}.resp.cis;
        xs = Rs{1,cn}.xs-Rs{1,cn}.xs(cis(1,2)); N = length(xs);
        mspeed = mean(mean_speed_over_trials(:,1:N)); semspeed = std(mean_speed_over_trials(:,1:N))/sqrt(5);
        plot(xs,mspeed);
        shadedErrorBar(xs,mspeed,semspeed);
%         changePosition(gca,[0.1 0.15 -0.05 -0.15]);
%         if cn == 1 || cn ==3
%             put_axes_labels(gca,{'',[0 0 0]},{{'Speed (cm/sec)'},[0 00 0]});
%         end
        xbTxt = [2.5 7.5 12.5]-3; ybTxt = 31;
%         text(xbTxt(1),ybTxt+5,cTxt{cn},'FontSize',5);
        ylim([0 25]);
        box off;
        format_axes(gca);
        plot([0 0],[0 30],'m','linewidth',0.25);
%         if cn == 3
%             hx = xlabel('Time (sec)');
%         end
        if ismember(cn,[2 3 5 6])
            set(gca,'YTick',[])
        end
        format_axes(gca);
%         if ismember(cn,[1 2 3])
%             set(gca,'XTick',[])
%         end
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('speeds_345'),600);
   
    %%
    break;
end

%% Overlap Indices ImageSC
while 1
    ntrials = 50;
    sic = {[Ab_On Abs_On]; [Ab_Off Abs_Off];[Ar_On ArL_On Ars_On]; [Ar_Off ArL_Off Ars_Off];[Lb Lbs];ArL_L};%;M_On;M_Off};
%     ind = 3;
    all_resp = [];
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
    %%
%     ntrials = 50;
%     si = [Ab_On Abs_On Ar_On ArL_On Ars_On];
%     props1 = get_props_Rs(o.Rs(:,si),ntrials);
%     all_resp = props1.vals;
    resp_OI = good_FR;
    [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(resp_OI,0.5,0.05);
%     mOI = mCI; semOI = semCI;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = {'B-A-Act','NB-A-Act','B-A-Sup','NB-A-Sup'}; 
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
    %
%     maxYss = [maxYs(9) maxYs(2) maxYs(7) maxYs(8) maxYs(7) maxYs(8) maxYs(7) maxYs(8) maxYs(9) maxYs(2) maxYs(11)];
    for sel_row = 1:sz
        sel_row
        for rr = 1:size(OI_mat,1)
            for cc = 1:size(OI_mat,2)
                if isnan(OI_mat(rr,cc,1))
                    continue;
                end
                if h_vals(rr,cc) == 1
        %             text(cc,rr,'*','Color','r','FontSize',12);
                end
            end
        end
    %     sel_row = 4;
        OIs = squeeze(OI_mat(sel_row,:,:))'; OIs(:,sel_row) = [];
        [within,dvn,xlabels1] = make_within_table({'Cond'},[size(OI_mat,1)-1]);
        dataT = make_between_table({OIs},dvn);
        ra = RMA(dataT,within);
        ra.ranova
        [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
        xdata = 1:size(OI_mat,1); 
        xdata(sel_row) = [];
%         xdata = xdataG;
        h(h==1) = 0;
        hf = get_figure(5,[8 7 1.95 1.5]);
        tcolors = mData.colors(setdiff(1:size(OI_mat,1),sel_row));
        [hbs,maxY] = plot_bars_p_table(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k','ptable',extras.pvalsTable,...
            'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',4,'barWidth',0.45);
        ylims = ylim;
        format_axes(gca);
        set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
        xticks = xdata; xticklabels = txl;
        set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
        changePosition(gca,[0.09 0.01 -0.03 -0.55]); put_axes_labels(gca,{[],[0 0 0]},{sprintf('Overlap Index of %s',txl{sel_row}),[-1.45 0.1 0]});
        ha = gca; ptable = extras.pvalsTable;
        display_p_table_img(ha,hbs,[0 0.23 0 0.37],ptable);ytickangle(0);
        save_pdf(hf,mData.pdf_folder,sprintf('OI_bar_%d.pdf',sel_row),600);
        maxYs(sel_row) = maxY;
    end
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
    set(hf,'Position',[7 3 1.75 1]);
    set(H,'linewidth',1);
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel('Eucledian Distance');changePosition(hx,[0 -0.1 0]);
    changePosition(gca,[0.03 0.0 0.05 -0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);
    %%
    break;
end
