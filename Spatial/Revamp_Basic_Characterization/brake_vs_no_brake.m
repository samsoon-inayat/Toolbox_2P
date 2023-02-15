function brake_vs_no_brake

%% Load Data
while 1
    mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
    ei = evalin('base','ei'); 
    
    selContexts = [1 4 6 2 7 2 7 3 4 5 3 4 5 0 0 3 4 5 3 4 5 3 5 2 3 4 5 7 3 4 5 2 3 4 5 7];
    rasterNames = {'light22T','light22T','light22T','airOnsets22T','airOnsets22T','airOffsets22T',...
                    'airOffsets22T','airOnsets22T','airOnsets22T','airOnsets22T',...
                    'airOffsets22T','airOffsets22T','airOffsets22T','motionOnsets','motionOffsets',...
                    'airD','airD','airD','airIT','airIT','airIT','light22T','light22T','airOffsets22_C','airOffsets22_C','airOffsets22_C','airOffsets22_C','airOffsets22_C'....
                    'airOnsets22_C','airOnsets22_C','airOnsets22_C','airOnsets22P','airOnsets22P','airOnsets22P','airOnsets22P','airOnsets22P'
                    };
    rasterNamesTxt = {'1-L','4-L','6-L','2-AOn','7-AOn','2-AOff','7-AOff',...
        '3-AOn','4-AOn','5-AOn','3-AOff','4-AOff','5-AOff','MOn','MOff',...
        '3-D','4-D','5-D','3-T','4-T','5-T','3-Lc','5-Lc','2-AOffc','3-AOffc','4-AOffc','5-AOffc','7-AOffc',...
        '3-AOnc','4-AOnc','5-AOnc','2-AOnp','3-AOnp','4-AOnp','5-AOnp','7-AOnp'};
    tic
    o = get_data(ei,selContexts,rasterNames);
    toc
    
    Lb = 1; ArL_L = 2; Lbs = 3; Ab_On = 4; Abs_On = 5; Ab_Off = 6; Abs_Off = 7; 
    Ar_On = 8; ArL_On = 9; Ars_On = 10; Ar_Off = 11; ArL_Off = 12; Ars_Off = 13;
    M_On = 14; M_Off = 15; Ar_D = 16; ArL_D = 17; Ars_D = 18; Ar_T = 19; ArL_T = 20; Ars_T = 21;
    Ar_L = 22; Ars_L = 23; Ab_Offc = 24; Ar_Offc = 25; ArL_Offc = 26; Ars_Offc = 27; Abs_Offc = 28;
    Ar_Onc = 29; ArL_Onc = 30; Ars_Onc = 31;
    Ab_Onp = 32; Ar_Onp = 33; ArL_Onp = 34; Ars_Onp = 35; Abs_Onp = 36;
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
    cellNums = {[436 72 352 92 168 387 436];
        [168 98 70 59 328 558];
        [191 140 39 17 8 66 193 140];
        [575 580 581 373 532];
        [567 142 66 96 54 429 235 164];
        [353 270 42 215 353 439];
        [552 567 165 103]};
    an = 1; cn = 3;
    ntrials = 50;
    asi = [Lb ArL_L Ab_On Abs_Off Ar_On Ar_Off ];
    si = asi(cn);
%     si = [Lb];
    Rs = o.Rs(:,si);
    props1 = get_props_Rs(Rs,ntrials);
% %     plotRasters_simplest(Rs{an},find(props1.vals{an}))
%     plotRasters_simplest(Rs{an},find(props1.good_FR_and_untuned{an}))
%     break;
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

%% Show sample rasters
while 1
    rtype = {'B_L','NB_L','B_AOn','B_AOff','NB_AOn','NB_AOff','NB-A'};
    cellNums = {[436 72 352 92 168 387 436];
        [168 98 70 59 328 558];
        [191 140 39 17 8 66 193 140];
        [575 580 581 373 532];
        [567 142 66 96 54 429 235 164];
        [353 270 42 215 353 439];
        [552 567 165 103]};
    an = 1; cn = 1;
    ntrials = 50;
    asi = [Lb ArL_L Ab_On Abs_Off Ar_On Ar_Off];
%     asi = [Lb ArL_L Ab_On Abs_Off Ar_On Ar_Off Ar_D Ar_T];
    si = asi(cn);
%     si = [Lb];
    Rs = o.Rs(:,si);
    props1 = get_props_Rs(Rs,ntrials);
% %     plotRasters_simplest(Rs{an},find(props1.vals{an}))
%     plotRasters_simplest(Rs{an},find(props1.good_FR_and_untuned{an}))
%     break;
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
    an = 1; cn = 3;
    ntrials = 50;
    asi = [Lb ArL_L Ab_On Abs_Off Ar_On Ar_Off];
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
    si = [Ab_On Abs_On Ar_On ArL_On Ars_On ];
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    good_FR = props1.vals;
    all_zMIs = props1.zMI;
    zMIs = [];
    for rr = 1:size(all_zMIs,1)
        for cc = 1:size(all_zMIs,2)
            resp = good_FR{rr,cc};
            tzmis = all_zMIs{rr,cc};
            zMIs(rr,cc) = nanmean(tzmis(resp));
        end
    end
%     azMIs = zMIs
    [within,dvn,xlabels] = make_within_table({'Cond'},[2]);
    dataT = make_between_table({zMIs},dvn);
    ra = RMA(dataT,within);
    ra.ranova
%     ra.mauchly
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Event','hsd'},[1 1 1]);
    xdata = make_xdata([2 2 2 1 3 3 2],[1 1.5]);
    xdata = make_xdata([2],[1 1.5]);
    ptab = 0;
    if ptab h(h==1) = 0; end
    hf = get_figure(5,[8 7 3.5 2]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    

    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'BaseValue',0.01,'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'yspacing',0.5);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rasterNamesTxt(si);
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 round(maxY)]); xtickangle(45);
    if ptab
    changePosition(gca,[0.08 0.01 0.0 0]); put_axes_labels(gca,{[],[0 0 0]},{'Avg. zMI',[0 0 0]});
    ha = gca; ptable = extras.pvalsTable;
    display_p_table(ha,hbs,[0 0.0 0 0],ptable);
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

%% compare percent responsive cells pooled light ON brake vs no-brake exc and inh and Control Light (from Condition 3 and 5 where there is no light stimulus)
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    allsic = {{[Lb Lbs];[ArL_L];[Ar_L];[Ars_L]};
    };
%     allsic = {{[ArL_L];[Ar_L];[Ars_L]};
%     };
    allsic = {{[ArL_L];[Ar_L];[Ars_L]};
        };
%     allsic = {{[Lb];[Lbs];[ArL_L]};};
    
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
    [within,dvn,xlabels] = make_within_table({'Type','Cond'},[2,length(sic)]);
    dataT = make_between_table({per_active},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
  
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


%% compare percent responsive cells pooled light ON brake vs no-brake exc and inh
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
  
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Type_by_Cond','hsd'},[1.5 1 1]);
	xdata = make_xdata([2 2],[1 1.5]);
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
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 5 10]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.5 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
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
  
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Type_by_Cond','hsd'},[1.5 1 1]);
    h(2)= 0;
	xdata = make_xdata([2 2],[1 1.5]);
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
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 35 70]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.5 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{'Trials (%)'},[0 0 0]});
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

%% compare percent responsive cells pooled air on or off brake vs no-brake exc and inh
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    event_type = {'Air ON','Air OFF','Light ON'};
%     sic = {[Lb Lbs Ab_On Abs_On];[Ab_Off Abs_Off];ArL_L;[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off]};%;M_On;M_Off};
    allsic = {{[Ab_On Abs_On];[Ar_On ArL_On Ars_On]}};
%     allsic = {{[Ab_Off Abs_Off];[Ar_Off ArL_Off Ars_Off]}};
%     allsic = {{[Ab_Offc Abs_Offc];[Ar_Offc ArL_Offc Ars_Offc]}};
%     ind = 3;
%     sic = allsic{1:6};
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
            gFR = props1.good_FR;
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
    print_for_manuscript(ra)
    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
    
  
    %%
    allrp = all_resp;
    per_active1 = 100*exec_fun_on_cell_mat(allrp,'sum')./exec_fun_on_cell_mat(allrp,'length');
    [within,dvn,xlabels] = make_within_table({'Cond'},[2]);
    dataT = make_between_table({per_active1},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    
   
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
    allsic = {{[Ab_Off Abs_Off];[Ar_Off ArL_Off Ars_Off]};};
    allsic = {{[Ab_Offc Abs_Offc];[Ar_Offc ArL_Offc Ars_Offc]}};
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
  
    
   
    %%
    break;
end


%% compare percent zMI pooled air on or off brake vs no-brake
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    event_type = {'Air ON','Air OFF','Light ON'};
%     sic = {[Lb Lbs Ab_On Abs_On];[Ab_Off Abs_Off];ArL_L;[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off]};%;M_On;M_Off};
    allsic = {{[Ab_On Abs_On];[Ar_On ArL_On Ars_On]};
    };
%     allsic = {{[Ab_Off Abs_Off];[Ar_Off ArL_Off Ars_Off]};};
%     allsic = {{[Ab_Offc Abs_Offc];[Ar_Offc ArL_Offc Ars_Offc]}};
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
            gFR = props1.vals; rf = props1.HaFD;
            all_gFR(:,ii) = mean(exec_fun_on_cell_mat(rf,'mean',gFR),2);
            gFR = props1.exc; rf = props1.HaFD;
            all_exc(:,ii) = mean(exec_fun_on_cell_mat(rf,'mean',gFR),2);
            gFR = props1.inh;rf = props1.HaFD;
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


%% compare percent responsive cells pooled air on or off brake vs no-brake exc and inh (and controls) big ANOVA test
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    event_type = {'Air ON','Air OFF','Light ON'};
    sic = {[Ab_On Abs_On];[Ab_Off Abs_Off];[Ab_Offc Abs_Offc];[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off];[Ar_Offc ArL_Offc Ars_Offc]};
    clear all_gFR all_exc all_inh all_gV 
    all_exc_inh = [];
    for ii = 1:length(sic)
        sit = sic{ii};
        tRs = o.Rs(:,sit);
        props1 = get_props_Rs(tRs,ntrials);
        gFR = props1.vals;
        gFR = cell_list_op(gFR,[],'or',1);
        all_gV(:,ii) = gFR;
        gFR = props1.good_FR;
        gFR = cell_list_op(gFR,[],'or',1);
        all_gFR(:,ii) = gFR;
        gFR = props1.exc;
        gFR = cell_list_op(gFR,[],'or',1);
        all_exc(:,ii) = gFR;
        gFR = props1.inh;
        gFR = cell_list_op(gFR,[],'or',1);
        all_inh(:,ii) = gFR;
        all_exc_inh = [all_exc_inh all_exc(:,ii) all_inh(:,ii)];
    end
    per_active = find_percent(all_exc_inh);
    [within,dvn,xlabels] = make_within_table({'Cond','ST','CT'},[2,3,2]);
    dataT = make_between_table({per_active},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    print_for_manuscript(ra)
    %%
    break;
end


%% Three-Way RM-ANOVA Most Recent compare all variables pooled air on or off brake vs no-brake exc and inh (and controls) big ANOVA test
while 1
    ntrials = 50; 
%     si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
%     event_type = {'2-B-AOn','7-B-AOn','2-B-AOff','7-B-AOff','2-B-Arb','7-B-Arb','3-NB-AOn','4-NB-AOn','5-NB-AOn','3-NB-AOff',...
%         '4-NB-AOff','5-NB-AOff','3-NB-Arb','4-NB-Arb','5-NB-Arb','1-B-L','6-B-L','3-NB-A','5-NB-A','4-NB-AL'};
%     
    event_type = {'B-AOn-Exc','B-AOn-Inh','B-AOff-Exc','B-AOff-Inh','B-Arb-Exc','B-Arb-Inh','NB-AOn-Exc','NB-AOn-Inh','NB-AOff-Exc','NB-AOff-Inh','NB-Arb-Exc','NB-Arb-Inh'};
    sic = {[Ab_On Abs_On];[Ab_Off Abs_Off];[Ab_Offc Abs_Offc];[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off];[Ar_Offc ArL_Offc Ars_Offc]};
    
%     sic = {[Ab_On Abs_On Ab_Off Abs_Off Ab_Offc Abs_Offc];[Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off Ar_Offc ArL_Offc Ars_Offc]};
    
%     event_type = {'B-AOn','B-AOff','B-Arb','NB-AOn','NB-AOff','NB-Arb'};
    
%     sic = {[Ab_On];[Abs_On];[Ab_Off];[Abs_Off];[Ab_Offc];[Abs_Offc];[Ar_On];[ArL_On];[Ars_On];[Ar_Off];[ArL_Off];[Ars_Off];[Ar_Offc];[ArL_Offc];[Ars_Offc];[Lb];[Lbs];[Ar_L];[Ars_L];[Ar_L]}; % for heat map after light sitmulus
    
%     event_type = {'2-B-AOn','7-B-AOn','2-B-AOff','7-B-AOff','3-NB-AOn','4-NB-AOn','5-NB-AOn','3-NB-AOff',...
%         '4-NB-AOff','5-NB-AOff','1-B-L','6-B-L','4-NB-AL','3-NB-A','5-NB-A'};
%     sic = {[Ab_On];[Abs_On];[Ab_Off];[Abs_Off];[Ar_On];[ArL_On];[Ars_On];[Ar_Off];[ArL_Off];[Ars_Off];[Lb];[Lbs];[ArL_L];[Ar_L];[Ars_L]}; % for heat map after light sitmulus
%     sic = {[Ab_On Abs_On];[Ar_On ArL_On Ars_On];[Ab_Off Abs_Off];[Ar_Off ArL_Off Ars_Off];[Ab_Offc Abs_Offc];[Ar_Offc ArL_Offc Ars_Offc]};
    clear all_gFR all_exc all_inh all_gV all_exc_inh
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
     %
%     avar = find_percent(all_exc_inh);
    avar = (all_exc_inh);
    [within,dvn,xlabels,withinD] = make_within_table({'Cond','ET','CT'},[2,3,2]); withinD3 = withinD;
    dataT = make_between_table({avar},dvn);
    ra = RMA(dataT,within,{0.05,{'lsd','hsd','bonferroni'}});
    ra.ranova
    print_for_manuscript(ra)
    %%
    avar1 = (all_gV);
    [within,dvn,xlabels,withinD] = make_within_table({'Cond','ET'},[2,3]); withinD3 = withinD;
    dataT = make_between_table({avar1},dvn);
    ra_cond_et = RMA(dataT,within,{0.05,{'lsd','hsd','bonferroni'}});
    ra_cond_et.ranova
    print_for_manuscript(ra_cond_et)
    %%
%     r_avar = [];
%     for ii = 1:5
%         r_avar = [r_avar;[ii*ones(size(avar,2)/2,1) [1;1;1;2;2;2] repmat([1;2;3],2,1) all_exc(ii,:)' all_inh(ii,:)']];
%     end
%     r_avarT = array2table(r_avar); r_avarT.Properties.VariableNames = {'Animal','BrakeState','EventType','Exc','Inh'};
%     writetable(r_avarT,'dataT.csv');
    %%
    [within,dvn,xlabels,withinD] = make_within_table({'ET','CT'},[3,2]);
    dataT_1 = make_between_table({avar(:,find(withinD3(:,1)==1))},dvn);
    ra_cond1 = RMA(dataT_1,within,{0.025,{'hsd','bonferroni'}});
    ra_cond1.ranova
    print_for_manuscript(ra_cond1);
    
    dataT_2 = make_between_table({avar(:,find(withinD3(:,1)==2))},dvn);
    ra_cond2 = RMA(dataT_2,within,{0.025,{'hsd','bonferroni'}});
    ra_cond2.ranova
    print_for_manuscript(ra_cond2);
    %%
    alpha = 0.025/3;
    [within,dvn,xlabels,withinD] = make_within_table({'CT'},[2]);
    dataT_1 = make_between_table({avar(:,find(withinD3(:,1)==1 & withinD3(:,2) == 1))},dvn);
    ra_cond1_e1 = RMA(dataT_1,within,{alpha,{'hsd','bonferroni'}});
    ra_cond1_e1.ranova
    print_for_manuscript(ra_cond1_e1);
    
    dataT_1 = make_between_table({avar(:,find(withinD3(:,1)==1 & withinD3(:,2) == 2))},dvn);
    ra_cond1_e2 = RMA(dataT_1,within,{alpha,{'hsd','bonferroni'}});
    ra_cond1_e2.ranova
    print_for_manuscript(ra_cond1_e2);
    
    dataT_1 = make_between_table({avar(:,find(withinD3(:,1)==1 & withinD3(:,2) == 3))},dvn);
    ra_cond1_e3 = RMA(dataT_1,within,{alpha,{'hsd','bonferroni'}});
    ra_cond1_e3.ranova
    print_for_manuscript(ra_cond1_e3);
    
    
    dataT_1 = make_between_table({avar(:,find(withinD3(:,1)==2 & withinD3(:,2) == 1))},dvn);
    ra_cond2_e1 = RMA(dataT_1,within,{alpha,{'hsd','bonferroni'}});
    ra_cond2_e1.ranova
    print_for_manuscript(ra_cond2_e1);
    
    dataT_1 = make_between_table({avar(:,find(withinD3(:,1)==2 & withinD3(:,2) == 2))},dvn);
    ra_cond2_e2 = RMA(dataT_1,within,{alpha,{'hsd','bonferroni'}});
    ra_cond2_e2.ranova
    print_for_manuscript(ra_cond2_e2);
    
    dataT_1 = make_between_table({avar(:,find(withinD3(:,1)==2 & withinD3(:,2) == 3))},dvn);
    ra_cond2_e3 = RMA(dataT_1,within,{alpha,{'hsd','bonferroni'}});
    ra_cond2_e3.ranova
    print_for_manuscript(ra_cond2_e3);
    
    %%
    alpha = 0.025/4;
    [within,dvn,xlabels,withinD] = make_within_table({'ET'},[3]);
    dataT_1 = make_between_table({avar(:,[1 3 5])},dvn);
    ra_cond1_exc = RMA(dataT_1,within,{alpha,{'hsd','bonferroni'}});
    ra_cond1_exc.ranova
    print_for_manuscript(ra_cond1_exc);
    
    dataT_1 = make_between_table({avar(:,[2 4 6])},dvn);
    ra_cond1_inh = RMA(dataT_1,within,{alpha,{'hsd','bonferroni'}});
    ra_cond1_inh.ranova
    print_for_manuscript(ra_cond1_inh);
    
    dataT_1 = make_between_table({avar(:,[1 3 5]+6)},dvn);
    ra_cond2_exc = RMA(dataT_1,within,{alpha,{'hsd','bonferroni'}});
    ra_cond2_exc.ranova
    print_for_manuscript(ra_cond2_exc);
    
    dataT_1 = make_between_table({avar(:,[2 4 6]+6)},dvn);
    ra_cond2_inh = RMA(dataT_1,within,{alpha,{'hsd','bonferroni'}});
    ra_cond2_inh.ranova
    print_for_manuscript(ra_cond2_inh);
    
    
    
    %%
   
    %
    [within,dvn,xlabels,withinD] = make_within_table({'ET','Cond'},[3,2]);
    dataT_1 = make_between_table({avar(:,[1 7 3 9 5 11])},dvn);
    ra_exc = RMA(dataT_1,within,{0.025,{'hsd','bonferroni'}});
    ra_exc.ranova
    print_for_manuscript(ra_exc);
    
    dataT_2 = make_between_table({avar(:,[1 7 3 9 5 11]+1)},dvn);
    ra_inh = RMA(dataT_2,within,{0.025,{'hsd','bonferroni'}});
    ra_inh.ranova
    print_for_manuscript(ra_inh);
    %%
    
    [within,dvn,xlabels,withinD] = make_within_table({'Cond'},[2]);
    dataT_1 = make_between_table({avar(:,[1 7])},dvn);
    ra_exc_e1 = RMA(dataT_1,within,{0.05,{'hsd','bonferroni'}});
    ra_exc_e1.ranova
    print_for_manuscript(ra_exc_e1);
    
    dataT_1 = make_between_table({avar(:,[3 9])},dvn);
    ra_exc_e2 = RMA(dataT_1,within,{0.05,{'hsd','bonferroni'}});
    ra_exc_e2.ranova
    print_for_manuscript(ra_exc_e2);
    
     dataT_1 = make_between_table({avar(:,[5 11])},dvn);
    ra_exc_e3 = RMA(dataT_1,within,{0.05,{'hsd','bonferroni'}});
    ra_exc_e3.ranova
    print_for_manuscript(ra_exc_e3);
    
    dataT_1 = make_between_table({avar(:,[1 7]+1)},dvn);
    ra_inh_e1 = RMA(dataT_1,within,{0.05,{'hsd','bonferroni'}});
    ra_inh_e1.ranova
    print_for_manuscript(ra_inh_e1);
    
    dataT_1 = make_between_table({avar(:,[3 9]+1)},dvn);
    ra_inh_e2 = RMA(dataT_1,within,{0.05,{'hsd','bonferroni'}});
    ra_inh_e2.ranova
    print_for_manuscript(ra_inh_e2);
    
     dataT_1 = make_between_table({avar(:,[5 11]+1)},dvn);
    ra_inh_e3 = RMA(dataT_1,within,{0.05,{'hsd','bonferroni'}});
    ra_inh_e3.ranova
    print_for_manuscript(ra_inh_e3);
    
    %%
    [within,dvn,xlabels,withinD] = make_within_table({'CT','Cond'},[2,2]);
    dataT_1 = make_between_table({avar(:,[1 7 2 8])},dvn);
    ra_e1 = RMA(dataT_1,within,{0.05/3,{'hsd','bonferroni'}});
    ra_e1.ranova
    print_for_manuscript(ra_e1);
    
    dataT_1 = make_between_table({avar(:,[3 9 4 10])},dvn);
    ra_e2 = RMA(dataT_1,within,{0.05/3,{'hsd','bonferroni'}});
    ra_e2.ranova
    print_for_manuscript(ra_e2);
    
    dataT_1 = make_between_table({avar(:,[5 11 6 12])},dvn);
    ra_e3 = RMA(dataT_1,within,{0.05/3,{'hsd','bonferroni'}});
    ra_e3.ranova
    print_for_manuscript(ra_e3);
    
    %%
    [within,dvn,xlabels,withinD] = make_within_table({'Cond'},[2]);
    tempvar = avar(:,[1 7 2 8]);
    dataT_1 = make_between_table({tempvar(:,1:2)},dvn);
    ra_e1_exc = RMA(dataT_1,within);
    ra_e1_exc.ranova
    print_for_manuscript(ra_e1_exc);
    
    dataT_1 = make_between_table({tempvar(:,3:4)},dvn);
    ra_e1_inh = RMA(dataT_1,within);
    ra_e1_inh.ranova
    print_for_manuscript(ra_e1_inh);
    
    tempvar = avar(:,[3 9 4 10]);
    dataT_1 = make_between_table({tempvar(:,1:2)},dvn);
    ra_e2_exc = RMA(dataT_1,within);
    ra_e2_exc.ranova
    print_for_manuscript(ra_e2_exc);
    
    dataT_1 = make_between_table({tempvar(:,3:4)},dvn);
    ra_e2_inh = RMA(dataT_1,within);
    ra_e2_inh.ranova
    print_for_manuscript(ra_e2_inh);
    
    tempvar = avar(:,[5 11 6 12]);
    dataT_1 = make_between_table({tempvar(:,1:2)},dvn);
    ra_e3_exc = RMA(dataT_1,within);
    ra_e3_exc.ranova
    print_for_manuscript(ra_e3_exc);
    
    dataT_1 = make_between_table({tempvar(:,3:4)},dvn);
    ra_e3_inh = RMA(dataT_1,within);
    ra_e3_inh.ranova
    print_for_manuscript(ra_e3_inh);
    
    
   
    %% for normality test
    for ii = 1:size(avar,2)
        [nth(ii,1),ntp(ii,1)] = swtest(avar(:,ii));
    end
    %% for using RMAOV33 function downloaded from the internet
    % I found the same result as my function RMA
    clear avar_33
    ind = 1;
    for cc = 1:size(avar,2)
        for rr = 1:size(avar,1)
            avar_33(ind,:) = [avar(rr,cc) withinD(cc,:) rr];
            ind = ind + 1;
        end
    end
    RMAOV33(avar_33)
    %%
     inds = [1 4];
    good_FRV = all_exc(:,inds);
     [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(good_FRV,0.5,0.05);
    intcells = squeeze(all_CI_mat(1,2,:));    bcells = squeeze(all_CI_mat(1,1,:)); nbcells = squeeze(all_CI_mat(2,2,:));
    bcellsu = squeeze(uni(1,2,:)); nbcellsu = squeeze(uni(2,1,:));
    
    aVar = [bcellsu intcells nbcellsu];
    
    good_FRV = all_inh(:,inds);
    [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(good_FRV,0.5,0.05);
    intcells = squeeze(all_CI_mat(1,2,:));    bcells = squeeze(all_CI_mat(1,1,:)); nbcells = squeeze(all_CI_mat(2,2,:));
    bcellsu = squeeze(uni(1,2,:)); nbcellsu = squeeze(uni(2,1,:));
    
    aVar = [aVar bcellsu intcells nbcellsu];
    
    [within,dvn,xlabels] = make_within_table({'CT','PT'},[2,3]);
    dataT = make_between_table({aVar},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    print_for_manuscript(ra)
    %%
    break;
end

%% Conjunction Complementation of Exc and Inh cells
while 1
    inds_all = {[1 4],[2 5],[3 6]};
    aVar = [];
    for ii = 1:length(inds_all)
        [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(all_exc(:,inds_all{ii}),0.5,0.05);
        intcells = squeeze(all_CI_mat(1,2,:));    bcells = squeeze(all_CI_mat(1,1,:)); nbcells = squeeze(all_CI_mat(2,2,:));
        bcellsu = squeeze(uni(1,2,:)); nbcellsu = squeeze(uni(2,1,:));
        aVar = [aVar bcellsu intcells nbcellsu];
    end
    for ii = 1:length(inds_all)
        [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(all_inh(:,inds_all{ii}),0.5,0.05);
        intcells = squeeze(all_CI_mat(1,2,:));    bcells = squeeze(all_CI_mat(1,1,:)); nbcells = squeeze(all_CI_mat(2,2,:));
        bcellsu = squeeze(uni(1,2,:)); nbcellsu = squeeze(uni(2,1,:));
        aVar = [aVar bcellsu intcells nbcellsu];
    end
    
    [within,dvn,xlabels,withinD] = make_within_table({'CT','ET','PT'},[2,3,3]); withinD3 = withinD;
    dataT = make_between_table({aVar},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd'}});
    ra.ranova
    print_for_manuscript(ra)
    %%
    [within,dvn,xlabels,withinD] = make_within_table({'ET','PT'},[3,3]);
    dataT = make_between_table({aVar(:,find(withinD3(:,1)==1))},dvn);
    ra_exc = RMA(dataT,within,{0.025,{'hsd'}}); ra_exc.ranova
    print_for_manuscript(ra_exc);
    
    dataT = make_between_table({aVar(:,find(withinD3(:,1)==2))},dvn);
    ra_inh = RMA(dataT,within,{0.025,{'hsd'}}); ra_inh.ranova
    print_for_manuscript(ra_inh);
    
    
    %% no simple two way interaction significant
    alpha = 0.05/3;
    [within,dvn,xlabels,withinD] = make_within_table({'CT','PT'},[2,3]);
    dataT = make_between_table({aVar(:,find(withinD3(:,2)==1))},dvn);
    ra_e1 = RMA(dataT,within,{alpha,{'hsd'}}); ra_e1.ranova
    print_for_manuscript(ra_e1);
    
    dataT = make_between_table({aVar(:,find(withinD3(:,2)==2))},dvn);
    ra_e2 = RMA(dataT,within,{alpha,{'hsd'}}); ra_e2.ranova
    print_for_manuscript(ra_e2);
    
    dataT = make_between_table({aVar(:,find(withinD3(:,2)==3))},dvn);
    ra_e3 = RMA(dataT,within,{alpha,{'hsd'}}); ra_e3.ranova
    print_for_manuscript(ra_e3);
    
    %%
    alpha = 0.05/3;
    [within,dvn,xlabels,withinD] = make_within_table({'CT','ET'},[2,3]);
    dataT = make_between_table({aVar(:,find(withinD3(:,3)==1))},dvn);
    ra_p1 = RMA(dataT,within,{alpha,{'hsd'}}); ra_p1.ranova
    print_for_manuscript(ra_p1); % simple two-way interaction significant
    
    dataT = make_between_table({aVar(:,find(withinD3(:,3)==2))},dvn);
    ra_p2 = RMA(dataT,within,{alpha,{'hsd'}}); ra_p2.ranova
    print_for_manuscript(ra_p2);
    
    dataT = make_between_table({aVar(:,find(withinD3(:,3)==3))},dvn);
    ra_p3 = RMA(dataT,within,{alpha,{'hsd'}}); ra_p3.ranova
    print_for_manuscript(ra_p3);
    
    %%
    [within,dvn,xlabels,withinD] = make_within_table({'ET'},[3]);
    dataT = make_between_table({aVar(:,find(withinD3(:,3)==1&withinD3(:,1)==1))},dvn);
    ra_p1_exc = RMA(dataT,within,{alpha/2,{'hsd','bonferroni'}}); ra_p1_exc.ranova
    print_for_manuscript(ra_p1_exc); 
    
    [within,dvn,xlabels,withinD] = make_within_table({'ET'},[3]);
    dataT = make_between_table({aVar(:,find(withinD3(:,3)==1&withinD3(:,1)==2))},dvn);
    ra_p1_inh = RMA(dataT,within,{alpha/2,{'hsd','bonferroni'}}); ra_p1_inh.ranova
    print_for_manuscript(ra_p1_inh); 
    
    dataT = make_between_table({aVar(:,find(withinD3(:,3)==2&withinD3(:,1)==1))},dvn);
    ra_p2_exc = RMA(dataT,within,{alpha/2,{'hsd','bonferroni'}}); ra_p2_exc.ranova
    print_for_manuscript(ra_p2_exc); 
    
    [within,dvn,xlabels,withinD] = make_within_table({'ET'},[3]);
    dataT = make_between_table({aVar(:,find(withinD3(:,3)==2&withinD3(:,1)==2))},dvn);
    ra_p2_inh = RMA(dataT,within,{alpha/2,{'hsd','bonferroni'}}); ra_p2_inh.ranova
    print_for_manuscript(ra_p2_inh); 
    
    
     dataT = make_between_table({aVar(:,find(withinD3(:,3)==3&withinD3(:,1)==1))},dvn);
    ra_p3_exc = RMA(dataT,within,{alpha/2,{'hsd','bonferroni'}}); ra_p3_exc.ranova
    print_for_manuscript(ra_p3_exc); 
    
    [within,dvn,xlabels,withinD] = make_within_table({'ET'},[3]);
    dataT = make_between_table({aVar(:,find(withinD3(:,3)==3&withinD3(:,1)==2))},dvn);
    ra_p3_inh = RMA(dataT,within,{alpha/2,{'hsd','bonferroni'}}); ra_p3_inh.ranova
    print_for_manuscript(ra_p3_inh); 
    
    
    %%
    break;
end

%% Conjunction Complementation of both Exc and Inh pooled cells
while 1
    inds_all = {[1 4],[2 5],[3 6]};
    aVar = [];
    for ii = 1:length(inds_all)
        [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(all_gV(:,inds_all{ii}),0.5,0.05);
        intcells = squeeze(all_CI_mat(1,2,:));    bcells = squeeze(all_CI_mat(1,1,:)); nbcells = squeeze(all_CI_mat(2,2,:));
        bcellsu = squeeze(uni(1,2,:)); nbcellsu = squeeze(uni(2,1,:));
        aVar = [aVar bcellsu intcells nbcellsu];
    end
    
    [within,dvn,xlabels,withinD] = make_within_table({'ET','PT'},[3,3]); withinD3 = withinD;
    dataT = make_between_table({aVar},dvn);
    ra_et_pt = RMA(dataT,within,{0.05,{'hsd','bonferroni'}});
    ra_et_pt.ranova
    print_for_manuscript(ra_et_pt)
    %%
    break;
end



%% compare percent zMI pooled air on or off brake vs no-brake
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    event_type = {'Air ON','Air OFF','Light ON'};
%     sic = {[Lb Lbs Ab_On Abs_On];[Ab_Off Abs_Off];ArL_L;[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off]};%;M_On;M_Off};
    allsic = {{[Ab_On Abs_On];[Ar_On ArL_On Ars_On]};
    };
%     allsic = {{[Ab_Off Abs_Off];[Ar_Off ArL_Off Ars_Off]};};
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
            gFR = props1.vals; rf = props1.zMI_MC;
            all_gFR(:,ii) = mean(exec_fun_on_cell_mat(rf,'nanmean',gFR),2);
            gFR = props1.exc; rf = props1.zMI_MC;
            all_exc(:,ii) = mean(exec_fun_on_cell_mat(rf,'nanmean',gFR),2);
            gFR = props1.inh;rf = props1.zMI_MC;
            all_inh(:,ii) = mean(exec_fun_on_cell_mat(rf,'nanmean',gFR),2);
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
        'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0.5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-Exc','NB-Exc','B-Inh','NB-Inh'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.5 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{'Response','Fidelity (% trials)'},[0 0 0]});
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
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_allf3.pdf',ntrials),600);
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
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'Type_by_Cond','hsd'},[1.5 1 1]);
    h([2 5]) = 0;
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
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 60]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-A','NB-A','B-L','NB-L'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 15 30]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.5 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all1.pdf',ntrials),600);
    
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
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 60]); format_axes(gca);
xticks = xdata; xticklabels = {'B-A','NB-A','B-L','NB-L'};    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 15 30]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.5 -0.04]); %put_axes_labels(gca,{[],[0 0 0]},{{'Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all2.pdf',ntrials),600);
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
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 60]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-A','NB-A','B-L','NB-L'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 15 30]); xtickangle(45)
    changePosition(gca,[0.05 0.01 -0.5 -0.04]); %put_axes_labels(gca,{[],[0 0 0]},{{'Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all3.pdf',ntrials),600);
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

%% compare percent responsive cells pooled air on vs off with brake no brake
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
%     sic = {[Lb Lbs Ab_On Abs_On];[Ab_Off Abs_Off];ArL_L;[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off]};%;M_On;M_Off};
    allsic = {{[Ab_On Abs_On];[Ab_Off Abs_Off]};
    {[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off]};
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
    h(h==1) = 0;
	xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 1 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',10,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 39]); format_axes(gca);
    xticks = xdata; xticklabels = {'AOn','AOff','AOn','AOff'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 15 30]); xtickangle(45)
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
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 39]); format_axes(gca);
    xticks = xdata; xticklabels = {'AOn','AOff','AOn','AOff'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 15 30]); xtickangle(45)
    changePosition(gca,[0.09 0.01 -0.4 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
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
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 39]); format_axes(gca);
    xticks = xdata; xticklabels = {'AOn','AOff','AOn','AOff'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 15 30]); xtickangle(45)
    changePosition(gca,[0.09 0.01 -0.4 -0.04]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
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
    good_FR_grand = [all_exc(:,1),all_inh(:,1),all_exc(:,2),all_inh(:,2),all_exc(:,5),all_inh(:,5),all_exc(:,3),all_inh(:,3),all_exc(:,4),all_inh(:,4),all_exc(:,6),all_inh(:,6)];
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
    pavar = find_percent(good_FR_grand);
    [within,dvn,xlabels] = make_within_table({'Cond','ET','CT'},[2,3,2]);
    dataT = make_between_table({pavar},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    print_for_manuscript(ra)
    %%
    avar = [all_gV(:,[1 3]) all_gV(:,[2 4]) all_gV(:,[5 6])];
    pavar = find_percent(avar);
    [within,dvn,xlabels] = make_within_table({'Type','Cond'},[3,2]);
    dataT = make_between_table({pavar},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    print_for_manuscript(ra)
    
    %%
    avar = [cell_list_op(all_gV(:,[1 2 5]),[],'or',1) cell_list_op(all_gV(:,[3 4 6]),[],'or',1)];
    pavar = find_percent(avar);
    [within,dvn,xlabels] = make_within_table({'Type','Cond'},[3,2]);
    dataT = make_between_table({pavar},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    print_for_manuscript(ra)
    %% combined population type and exc and inh cell type RM anova test
    inds = [1 3];
    good_FRV = all_exc(:,inds);
     [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(good_FRV,0.5,0.05);
    intcells = squeeze(all_CI_mat(1,2,:));    bcells = squeeze(all_CI_mat(1,1,:)); nbcells = squeeze(all_CI_mat(2,2,:));
    bcellsu = squeeze(uni(1,2,:)); nbcellsu = squeeze(uni(2,1,:));
    
    aVar = [bcellsu intcells nbcellsu];
    
    good_FRV = all_inh(:,inds);
    [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(good_FRV,0.5,0.05);
    intcells = squeeze(all_CI_mat(1,2,:));    bcells = squeeze(all_CI_mat(1,1,:)); nbcells = squeeze(all_CI_mat(2,2,:));
    bcellsu = squeeze(uni(1,2,:)); nbcellsu = squeeze(uni(2,1,:));
    
    aVar = [aVar bcellsu intcells nbcellsu];
    
    [within,dvn,xlabels] = make_within_table({'CT','PT'},[2,3]);
    dataT = make_between_table({aVar},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    print_for_manuscript(ra)
    
     %% combined population type and pooled exc and inh cell type RM anova test for on off arb
    good_FRV = all_gV(:,[1 3]);
     [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(good_FRV,0.5,0.05);
    intcells = squeeze(all_CI_mat(1,2,:));    bcells = squeeze(all_CI_mat(1,1,:)); nbcells = squeeze(all_CI_mat(2,2,:));
    bcellsu = squeeze(uni(1,2,:)); nbcellsu = squeeze(uni(2,1,:));
    
    aVar = [bcellsu intcells nbcellsu];
    
    good_FRV = all_gV(:,[2 4]);
     [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(good_FRV,0.5,0.05);
    intcells = squeeze(all_CI_mat(1,2,:));    bcells = squeeze(all_CI_mat(1,1,:)); nbcells = squeeze(all_CI_mat(2,2,:));
    bcellsu = squeeze(uni(1,2,:)); nbcellsu = squeeze(uni(2,1,:));
    
    aVar = [aVar bcellsu intcells nbcellsu];
    
    good_FRV = all_gV(:,[5 6]);
     [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(good_FRV,0.5,0.05);
    intcells = squeeze(all_CI_mat(1,2,:));    bcells = squeeze(all_CI_mat(1,1,:)); nbcells = squeeze(all_CI_mat(2,2,:));
    bcellsu = squeeze(uni(1,2,:)); nbcellsu = squeeze(uni(2,1,:));
    
    aVar = [aVar bcellsu intcells nbcellsu];
    
    [within,dvn,xlabels] = make_within_table({'ET','PT'},[3,3]);
    dataT = make_between_table({aVar},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    print_for_manuscript(ra)
    %%
    % for different VENN Diagrams change good_FRV = all_exc(:,[1 3]) or
    % all_exc(:,[2 4]) or all_exc(:,[5 6]) same thing for all_inh
%     good_FRV = [cell_list_op(all_gV(:,[1:2 5]),[],'or',1) cell_list_op(all_gV(:,[3:4 6]),[],'or',1)]; 
    good_FRV = avar
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
%     good_FR = [all_exc(:,1),all_inh(:,1),all_exc(:,2),all_inh(:,2),all_exc(:,5),all_inh(:,5),all_exc(:,3),all_inh(:,3),all_exc(:,4),all_inh(:,4),all_exc(:,6),all_inh(:,6)]; % for comparison of air onset offset and arb
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(good_FR,0.5,0.05);
    mUni = nanmean(uni,3); mUni1 = tril(mUni,-1) + tril(mUni,-1)'; mUni2 = triu(mUni,1) + triu(mUni,1)'; mmUni = min(mUni(:)); MmUni = max(mUni(:));
    mOI = mCI; semOI = semCI;
%     mOI = mUni2;
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
    txl = {' B-AOn-Exc',' B-AOn-Inh',' B-AOff-Exc',' B-AOff-Inh','NB-AOn-Exc','NB-AOn-Inh','NB-AOff-Exc','NB-AOff-Inh'};
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
    for noth = 1
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
    
    end
    
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
%     txl = {' B-AOn-Exc',' B-AOn-Inh',' B-AOff-Exc',' B-AOff-Inh',' B-Arb-Exc',' B-Arb-Inh','NB-AOn-Exc','NB-AOn-Inh','NB-AOff-Exc','NB-AOff-Inh','NB-Arb-Exc','NB-Arb-Inh'};
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
    ind = 1;
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
    {[Ab_On Abs_On Ab_Off Abs_Off Ab_Offc Abs_Offc Lb Lbs];[Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off Ar_Offc ArL_Offc Ars_Offc ArL_L];[M_On M_Off]}};
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
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
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
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
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
    set(hf,'Position',[7 3 2 1]);
%     set(hf,'Position',[7 3 2 1]);
    set(H,'linewidth',0.5);
    txl = {' B-AOn-Exc',' B-AOn-Inh',' B-AOn-Unt',' B-AOff-Exc',' B-AOff-Inh',' B-AOff-Unt','NB-AOn-Exc','NB-AOn-Inh','NB-AOn-Unt','NB-AOff-Exc','NB-AOff-Inh','NB-AOff-Unt'};
    txl = {' B-AOn-Exc',' B-AOn-Inh',' B-AOff-Exc',' B-AOff-Inh','NB-AOn-Exc','NB-AOn-Inh','NB-AOff-Exc','NB-AOff-Inh'};
%     txl = {' B-AOn',' B-AOff',' B-LOn','NB-AOn','NB-AOff','NB-LOn'};
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel({'Eucledian','Distance'});changePosition(hx,[0 -0.3 0]);
    changePosition(gca,[0.05 0.0 0.0 -0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);
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

%% agglomerative hierarchical clustering
while 1
    mOI1 = mOI;
%     mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1,@naneucdist);
%     tree = linkage(mOI1,'average','euclidean');
    tree = linkage(Di,'average');
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 1 1]);
%     set(hf,'Position',[7 3 2 1]);
    set(H,'linewidth',0.5);
    txl = {' B-AOn-Exc',' B-AOn-Inh',' B-AOff-Exc',' B-AOff-Inh',' B-LOn-Exc',' B-LOn-Inh','NB-AOn-Exc','NB-AOn-Inh','NB-AOff-Exc','NB-AOff-Inh','NB-LOn-Exc','NB-LOn-Inh'};
    txl = {' B-A',' B-L','NB-A','NB-L'};
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel({'Eucledian','Distance'});changePosition(hx,[0 -0.3 0]);
    changePosition(gca,[0.2 0.0 -0.07 -0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);
    %%
    break;
end


%% Brake vs No-Brake Air on and Air off and light ON Untuned Cells
while 1
    %%
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    sic = {[Ab_On Abs_On Ab_Off Abs_Off];[Lb Lbs];[Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off];[ArL_L]};
    clear all_gV all_gFR all_exc all_inh all_good_FR;
    for ii = 1:length(sic)
        sit = sic{ii};
        tRs = o.Rs(:,sit);
        props1 = get_props_Rs(tRs,ntrials);
        gFR = props1.vals;
        gFR = cell_list_op(gFR,[],'or',1);
        all_gV(:,ii) = gFR;
        gFR = props1.exc;
        gFR = cell_list_op(gFR,[],'or',1);
        all_exc(:,ii) = gFR;
        gFR = props1.inh;
        gFR = cell_list_op(gFR,[],'or',1);
        all_inh(:,ii) = gFR;
        gFR = props1.good_FR_and_untuned;
        gFR = cell_list_op(gFR,[],'or',1);
        all_gFR(:,ii) = gFR;
%         gFR = cell_list_op(props1.good_FR,props1.vals,'or'); %gFR = props1.good_FR_and_inh;
%         gFR = cell_list_op(gFR,[],'or',1);
%         all_good_FR(:,ii) = gFR;
        
    end
    good_FR = [all_gV all_gFR];%[all_exc(:,1),all_inh(:,1),all_gFR(:,1),all_exc(:,2),all_inh(:,2),all_gFR(:,2),all_exc(:,3),all_inh(:,3),all_gFR(:,3),all_exc(:,4),all_inh(:,4),all_gFR(:,4)];
    good_FR_bnb = [cell_list_op(all_gFR(:,1:2),[],'or',1) cell_list_op(all_gFR(:,3:4),[],'or',1)];
    %%
    tcolors = mData.colors;
    hf = get_figure(6,[10 7 1.5 1]);
    good_FRV = good_FR_bnb; %good_FRV = all_gFR(:,[1 3]);
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
    changePosition(gca,[0.1 0.01 -0.3 0]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all3.pdf',ntrials),600);
    
      
    %%
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(all_gFR,0.5,0.05);
    mOI = mCI; semOI = semCI;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
%     mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
%     txl = {'B-A-On-Exc','B-A-On-Inh','B-A-On-Unt','B-A-Off-Exc','B-A-Off-Inh','B-A-Off-Unt','NB-A-On-Exc','NB-A-On-Inh','NB-A-On-Unt','NB-A-Off-Exc','NB-A-Off-Inh','NB-A-Off-Unt'};%,'M-On','M-Off'};
    txl = {'Air','Light'}; txl = repmat(txl,1,4);
%     txl = {'Air-On','Air-Off'};txl = repmat(txl,1,2);
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
%     imAlpha(mask1 == 1) = 0;
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
    changePosition(gca,[0.02 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    colormap jet
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
    set(hf,'Position',[7 3 2 1.5]);
%     set(hf,'Position',[7 3 2 1]);
    set(H,'linewidth',0.5);
    txl = {' B-AOn-Exc',' B-AOn-Inh',' B-AOff-Exc',' B-AOff-Inh',' B-LOn-Exc',' B-LOn-Inh','NB-AOn-Exc','NB-AOn-Inh','NB-AOff-Exc','NB-AOff-Inh','NB-LOn-Exc','NB-LOn-Inh'};
    txl = {' Tun-B-A',' Tun-B-L','Tun-NB-A','Tun-NB-L',' Unt-B-A',' Unt-B-L','Unt-NB-A','Unt-NB-L'};
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel({'Eucledian','Distance'});changePosition(hx,[0 -0.3 0]);
    changePosition(gca,[0.05 0.0 0.0 -0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);
    %%
    break;
end

%% Brake vs No-Brake Air on and Air off and light ON Untuned Cells
while 1
    %%
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    sic = {[Ab_On Abs_On];[Ab_Off Abs_Off];[Ar_On ArL_On Ars_On];[Ar_Off ArL_Off Ars_Off]};
    clear all_gV all_gFR all_exc all_inh all_good_FR;
    for ii = 1:length(sic)
        sit = sic{ii};
        tRs = o.Rs(:,sit);
        props1 = get_props_Rs(tRs,ntrials);
        gFR = props1.vals;
        gFR = cell_list_op(gFR,[],'or',1);
        all_gV(:,ii) = gFR;
        gFR = props1.exc;
        gFR = cell_list_op(gFR,[],'or',1);
        all_exc(:,ii) = gFR;
        gFR = props1.inh;
        gFR = cell_list_op(gFR,[],'or',1);
        all_inh(:,ii) = gFR;
        gFR = props1.good_FR_and_untuned;
        gFR = cell_list_op(gFR,[],'or',1);
        all_gFR(:,ii) = gFR;
%         gFR = cell_list_op(props1.good_FR,props1.vals,'or'); %gFR = props1.good_FR_and_inh;
%         gFR = cell_list_op(gFR,[],'or',1);
%         all_good_FR(:,ii) = gFR;
        
    end
    good_FR = [all_gV all_gFR];%[all_exc(:,1),all_inh(:,1),all_gFR(:,1),all_exc(:,2),all_inh(:,2),all_gFR(:,2),all_exc(:,3),all_inh(:,3),all_gFR(:,3),all_exc(:,4),all_inh(:,4),all_gFR(:,4)];
    good_FR_bnb = [cell_list_op(all_gFR(:,1:2),[],'or',1) cell_list_op(all_gFR(:,3:4),[],'or',1)];
    
    %%
    a3c = [all_exc(:,1) all_inh(:,1) all_exc(:,2) all_inh(:,2) all_exc(:,3) all_inh(:,3) all_exc(:,4) all_inh(:,4)];
    a3 = 100*exec_fun_on_cell_mat(a3c,'sum')./exec_fun_on_cell_mat(a3c,'length');
    [within,dvn,xlabels] = make_within_table({'B','E','C'},[2,2,2]);
    dataT = make_between_table({a3},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'B_E_C','hsd'},[1.5 1 1]);
%     h(h==1) = 0;
	xdata = make_xdata([2 2 2 2],[1 1.5]);
    hf = get_figure(5,[8 7 5 2]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',7,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
%     set(hbs(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); set(hbs(3),'FaceColor',tcolors{2},'FaceAlpha',0.75);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0,maxY+6,sprintf('Any Condition (%d\x00B1%d%%),   All Conditions (%d%%)',round(mra),round(semra),round(mrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'B-On-Exc','B-On-Inh','B-Off-Exc','B-Off-Inh','NB-On-Exc','NB-On-Inh','NB-Off-Exc','NB-Off-Inh'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.1 0.01 -0.3 0]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    
    %%
    tcolors = mData.colors;
    hf = get_figure(6,[10 7 1.5 1]);
    good_FRV = good_FR_bnb; good_FRV = all_gFR(:,[2 4]);
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
    changePosition(gca,[0.1 0.01 -0.3 0]); put_axes_labels(gca,{[],[0 0 0]},{{'Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('active_cells_across_conditions_%d_all3.pdf',ntrials),600);
    
      
    %%
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(all_gFR,0.5,0.05);
    mOI = mCI; semOI = semCI;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
%     mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
%     txl = {'B-A-On-Exc','B-A-On-Inh','B-A-On-Unt','B-A-Off-Exc','B-A-Off-Inh','B-A-Off-Unt','NB-A-On-Exc','NB-A-On-Inh','NB-A-On-Unt','NB-A-Off-Exc','NB-A-Off-Inh','NB-A-Off-Unt'};%,'M-On','M-Off'};
    txl = {'Air','Light'}; txl = repmat(txl,1,4);
%     txl = {'Air-On','Air-Off'};txl = repmat(txl,1,2);
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
%     imAlpha(mask1 == 1) = 0;
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
    changePosition(gca,[0.02 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',ntrials),600);
   
    %%
    break;
end


%% compare the firing rate
while 1
    ntrials = 50; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    si = [Lb Ab_On Ab_Off Ar_On Ar_Off ArL_On ArL_L ArL_Off Ars_On Ars_Off Lbs Abs_On Abs_Off];
    Rs = o.Rs(:,si); Rs_MC = o.RsMC(:,si);
    props1 = get_props_Rs(Rs,ntrials);
%     good_FR = props1.vals;
    all_zMIs = props1.mean_FR;
    zMIs = exec_fun_on_cell_mat(all_zMIs,'nanmean');
%     azMIs = zMIs
    [within,dvn,xlabels] = make_within_table({'Cond'},[13]);
    dataT = make_between_table({zMIs},dvn);
    ra = RMA(dataT,within);
    ra.ranova
%     ra.mauchly
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = make_xdata([13],[1 1.5]);
%     xdata = make_xdata([2],[1 1.5]);
    ptab = 0;
    if ptab h(h==1) = 0; end
    hf = get_figure(5,[8 7 3.5 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;

    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.02,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',-0.01);
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = rasterNamesTxt(si);
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45);
    changePosition(gca,[0.01 -0.01 0.06 0.01]); put_axes_labels(gca,{[],[0 0 0]},{'Firing Rate (A.U.)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('zMIs_good_FR_all_Conditions.pdf'),600);
    %%
    break;
end
