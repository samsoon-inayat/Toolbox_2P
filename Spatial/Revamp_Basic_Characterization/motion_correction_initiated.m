function motion_correction_initiated


%% compare percent responsive cells
while 1
    ntrials = 40; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    si = [Lb Lbs Ab_On Abs_On Ab_Off Abs_Off ArL_L Ar_On ArL_On Ars_On Ar_Off ArL_Off Ars_Off M_On M_Off];
    props1 = get_props_Rs(o.Rs,ntrials);
    good_FR = props1.good_FR_and_tuned(:,si);
    good_FR = props1.good_zMI_MC(:,si);
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
    ntrials = 50;
    si = [Ab_On Abs_On Ar_On ArL_On Ars_On];
    props1 = get_props_Rs(o.Rs(:,si),ntrials);
    all_resp = props1.good_FR_and_zMI;
    all_resp1 = props1.good_FR_and_zMI_MC;
    resp_OI = [all_resp all_resp1];
%     [OI,mOI,semOI,OI_mat,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat] = get_overlap_index(resp_OI,0.5,0.05);
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(resp_OI,0.5,0.05);
    mOI = mCI; semOI = semCI;
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = {'B-A-Act','NB-A-Act','B-A-Sup','NB-A-Sup'}; 
    txl = rasterNamesTxt([si si]);
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 3.5 3.5]);
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
