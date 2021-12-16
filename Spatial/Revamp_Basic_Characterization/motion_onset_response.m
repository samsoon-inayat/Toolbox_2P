function firing_rate_motion_vs_rest
while 1
    titles = {'MOn-T','MOff-T'};
    si = [MOn_T MOff_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    %%
    break;
end
%% See onsets and offsets
while 1
    an = 3; cn = 1;
    motionOnsets = Rs{an,cn}.onsets;
    motionOffsets = Rs{an,cn}.offsets;
    hf = figure(1000);clf;
    set(gcf,'color','w'); set(gcf,'Position',[10 4 1.25 1]);
    display_with_air_puff(ei{an}.b,motionOnsets,motionOffsets);
    xlim([270 330]/60);ylim([0 1.1]);
    changePosition(gca,[-0.07 0.15 0 -0.2]); box off;
    ax = gca; ax.YAxis.Visible = 'off';
    xlabel('Time (min)');
    set(gca,'FontSize',6,'FontWeight','Normal','TickDir','out');
    save_pdf(hf,mData.pdf_folder,sprintf('motionOnset_trials'),600);
    break;
end
%% Show sample rasters
an = 3; cn = 1;
figure(2000);clf;imagesc(Rs{an,cn}.speed);colorbar;
plotRasters_simplest(Rs{an,cn})
%% population vector and correlation
while 1
    selected_property = 'tuned';
    an = 4;
    resp = props1.vals;
%     eval(cmdTxt);
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 2],...
        'spaceRowsCols',[0.01 0.05],'rightUpShifts',[0.16 0.11],'widthHeightAdjustment',...
        [-130 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 1.7 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    colormap parula
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_motion_%s.pdf',selected_property),600);

    % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 2],...
        'spaceRowsCols',[0.01 0.05],'rightUpShifts',[0.16 0.25],'widthHeightAdjustment',...
        [-130 -320]);    set(gcf,'color','w');    set(gcf,'Position',[10 8 1.7 0.8]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    colormap parula
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_motion_%s.pdf',selected_property),600);
    %%
    break;
end

%% number of trials
while 1
    for rr = 1:size(Rs,1)
        for cc = 1
            tRs = Rs{rr,cc};
            num_trials(rr,cc) = size(tRs.sp_rasters1,1);
        end
    end
    %%
    break;
end

%% average % of trials in which the cell responded
while 1
    mean_N_trials_resp = exec_fun_on_cell_mat(props1.N_Resp_Trials,'mean');
    [within,dvn,xlabels] = make_within_table({'Cond'},[length(si)]);
    dataT = make_between_table({mean_N_trials_resp},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'MOn','MOff'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[0.065 0.0 -0.4 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Trials (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('motion_trials.pdf'),600);
    %%
    break;
end
%% Percentage of Responsive Cells
while 1
    resp = props1.vals;
    perc_r = 100*exec_fun_on_cell_mat(resp,'sum')./exec_fun_on_cell_mat(resp,'length');
    within = make_within_table({'Cond'},2);
    dataT = make_between_table({perc_r},{'M_On','M_Off'});
    ra = RMA(dataT,within);
    ra.ranova
    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w'); hold on;
    tcolors = mData.dcolors(17:20);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.05);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
%     hatch(hbs(2),135,'k','-',2,0.1);
    xticks = [xdata(1:end)]; xticklabels = {'Onset','Offset'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30)
    changePosition(gca,[0.2 0.03 -0.5 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Motion Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('percentage_motion_responsive'),600);
    break;
end
%% Percentage of excitatory inhibitory responsive cells
while 1
    resp1 = props1.exc;     resp2 = props1.inh;
    presp1 = 100 * exec_fun_on_cell_mat(resp1,'sum')./exec_fun_on_cell_mat(resp1,'length');
    presp2 = 100 * exec_fun_on_cell_mat(resp2,'sum')./exec_fun_on_cell_mat(resp2,'length');
    
    [within,dvn,xlabels] = make_within_table({'MT','CT'},[2,2]);
    dataT = make_between_table({presp1(:,1),presp2(:,1),presp1(:,2),presp2(:,2)},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'MT_by_CT','hsd'},[1 1 1]);
    xdata = make_xdata([2,2],[1 2]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Exc','Sup'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[0.18 0.0 -0.2 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Motion Responsive','Cells (%)'},[0 0 0]});
    
    save_pdf(hf,mData.pdf_folder,'motion_responsive_exc_inh',600);
    %%
    break;
end


%% Overlap Indices ImageSC 
while 1
    ntrials = 50;
    si = [Lb_T ArL_L_T Lbs_T Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
    props1A = get_props_Rs(o.Rs,ntrials);
    respO = [props1A.good_FR(:,si)];% resp_speed];
    resp = [resp_speed(:,1) props1.vals respO];% resp_speed];
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = [{'M-On','M-Off','Speed'} rasterNamesTxt(si)]; 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 7 2 2]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    text(0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean_motion.pdf',ntrials),600);
    %%
    break;
end

%% agglomerative hierarchical clustering zMIT>zMID
while 1
    mOI1 = mOI;
    mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1);
    tree = linkage(Di);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 2.5 1]);
    set(H,'linewidth',1);
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel({'Eucledian','Distance'});%changePosition(hx,[-0.051 0 0]);
    changePosition(gca,[0 0.0 0.07 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster_motion.pdf'),600);
    %%
    break;
end



%% Overlap Indices ImageSC for presentation
while 1
    ntrials = 50;
    si = [Lb_T ArL_L_T Lbs_T Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D];
    props1A = get_props_Rs(o.Rs,ntrials);
    respO = [props1A.good_FR(:,si)];% resp_speed];
    resp = [props1.vals respO];% resp_speed];
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = [{'M-On','M-Off'} rasterNamesTxt(si)]; 
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
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    text(-1,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean_motion.pdf',ntrials),600);
    %%
    break;
end

%% agglomerative hierarchical clustering 
while 1
    mOI1 = mOI;
    mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1);
    tree = linkage(Di);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 2.5 1]);
    set(H,'linewidth',1);
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel({'Eucledian','Distance'});%changePosition(hx,[-0.051 0 0]);
    changePosition(gca,[0 0.0 0.07 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster_motion.pdf'),600);
    %%
    break;
end

