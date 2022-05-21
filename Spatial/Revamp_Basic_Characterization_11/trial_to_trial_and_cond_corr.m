function trial_to_trial_corr

%% light sensory trial to trial correlation
while 1
    si = [Lb_T Lbs_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    exc = props1.good_FR_and_exc;     sup = props1.good_FR_and_inh;    com = props1.good_FR_and_untuned;
    trials = mat2cell([1:10]',ones(size([1:10]')));
    parfor ii = 1:size(Rs,2)
        outTrialsExc{ii} = find_population_vector_corr_remap_trials(Rs(:,ii),exc(:,ii),trials);
        outTrialsCom{ii} = find_population_vector_corr_remap_trials(Rs(:,ii),com(:,ii),trials);
    end
    %%
    break;
end

%% light trial to trial corr
while 1
    meancorr_trials = [];
    for ii = 1:2
        toutTrials = outTrialsExc{ii};
        meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(toutTrials.adj_SP_corr_diag,'nanmean')];
    end
    for ii = 1:2
        toutTrials = outTrialsCom{ii};
        meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(toutTrials.adj_SP_corr_diag,'nanmean')];
    end
    [within,dvn,xlabels] = make_within_table({'CT','Cond','TP'},[2,2,9]);
    dataT = make_between_table({meancorr_trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors([1 3]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Exc','Com'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.5 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Trial Pair','Correlation(%)'},[0 0 0]});
    ht = title('Lb and Lb* (pooled)'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('light_trial_to_trial_corr.pdf'),600);
    %%
    break;
end

%% light rest across condition corr
while 1
    %%
    si = [Lb_T Lbs_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    exc = props1.good_FR_and_exc;     com = props1.good_FR_and_untuned;
    respE_OR = cell_list_op(exc,[],'or'); outE = find_population_vector_corr_remap(Rs,mR,respE_OR);
    respC_OR = cell_list_op(com,[],'or'); outC = find_population_vector_corr_remap(Rs,mR,respC_OR);
    disp('Done');
    si = [ArL_L_T];
    Rsm = o.Rs(:,si);mRm = o.mR(:,si);
    ntrials = 50;
    props1m = get_props_Rs(Rsm,ntrials);
    excm = props1m.good_FR_and_exc; supm = props1m.good_FR_and_inh;    comm = props1m.good_FR_and_untuned;
    disp('Done');
    n = 0;
    %%
    break;
end
    
%% light across conditions corr temporal
while 1
    meancorr_trials = [];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outE.adj_SP_corr_diag,'nanmean')];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outC.adj_SP_corr_diag,'nanmean')];
    [within,dvn,xlabels] = make_within_table({'CT'},[2]);
    dataT = make_between_table({meancorr_trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors([1 3]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Exc','Com'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.5 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Temporal','Correlation(%)'},[0 0 0]});
    ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('Lb_Lbs_corr.pdf'),600);
    %%
    break;
end

%% light across conditions corr PV
while 1
    meancorr_trials = [];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outE.adj_PV_corr_diag,'nanmean')];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outC.adj_PV_corr_diag,'nanmean')];
    [within,dvn,xlabels] = make_within_table({'CT'},[2]);
    dataT = make_between_table({meancorr_trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors([1 3]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Exc','Com'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.5 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Temporal','Correlation(%)'},[0 0 0]});
    ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('Lb_Lbs_corr_PV.pdf'),600);
    %%
    break;
end

%% light across conditions corr RR
while 1
    meancorr_trials = [];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outE.adj_RR_SP,'nanmean')];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outC.adj_RR_SP,'nanmean')];
    [within,dvn,xlabels] = make_within_table({'CT'},[2]);
    dataT = make_between_table({meancorr_trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors([1 3]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Exc','Com'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.5 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Temporal','Correlation(%)'},[0 0 0]});
    ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('Lb_Lbs_corr_PV.pdf'),600);
    %%
    break;
end



%% air sensory trial to trial correlation
while 1
    titles = {'Ab','Ab*'};
    si = [Ab_T Abs_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    exc = props1.good_FR_and_exc;     sup = props1.good_FR_and_inh;    com = props1.good_FR_and_untuned;
    trials = mat2cell([1:10]',ones(size([1:10]')));
    parfor ii = 1:size(Rs,2)
        outTrialsExc{ii} = find_population_vector_corr_remap_trials(Rs(:,ii),exc(:,ii),trials);
        outTrialsSup{ii} = find_population_vector_corr_remap_trials(Rs(:,ii),sup(:,ii),trials);
        outTrialsCom{ii} = find_population_vector_corr_remap_trials(Rs(:,ii),com(:,ii),trials);
    end
    
    %%
    break;
end

%% air trial to trial corr
while 1
    meancorr_trials = [];
    for ii = 1:2
        toutTrials = outTrialsExc{ii};
        meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(toutTrials.adj_SP_corr_diag,'nanmean')];
    end
    for ii = 1:2
        toutTrials = outTrialsSup{ii};
        meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(toutTrials.adj_SP_corr_diag,'nanmean')];
    end
    for ii = 1:2
        toutTrials = outTrialsCom{ii};
        meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(toutTrials.adj_SP_corr_diag,'nanmean')];
    end
    [within,dvn,xlabels] = make_within_table({'CT','Cond','TP'},[3,2,9]);
    dataT = make_between_table({meancorr_trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([3],[1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Exc','Sup','Com'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.4 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Trial Pair','Correlation(%)'},[0 0 0]});
    ht = title('Ab and Ab* (pooled)'); changePosition(ht,[-1.25 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('air_trial_to_trial_corr.pdf'),600);
    %%
    break;
end

%% air rest across condition corr
while 1
    %%
    titles = {'Ab','Ab*'};
    si = [Ab_T Abs_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    exc = props1.good_FR_and_exc;     sup = props1.good_FR_and_inh;    com = props1.good_FR_and_untuned;
    respE_OR = cell_list_op(exc,[],'or'); outE = find_population_vector_corr_remap(Rs,mR,respE_OR);
    respS_OR = cell_list_op(sup,[],'or'); outS = find_population_vector_corr_remap(Rs,mR,respS_OR);
    respC_OR = cell_list_op(com,[],'or'); outC = find_population_vector_corr_remap(Rs,mR,respC_OR);
    disp('Done');
    n = 0;
    %%
    break;
end
    
%% air across conditions corr
while 1
    meancorr_trials = [];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outE.adj_SP_corr_diag,'nanmean')];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outS.adj_SP_corr_diag,'nanmean')];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outC.adj_SP_corr_diag,'nanmean')];
    [within,dvn,xlabels] = make_within_table({'CT'},[3]);
    dataT = make_between_table({meancorr_trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([3],[1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Exc','Sup','Com'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.4 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Temporal','Correlation(%)'},[0 0 0]});
    ht = title('Across Ab and Ab*'); changePosition(ht,[-1.25 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('Ab_Abs_corr.pdf'),600);
    %%
    break;
end


%% Overlap Indices ImageSC
while 1
    ntrials = 50;
    si_nb = [Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
    props_nb = get_props_Rs(o.Rs(:,si_nb),ntrials);
    resp_nb = props_nb.good_FR;
    si = [Ab_T Abs_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    exc = props1.good_FR_and_exc;     sup = props1.good_FR_and_inh;    com = props1.good_FR_and_untuned;
    
    respE_OR = cell_list_op(exc,[],'or');
    respS_OR = cell_list_op(sup,[],'or');
    respC_OR = cell_list_op(com,[],'or'); 
    
    resp = [respE_OR(:,1) respS_OR(:,1) respC_OR(:,1) resp_nb];
%     resp = [respE_OR(:,1) respS_OR(:,1) respC_OR(:,1)];
%     resp = [exc(:,1), sup(:,1), com(:,1)];
%     
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = [{'Exc','Sup','Com'} rasterNamesTxt(si_nb)]; 
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
    text(-0.3,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_exc_nb.pdf'),600);
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
    [H,T,TC] = dendrogram(tree,'Orientation','right','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 1.25 1.25]);
    set(H,'linewidth',1);
    set(gca,'yticklabels',txl(TC));ytickangle(30);
    format_axes(gca);
    hx = xlabel('Eucledian Distance');%changePosition(hx,[-0.051 0 0]);
    changePosition(gca,[0 0.0 0.01 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster_exc_nb.pdf'),600);
    %%
    break;
end



%% Overlap Indices ImageSC light
while 1
    ntrials = 50;
    si_nb = [Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
    props_nb = get_props_Rs(o.Rs(:,si_nb),ntrials);
    resp_nb = props_nb.good_FR;
    resp = [respE_OR(:,1) respC_OR(:,1) resp_nb];
    resp = [respE_OR(:,1) respC_OR(:,1) resp_nb excm supm comm];
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = [{'Exc','Com'} rasterNamesTxt(si_nb) {'ExcM','SupM','ComM'}]; 
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
    text(-0.3,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.0 0 -0.04 0]);
    hc = putColorBar(gca,[0.03 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_exc_nb_light.pdf'),600);
    %%
    break;
end

%% light agglomerative hierarchical clustering
while 1
    mOI1 = mOI;
    mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1);
    tree = linkage(Di);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','right','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 1.25 1.25]);
    set(H,'linewidth',1);
    set(gca,'yticklabels',txl(TC));ytickangle(30);
    format_axes(gca);
    hx = xlabel('Eucledian Distance');%changePosition(hx,[-0.051 0 0]);
    changePosition(gca,[0 0.0 -0.05 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster_exc_nb_light.pdf'),600);
    %%
    break;
end




%% Overlap Indices ImageSC sensory responses immobility versus locomotion
while 1
    ntrials = 50;
    si_nb = [Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
    props_nb = get_props_Rs(o.Rs(:,si_nb),ntrials);
    resp_nb = props_nb.good_FR;
    
    si = [Lb_T Lbs_T Ab_T Abs_T ArL_L_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    exc = props1.good_FR_and_exc;     sup = props1.good_FR_and_inh;    com = props1.good_FR_and_untuned;
    
    resp = [exc(:,1:4) sup(:,1:4) com(:,1:4) resp_nb exc(:,5) sup(:,5) com(:,5)];
    
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = [rasterNamesTxt([si(1:4) si(1:4) si(1:4) si_nb si(5) si(5) si(5)])]; 
    for ii = 1:12
        if ii >=1 && ii <=4
            txl{ii} = sprintf('%s-E',txl{ii});
        end
        if ii >=5 && ii <=8
            txl{ii} = sprintf('%s-S',txl{ii});
        end
        if ii >=9 && ii <=12
            txl{ii} = sprintf('%s-C',txl{ii});
        end
    end
    ii = 19; txl{ii} = sprintf('%s-E',txl{ii}); ii = 20; txl{ii} = sprintf('%s-S',txl{ii}); ii = 21; txl{ii} = sprintf('%s-C',txl{ii}); 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(5,[8 2 5 5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    text(-0.3,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.0 0 -0.04 0]);
    hc = putColorBar(gca,[0.1 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_exc_nb_light_air_motion.pdf'),600);
    %%
    break;
end

%% light agglomerative hierarchical clustering
while 1
    mOI1 = mOI;
    mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1);
    tree = linkage(Di);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 5 3]);
    set(H,'linewidth',1);
    set(gca,'xticklabels',txl(TC));xtickangle(30);
    format_axes(gca);
    hx = ylabel('Eucledian Distance');%changePosition(hx,[-0.051 0 0]);
    changePosition(gca,[0 0.0 -0.05 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster_exc_nb_light_air_motion.pdf'),600);
    %%
    break;
end



%% spatial tuned untuned zMI
while 1
    si = [Ar_t_D ArL_t_D Ars_t_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    gauss = props1.good_FR_and_Gauss; n_gauss = props1.good_FR_and_notGauss;
    zMIs = props1.zMI;
    
    for rr = 1:size(zMIs,1)
        for cc = 1:size(zMIs,2)
            tzMI = zMIs{rr,cc};
            tgauss = gauss{rr,cc}; t_n_gauss = n_gauss{rr,cc};
            all_zMIsG(rr,cc) = nanmean(tzMI(tgauss));
            all_zMIsnG(rr,cc) = nanmean(tzMI(t_n_gauss));
        end
    end
    [within,dvn,xlabels] = make_within_table({'CT','Cond'},[2,3]);
    dataT = make_between_table({all_zMIsG,all_zMIsnG},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'gT','gU'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.5 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Mutual Information','(z-score)'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('spatial_across_conditions_zMI.pdf'),600);

    
    %%
    break;
end

%% find spatial trial to trial correlation
while 1
    si = [Ar_t_D ArL_t_D Ars_t_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    gauss = props1.good_FR_and_Gauss_loose; n_gauss = props1.good_FR_and_notGauss_loose;
    trials = mat2cell([1:10]',ones(size([1:10]')));
    parfor ii = 1:size(Rs,2)
        outTrialsG{ii} = find_population_vector_corr_remap_trials(Rs(:,ii),gauss(:,ii),trials);
        outTrials_nG{ii} = find_population_vector_corr_remap_trials(Rs(:,ii),n_gauss(:,ii),trials);
    end
    %%
    break;
end

%% spatial trial to trial corr
while 1
    meancorr_trials = [];
    for ii = 1:3
        toutTrials = outTrialsG{ii};
        meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(toutTrials.adj_SP_corr_diag,'nanmean')];
    end
    for ii = 1:3
        toutTrials = outTrials_nG{ii};
        meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(toutTrials.adj_SP_corr_diag,'nanmean')];
    end
    [within,dvn,xlabels] = make_within_table({'CT','Cond','TP'},[2,3,9]);
    dataT = make_between_table({meancorr_trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1 2]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'gT','gU'};%{'Trials 1-2','Trials 2-3','Trials 3-4','Trials 4-5','Trials 5-6','Trials 6-7','Trials 7-8','Trials 8-9','Trials 9-10'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.5 -0.07]); put_axes_labels(gca,{[],[0 0 0]},{{'Trial Pair','Correlation'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('spatial_trial_to_trial_corr.pdf'),600);
     %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([9 9 9 9 9 9],[1 2]);
    hf = get_figure(5,[8 7 2 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = [colors colors];
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Trials 1-2','Trials 2-3','Trials 3-4','Trials 4-5','Trials 5-6','Trials 6-7','Trials 7-8','Trials 8-9','Trials 9-10'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.1 0.0 -0.0 -0.07]); put_axes_labels(gca,{[],[0 0 0]},{{'Avg. Trial-to-Trial','Correlation(%)'},[0 0 0]});
    ht = title('Pooled across conditions and cell types'); changePosition(ht,[-1 0.01 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('spatial_trial_to_trial_corr.pdf'),600);
    %%
    break;
end

%% find spatial across condition corr
while 1
    %%
    si = [Ar_t_D ArL_t_D Ars_t_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    gauss = props1.good_FR_and_Gauss; n_gauss = props1.good_FR_and_notGauss;
    respE_OR = cell_list_op(gauss,[],'or'); outT = find_population_vector_corr_remap(Rs,mR,respE_OR);
    respC_OR = cell_list_op(n_gauss,[],'or'); outU = find_population_vector_corr_remap(Rs,mR,respC_OR);
    respR_OR = cell_list_op(props1.good_FR,[],'or'); outR = find_population_vector_corr_remap(Rs,mR,respR_OR);
    disp('Done');
    n = 0;
    %%
    break;
end
    
%% spatial across conditions corr
while 1
    meancorr_trials = [];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outT.adj_SP_corr_diag,'nanmean')];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outU.adj_SP_corr_diag,'nanmean')];
    [within,dvn,xlabels] = make_within_table({'CT','Cond'},[2,2]);
    dataT = make_between_table({meancorr_trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'gT','gU'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.5 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Spatial','Correlation(%)'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('spatial_across_conditions_corr.pdf'),600);
    %%
    break;
end

%% spatial across conditions corr PV
while 1
    meancorr_trials = [];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outR.adj_cell_pos_shift,'nanmean')];
    [within,dvn,xlabels] = make_within_table({'CT','Cond'},[2,2]);
    dataT = make_between_table({meancorr_trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1]);
    hf = get_figure(5,[8 7 1.25 1]);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'gT','gU'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.5 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Pop. Vec. ','Correlation(%)'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('spatial_across_conditions_corr_PV.pdf'),600);
    %%
    break;
end


%% spatial across conditions corr PV
while 1
    meancorr_trials = [];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outT.adj_PV_corr_diag,'nanmean')];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outU.adj_PV_corr_diag,'nanmean')];
    [within,dvn,xlabels] = make_within_table({'CT','Cond'},[2,2]);
    dataT = make_between_table({meancorr_trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1]);
    hf = get_figure(5,[8 7 1.25 1]);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'gT','gU'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.5 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Pop. Vec. ','Correlation(%)'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('spatial_across_conditions_corr_PV.pdf'),600);
    %%
    break;
end

%% spatial across conditions corr SP _ cell ordering difference
while 1
    meancorr_trials = [];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outR.adj_PV_corr_diag,'nanmean')];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outR.adj_PV_corr_diag_ind,'nanmean')];
    [within,dvn,xlabels] = make_within_table({'C_ord','Cond'},[2,2]);
    dataT = make_between_table({meancorr_trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'C_ord','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1 2]);
    hf = get_figure(5,[8 7 1.25 1]);
    tcolors = colors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Same Order','Peak Order'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.5 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Spatial','Correlation(%)'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('spatial_across_conditions_corr_PV_ind.pdf'),600);
    %%
    break;
end


%% spatial across conditions corr RR
while 1
    meancorr_trials = [];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outT.adj_RR_SP,'nanmean')];
    meancorr_trials = [meancorr_trials exec_fun_on_cell_mat(outU.adj_RR_SP,'nanmean')];
    [within,dvn,xlabels] = make_within_table({'CT','Cond'},[2,2]);
    dataT = make_between_table({meancorr_trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([2],[1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.025,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
%     maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[0.8 0.93]); format_axes(gca);
    xticks = xdata; xticklabels = {'gT','gU'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.2 0.0 -0.5 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'Rate Remapping','Score'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('spatial_across_conditions_corr_RR.pdf'),600);
    %%
    break;
end

%% Overlap Indices ImageSC spatial
while 1
    ntrials = 50;
    si = [Ab_T Abs_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    propsB = get_props_Rs(Rs,ntrials);
    exc = propsB.good_FR_and_exc;     sup = propsB.good_FR_and_inh;    com = propsB.good_FR_and_untuned;
    respE_OR = cell_list_op(exc,[],'or');     respS_OR = cell_list_op(sup,[],'or');     respC_OR = cell_list_op(com,[],'or');
    
    si_nb = [Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
    props_nb = get_props_Rs(o.Rs(:,si_nb),ntrials);
    resp_nb = props_nb.good_FR;
    
    si = [Ar_t_D ArL_t_D Ars_t_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    propsAt = get_props_Rs(Rs,ntrials);
    gauss = propsAt.good_FR_and_Gauss; n_gauss = propsAt.good_FR_and_notGauss;
    respG_OR = cell_list_op(gauss,[],'or'); respNG_OR = cell_list_op(n_gauss,[],'or');
    
    resp = [respE_OR(:,1) respS_OR(:,1) respC_OR(:,1) gauss n_gauss];
    resp = [gauss n_gauss resp_nb(:,4:6)];
    [OI,mOI,semOI,OI_mat,p_vals,h_vals] = get_overlap_index(resp,0.5,0.05);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:);semOI(:)]);    
    minI = min([mOI(:);semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = [{'Exc','Sup','Com'} {'gT','gT','gT','gU','gU','gU'}]; 
    txl = [{'gT3','gT4','gT5','gU3','gU4','gU5'} rasterNamesTxt(si_nb(4:6))]; 
%     mOI = mOI .* mask;
    imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
    imAlpha(mask1 == 1) = 0;
%     ff = makeFigureRowsCols(2020,[10 4 6 1.5],'RowsCols',[1 2],'spaceRowsCols',[0.1 0.01],'rightUpShifts',[0.05 0.13],'widthHeightAdjustment',[-240 -150]);
    hf = get_figure(6,[8 7 1.5 1.5]);
    %
%     axes(ff.h_axes(1));
    im1 = imagesc(mOI,[minI,maxI]);    im1.AlphaData = imAlpha;
    set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
    format_axes(gca);
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    text(-0.3,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[0.03 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%d',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.06 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_tu_spatial.pdf'),600);
    %%
    break;
end

%% spatial agglomerative hierarchical clustering
while 1
    mOI1 = mOI;
    mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1);
    tree = linkage(Di);
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','right','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 1.25 1.25]);
    set(H,'linewidth',1);
    set(gca,'yticklabels',txl(TC));ytickangle(30);
    format_axes(gca);
    hx = xlabel('Eucledian Distance');%changePosition(hx,[-0.051 0 0]);
    changePosition(gca,[0 0.0 -0.05 0.05]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster_tu_spatial.pdf'),600);
    %%
    break;
end


