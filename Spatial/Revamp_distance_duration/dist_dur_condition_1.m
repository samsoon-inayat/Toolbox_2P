function dist_dur

%%
RsDt = o.Rs(:,[Ar_t_D]);  RsTt = o.Rs(:,[Ar_t_T]);
RsDi = o.Rs(:,[Ar_i_D]);  RsTi = o.Rs(:,[Ar_i_T]);
[dzMI_FD,dzMI_FT] = get_zMI_comp_dist_time(RsDt,RsTt,RsDi,RsTi);

while 0
for ii = 1:5
    tRs = RsTt{ii};
    speed = tRs.speed';
    rasters = tRs.sp_rasters;
    rasters = permute(rasters,[2 1 3]);
    firing = reshape(rasters,size(rasters,1)*size(rasters,2),size(rasters,3));
    speed = reshape(speed,size(rasters,1)*size(rasters,2),1);
    corrfs = corr(firing,speed,'rows','complete');
    speedResp{ii,1} = corrfs > 0.2;
    figure(100);clf;hist(corrfs,50);pause(0.5);
    tRs = RsTi{ii};
    speed = tRs.speed';
    rasters = tRs.sp_rasters;
    rasters = permute(rasters,[2 1 3]);
    firing = reshape(rasters,size(rasters,1)*size(rasters,2),size(rasters,3));
    speed = reshape(speed,size(rasters,1)*size(rasters,2),1);
    corrfs = corr(firing,speed,'rows','complete');
    speedRespI{ii,1} = corrfs > 0.2;
    figure(101);clf;hist(corrfs,50);pause(0.5);
end
end

%
respfids = {[40 100]};

raster_types = {'RsTt','RsDt','RsTi','RsDi'};
% raster_types = {'RsTt','RsTi'};

clear props
for ii = 1:length(respfids)
    trialsR = respfids{ii};
    for jj = 1:length(raster_types)
        cmdTxt = sprintf('props{ii,jj} = get_props_Rs(%s,trialsR);',raster_types{jj});
        eval(cmdTxt);
    end
end
% find dis, dur, and mix cells
while 1
    for ii = 1:length(respfids)
        for jj = 1:length(raster_types)
            cell_list{ii,jj} = props{ii,jj}.good_FR_and_tuned;
        end
    end

    clear FD_Dis_comp FD_Dur_comp FD_conj FT_Dis_comp FT_Dur_comp FT_conj
    for ii = 1:length(respfids)
        FD_Dur = cell_list_op(props{ii,1}.vals,props{ii,1}.good_FR,'and'); FD_Dis = cell_list_op(props{ii,2}.vals,props{ii,2}.good_FR,'and');
        FT_Dur = cell_list_op(props{ii,3}.vals,props{ii,3}.good_FR,'and'); FT_Dis = cell_list_op(props{ii,4}.vals,props{ii,4}.good_FR,'and');
        
        cellP1 = FD_Dis; cellP2 = FD_Dur;
        FD_Dis_comp{ii} = cell_list_op(cellP1,cell_list_op(cellP2,[],'not'),'and');
        FD_Dur_comp{ii} = cell_list_op(cell_list_op(cellP1,[],'not'),cellP2,'and');
        FD_conj{ii} = cell_list_op(cellP1,cellP2,'and');
        
        cellP1 = FT_Dis; cellP2 = FT_Dur;
        FT_Dis_comp{ii} = cell_list_op(cellP1,cell_list_op(cellP2,[],'not'),'and');
        FT_Dur_comp{ii} = cell_list_op(cell_list_op(cellP1,[],'not'),cellP2,'and');
        FT_conj{ii} = cell_list_op(cellP1,cellP2,'and');
%         an = 4; cn = 3; respC = FT_Dur_comp{2}; tempCL = respC{an,cn};
    end
break;
end

%% %% choosing 1 condition and comparing across cell types, response fidelity, trials and intertrials
while 1
    %%
    cni = 1;
    rfi = 1;
    resp = [FD_Dur_comp{rfi}(:,cni) FD_Dis_comp{rfi}(:,cni) FD_conj{rfi}(:,cni) FT_Dur_comp{rfi}(:,cni) FT_Dis_comp{rfi}(:,cni) FT_conj{rfi}(:,cni)];
%     rfi = 2;
%     resp = [resp FD_Dur_comp{rfi}(:,cni) FD_Dis_comp{rfi}(:,cni) FD_conj{rfi}(:,cni) FT_Dur_comp{rfi}(:,cni) FT_Dis_comp{rfi}(:,cni) FT_conj{rfi}(:,cni)];

    per_resp = 100*exec_fun_on_cell_mat(resp,'sum')./exec_fun_on_cell_mat(resp,'length');

    [within,dvn,xlabels] = make_within_table({'TI','CT'},[2,3]);
    dataT = make_between_table({per_resp},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova

    any_cells = cell_list_op(resp,[],'or',1);
    per_active_any = 100*exec_fun_on_cell_mat(any_cells,'sum')./exec_fun_on_cell_mat(any_cells,'length');
    any_cells = cell_list_op(resp,[],'and',1);
    per_active_all = 100*exec_fun_on_cell_mat(any_cells,'sum')./exec_fun_on_cell_mat(any_cells,'length');

    [mra,semra] = findMeanAndStandardError(per_active_any);
    [mrall,semrall] = findMeanAndStandardError(per_active_all);
    
    %%
   
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_CT','hsd'},[1.5 1 1]);
     h(h==1) = 0;
     p = ones(size(p)); 
     ind = ismember(combs,[1 2],'rows');  p(ind) = ra.MC.hsd.CT_by_TI.pValue(1); 
     ind = ismember(combs,[1 3],'rows');  p(ind) = ra.MC.hsd.CT_by_TI.pValue(2); 
    ind = ismember(combs,[2 3],'rows');  p(ind) = ra.MC.hsd.CT_by_TI.pValue(4); 
    ind = ismember(combs,[4 5],'rows'); p(ind) = ra.MC.hsd.CT_by_TI.pValue(7); 
    ind = ismember(combs,[4 6],'rows'); p(ind) = ra.MC.hsd.CT_by_TI.pValue(8); 
    ind = ismember(combs,[5 6],'rows'); p(ind) = ra.MC.hsd.CT_by_TI.pValue(10); 
    h = p<0.05;
    xdata = make_xdata([3 3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = repmat(mData.dcolors(1:3),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    make_bars_hollow(hbs(4:end))
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%)',round(mra),round(semra)),'FontSize',6);
%     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dur','Dis','Mix'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.1 0.02]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_TI_CT.pdf'),600);
    %%


    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1.5 1 1]);
    xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(1:3);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',4,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dur','Dis','Mix'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.4 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_CT.pdf'),600);
    %%
    break;
end

%% check the difference in zMI for Dis and Dur for the different cell types DurC, DisC, and DDM
while 1
    %%
    mean_dzMI = [];
    FD_Prop = dzMI_FD.diff_T_D; FT_Prop = dzMI_FT.diff_T_D; 
%     FD_Prop = dzMI_FD.rs.diff_T_D; FT_Prop = dzMI_FT.rs.diff_T_D; 
%     FD_Prop = dzMI_FD.HaFD.diff_T_D; FT_Prop = dzMI_FT.HaFD.diff_T_D; 
%     FD_Prop = dzMI_FD.HiFD.diff_T_D; FT_Prop = dzMI_FT.HiFD.diff_T_D; 
    for rfi = 1
        TD = FD_Prop;
        cell_resp = FD_Dur_comp{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp)];
        cell_resp = FD_Dis_comp{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp)];
        cell_resp = FD_conj{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp)];

        TD = FT_Prop;
        cell_resp = FT_Dur_comp{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp)];
        cell_resp = FT_Dis_comp{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp)];
        cell_resp = FT_conj{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp)];
    end
    
    [within,dvn,xlabels] = make_within_table({'TI','CT'},[2,3]);
    dataT = make_between_table({mean_dzMI(2:4,:)},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova
    %%
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1.5 1 1]);
    xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(7:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1)-0.5 maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dur','Dis','Mix'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.4 0]); put_axes_labels(gca,{[],[0 0 0]},{'Difference',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('dzMI_CT.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI','hsd'},[1.5 1 1]);
    xdata = make_xdata([2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(10:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.04,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.025);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'T','I'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.5 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]})
    %%
    break;
end

%% Overlap Indices ImageSC
while 1
    
    rfi = 1;
    respAll = [FD_Dur_comp{rfi} FD_Dis_comp{rfi} FD_conj{rfi} FT_Dur_comp{rfi} FT_Dis_comp{rfi} FT_conj{rfi}];
%     respAll = [FD_Dur_comp{rfi} FD_Dis_comp{rfi} FD_conj{rfi}];
%     respAll = [FT_Dur_comp{rfi} FT_Dis_comp{rfi} FT_conj{rfi}];
    [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(respAll,0.5,0.05);
    mOI = mCI; semOI = semCI;
%     mOI = mean(uni,3); semOI = std(uni,[],3)/sqrt(5);
    sz = size(mOI,1);
%     mOI = OI_mat(:,:,4);
    oM = ones(size(mOI));
    
    mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
    maxI = max([mOI(:)]);%;semOI(:)]);    
    minI = min([mOI(:)]);%semOI(:)]);
    
    mask = tril(NaN(size(mOI)),0); mask(mask==0) = 1; 
    txl = {'T-Dur','T-Dis','T-Mix','I-Dur','I-Dis','I-Mix'};
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
    for rr = 1:size(mCI,1)
        for cc = 1:size(mCI,1)
            if rr == cc
                text(rr-0.25,cc,sprintf('%.0f',mCI(rr,cc)),'FontSize',5);
            end
        end
    end
    set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
%     set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
    set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',txl,'yticklabels',txl,'Ydir','reverse'); xtickangle(45);
%     text(-0.5,sz+1.1,'Average Overlap Index (N = 5 animals)','FontSize',5); 
    set(gca,'Ydir','normal');ytickangle(20);
    box on
    changePosition(gca,[-0.01 0 -0.04 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.05 0.09 0.05]);
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_%d_mean.pdf',rfi),600);
    %%
    break;
end
%% agglomerative hierarchical clustering
while 1
    %%
    mOI1 = mOI;
%     mOI1(isnan(mOI1)) = 1;
    Di = pdist(mOI1,@naneucdist);
%     tree = linkage(mOI1,'average','euclidean');
    tree = linkage(Di,'average');
    figure(hf);clf
    [H,T,TC] = dendrogram(tree,'Orientation','top','ColorThreshold','default');
    hf = gcf;
    set(hf,'Position',[7 3 3.5 1.5]);
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
    Rs = [RsDt(:,1) RsTi(:,1)];
    for ii = 1:length(ei)
        b1 = ei{ii}.b;
        for jj = 1:10
            alds(ii,jj) = b1.dist(b1.stim_r(jj+10)) - b1.dist(b1.air_puff_r(jj+20));
        end
    end
    ald = round(mean(alds(:)));
     ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 2],'spaceRowsCols',[0 0.01],'rightUpShifts',[0.03 0.13],'widthHeightAdjustment',[-60 -350]);
    set(gcf,'color','w'); set(gcf,'Position',[5 5 2 1]);
    [Y,E] = discretize(1:49,3);
    all_speeds = []; cTxt = {'3-Trials','3-Intertrials','4-Trials','4-Intertrials','5-Trials','5-Intertrials'}; 
    for cn = 1:2
        mean_speed_over_trials = [];
        aThisSpeed = [];
        for an = 1:size(Rs,1)
            thisSpeed = nanmean(Rs{an,cn}.speed);
            mean_speed_over_trials(an,:) = thisSpeed;
            for bb = 1:3
                aThisSpeed(an,bb) = nanmean(thisSpeed(Y==bb));
            end
        end
        all_speeds = [all_speeds aThisSpeed];
        axes(ff.h_axes(1,cn));
        hold on;
        xs = Rs{1,cn}.xs; N = length(xs);
        mspeed = mean(mean_speed_over_trials(:,1:N)); semspeed = std(mean_speed_over_trials(:,1:N))/sqrt(5);
        plot(xs,mspeed);
        shadedErrorBar(xs,mspeed,semspeed);
        changePosition(gca,[0.1 0.15 -0.05 -0.15]);
        if cn == 1
            put_axes_labels(gca,{'',[0 0 0]},{{'Speed (cm/sec)'},[0 -0.5 0]});
        end
        if ismember(cn,[1 3 5])
            if cn == 3
                plot([ald ald],[0 30],'m:','linewidth',0.5);
            end
%             plot([50 50],[0 30],'k--','linewidth',0.25);
%             plot([100 100],[0 30],'k--','linewidth',0.25);
%             bTxt = {'dB1','dB2','dB3'}; 
%             xbTxt = [25 75 125]-7; ybTxt = 31;
%             for ii = 1:length(bTxt)
%                 text(xbTxt(ii),ybTxt,bTxt{ii},'FontSize',5);
%             end
            xlabel('Distance (cm)');
        else
%             plot([5 5],[0 30],'k--','linewidth',0.25);
%             plot([10 10],[0 30],'k--','linewidth',0.25);
%             bTxt = {'tB1','tB2','tB3'}; xbTxt = [2.5 7.5 12.5]-1; ybTxt = 31;
%             for ii = 1:length(bTxt)
%                 text(xbTxt(ii),ybTxt,bTxt{ii},'FontSize',5);
%             end
            set(gca,'Ytick',[]);
            xlabel('Time (sec)');
        end
        text(xbTxt(1),ybTxt+5,cTxt{cn},'FontSize',5);
        ylim([0 30]);
        box off;
        format_axes(gca);
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('speeds_345'),600);

    %%
    all_speeds_Trials = all_speeds(:,[1:3]+0);
    [within,dvn,xlabels] = make_within_table({'Bins'},[3]);
    dataT = make_between_table({all_speeds_Trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Bins','hsd'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors(7:9);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 37]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'yticklabels',[]); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.35 -0.05])
%     put_axes_labels(gca,{[],[0 0 0]},{'Avg. Speed (cm/sec)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('avg_speed_anova_pooled.pdf'),600);
    %%
    all_speeds_Trials = all_speeds(:,[1:3 7:9 13:15]+3);
    [within,dvn,xlabels] = make_within_table({'Conds','Bins'},[3,3]);
    dataT = make_between_table({all_speeds_Trials},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Bins','hsd'},[1 1 1]);
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.colors([10 11 12]);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    ylims = ylim;
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 37]);
    format_axes(gca);
    xticks = xdata; xticklabels = {'tB1','tB2','tB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.02 -0.35 -0.05])
    put_axes_labels(gca,{[],[0 0 0]},{'Avg. Speed (cm/sec)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('avg_speed_anova_pooled_IT.pdf'),600);
    %%
    break;
end