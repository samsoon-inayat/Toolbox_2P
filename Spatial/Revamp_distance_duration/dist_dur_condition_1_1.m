function dist_dur

%%
clear dis_cells_T dur_cells_T dis_cells_I dur_cells_I
rfi = 1;
for ani = 1:5
    thisan = {dzMI_FD.diff_D_T{ani} > 0};
    thisan = cell_list_op(FD_conj{rfi}(ani),thisan,'and');
    dis_cells_T(ani,1) = cell_list_op(FD_Dis_comp{rfi}(ani),thisan,'or');
    thisan = {dzMI_FD.diff_D_T{ani} < 0};
    thisan = cell_list_op(FD_conj{rfi}(ani),thisan,'and');
    dur_cells_T(ani,1) = cell_list_op(FD_Dur_comp{rfi}(ani),thisan,'or');
    
    thisan = {dzMI_FT.diff_T_D{ani} > 0};
    thisan = cell_list_op(FT_conj{rfi}(ani),thisan,'and');
    dur_cells_I(ani,1) = cell_list_op(FT_Dur_comp{rfi}(ani),thisan,'or');
    thisan = {dzMI_FT.diff_T_D{ani} < 0};
    thisan = cell_list_op(FT_conj{rfi}(ani),thisan,'and');
    dis_cells_I(ani,1) = cell_list_op(FT_Dis_comp{rfi}(ani),thisan,'or');
end

resp = [dis_cells_T dur_cells_T dis_cells_I dur_cells_I];

%% percentage of cells
while 1
per_cells = find_percent(resp);
[within,dvn,xlabels] = make_within_table({'TI','CT'},[2,3]);
dataT = make_between_table({per_cells},dvn);
ra = RMA(dataT,within,{'hsd'});
ra.ranova

[mr,semr] = findMeanAndStandardError(per_cells);

any_cells = cell_list_op(resp,[],'or',1);
per_active_any = 100*exec_fun_on_cell_mat(any_cells,'sum')./exec_fun_on_cell_mat(any_cells,'length');
any_cells = cell_list_op(resp,[],'and',1);
per_active_all = 100*exec_fun_on_cell_mat(any_cells,'sum')./exec_fun_on_cell_mat(any_cells,'length');

[mra,semra] = findMeanAndStandardError(per_active_any);
[mrall,semrall] = findMeanAndStandardError(per_active_all);

%%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_CT','hsd'},[1.5 1 1]);
    xdata = make_xdata([3 3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = repmat(mData.colors(1:3),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1)-0.5 maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dis','Dur'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.4 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('perc_cells.pdf'),600);
    %%
    break;
end
%% check the difference in zMI for Dis and Dur for the different cell types DurC, DisC, and DDM
while 1
    %%
    mean_dzMI = [];
    FD_Prop = dzMI_FD.diff_D_T; FT_Prop = dzMI_FT.diff_D_T; 
%     FD_Prop = dzMI_FD.rs.diff_T_D; FT_Prop = dzMI_FT.rs.diff_T_D; 
%     FD_Prop = dzMI_FD.HaFD.diff_T_D; FT_Prop = dzMI_FT.HaFD.diff_T_D; 
%     FD_Prop = dzMI_FD.HiFD.diff_T_D; FT_Prop = dzMI_FT.HiFD.diff_T_D; 
        TD = FD_Prop;
        cell_resp = dis_cells_T;
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp)];
        cell_resp = dur_cells_T;
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp)];
        TD = FT_Prop;
        cell_resp = dis_cells_I;
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp)];
        cell_resp = dur_cells_I;
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp)];
   
    [within,dvn,xlabels] = make_within_table({'TI','CT'},[2,2]);
    dataT = make_between_table({mean_dzMI},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova
    %%
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_CT','hsd'},[1.5 1 1]);
    xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(7:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1)-0.5 maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dis','Dur'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.4 0]); put_axes_labels(gca,{[],[0 0 0]},{'Difference',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('dzMI_CT.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1.5 1 1]);
    xdata = make_xdata([2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(10:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.025);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1)-0.25 maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dis','Dur'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.5 0]); put_axes_labels(gca,{[],[0 0 0]},{'Difference zMI',[0 0 0]})
    %%
    break;
end
%% Overlap Indices ImageSC
while 1
    si = [Ab_On Abs_On]; props_ON = get_props_Rs(o.Rs(:,si)); resp_ON = cell_list_op(props_ON.vals,[],'or',1);
    si = [Ab_Off Abs_Off]; props_OFF = get_props_Rs(o.Rs(:,si)); resp_OFF = cell_list_op(props_OFF.vals,[],'or',1);
    respAll = [resp_ON resp resp_OFF];

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
    txl = {'A-On','Dis-T','Dur-T','Dis-I','Dur-I','A-Off'};
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
    changePosition(gca,[-0.01 0 -0.08 0]);
    hc = putColorBar(gca,[0.09 0.07 -0.11 -0.15],{sprintf('%.1f',minI),sprintf('%.1f',maxI)},6,'eastoutside',[0.1 0.07 0.09 0.09]);
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
    set(hf,'Position',[7 3 1 1]);
    set(H,'linewidth',1);
    set(gca,'xticklabels',txl(TC));xtickangle(45);
    format_axes(gca);
    hx = ylabel('Eucledian Distance');changePosition(hx,[0 -0.1 0]);
    changePosition(gca,[0.1 0.01 -0.01 0.0]);
    save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_cluster.pdf'),600);
    %%
    break;
end
%% Venn Diagram
while 1
    %%
    hf = get_figure(6,[10 7 1.5 1]);
    good_FRV = [dis_cells dur_cells];;
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
    set(HVenn(1),'FaceColor',tcolors{1},'FaceAlpha',0.75); set(HVenn(2),'FaceColor',tcolors{2},'FaceAlpha',0.75);
    ylims = ylim;
    save_pdf(hf,mData.pdf_folder,sprintf('conjunctive_cells_venn_diagram.pdf'),600);
    btxt = (sprintf('%0.1f%c%0.1f%%',m_p_gFR_C(1,1),pmchar,sem_p_gFR_C(1,1)));
    nbtxt = (sprintf('%0.1f%c%0.1f%%',m_p_gFR_C(2,2),pmchar,sem_p_gFR_C(2,2)));
    inttxt = (sprintf('%0.1f%c%0.1f%%',m_p_gFR_C(1,2),pmchar,sem_p_gFR_C(1,2)));
%     ht = title(sprintf('%s   %s   %s',btxt,inttxt,nbtxt));set(ht,'FontWeight','normal','FontSize',5);
    clc
    disp(btxt);
    disp(inttxt);
    disp(nbtxt);
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


%% distribution of percentage of cells on the belt OR intertrials
while 1
    %%
    propsD = get_props_Rs(RsDt,trialsR); propsT = get_props_Rs(RsTi,trialsR);
    
    bins = propsD.peak_location_bin;
    dis_cells_T_b = distribution_of_cells_ocross(dis_cells_T,bins);
    dur_cells_T_b = distribution_of_cells_ocross(dur_cells_T,bins);
    bins = propsT.peak_location_bin;
    dis_cells_I_b = distribution_of_cells_ocross(dis_cells_I,bins);
    dur_cells_I_b = distribution_of_cells_ocross(dur_cells_I,bins);
    
    resp = [dis_cells_T_b dur_cells_T_b dis_cells_I_b dur_cells_I_b];
    per_resp = find_percent(resp);
    [within,dvn,xlabels] = make_within_table({'TI','DT','Bins'},[2,2,2]);
    dataT = make_between_table({per_resp},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova
    %%
    propsDt = get_props_Rs(RsDt,trialsR); propsTi = get_props_Rs(RsTi,trialsR); propsDi = get_props_Rs(RsDi,trialsR); propsTi = get_props_Rs(RsTt,trialsR);

    FD_Prop = dzMI_FD.diff_D_T; FT_Prop = dzMI_FT.diff_D_T; 
    mean_dzMI = [];
    TD = repmat(FD_Prop,1,2);
    mean_dzMI = [mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dis_cells_T_b)];
%     mean_dzMI = [mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dur_cells_T_b)];
    TD = repmat(FT_Prop,1,2);
    mean_dzMI = [mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dis_cells_I_b)];
%     mean_dzMI = [mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dur_cells_I_b)];
    
    [within,dvn,xlabels] = make_within_table({'TI','Bins'},[2,2]);
    dataT = make_between_table({mean_dzMI},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova

    %%
    propsDt = get_props_Rs(RsDt,trialsR); propsTi = get_props_Rs(RsTi,trialsR); propsDi = get_props_Rs(RsDi,trialsR); propsTi = get_props_Rs(RsTt,trialsR);
    mean_dzMI = [];
    TD = repmat(propsDt.PWs,1,2);
    mean_dzMI = [mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dis_cells_T_b)];
%     mean_dzMI = [mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dur_cells_T_b)];
    TD = repmat(propsDi.PWs,1,2);
    mean_dzMI = [mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dis_cells_I_b)];
%     mean_dzMI = [mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dur_cells_I_b)];
    
    [within,dvn,xlabels] = make_within_table({'TI','Bins'},[2,2]);
    dataT = make_between_table({mean_dzMI},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova



    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_DT','hsd'},[1.5 1 1]);
    xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(10:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',1,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.025);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'B1','B2','B3','B4'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.5 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]})
    
    %%
    break;
end
