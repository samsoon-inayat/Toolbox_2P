function dist_dur

%%
RsDt = o.Rs(:,[Ar_t_D ArL_t_D Ars_t_D]);  RsTt = o.Rs(:,[Ar_t_T ArL_t_T Ars_t_T]);
RsDi = o.Rs(:,[Ar_i_D ArL_i_D Ars_i_D]);  RsTi = o.Rs(:,[Ar_i_T ArL_i_T Ars_i_T]);
[dzMI_FD,dzMI_FT] = get_zMI_comp_dist_time(RsDt,RsTt,RsDi,RsTi);

%%
% respfids = {[30 40],[50 60],[70 80],[90 100]};
respfids = {[30 60],[70 100],[30 100],[50 100]};

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
%% find dis, dur, and mix cells
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

%% find dis, dur, and mix cells --> pooled across conditions because condition had no main effect regardless of response fidelity
while 1
    for ii = 1:length(respfids)
        for jj = 1:length(raster_types)
            cell_list{ii,jj} = props{ii,jj}.good_FR_and_tuned;
        end
    end

    clear FD_Dis_comp FD_Dur_comp FD_conj FT_Dis_comp FT_Dur_comp FT_conj
    for ii = 1:length(respfids)
        FD_Dur = cell_list_op(cell_list_op(props{ii,1}.vals,props{ii,1}.good_FR,'and'),[],'or',1); FD_Dis = cell_list_op(cell_list_op(props{ii,2}.vals,props{ii,1}.good_FR,'and'),[],'or',1);
        FT_Dur = cell_list_op(cell_list_op(props{ii,3}.vals,props{ii,3}.good_FR,'and'),[],'or',1); FT_Dis = cell_list_op(cell_list_op(props{ii,4}.vals,props{ii,3}.good_FR,'and'),[],'or',1);

        for rr = 1:size(FD_Dis,1)
            for cc = 1:size(FD_Dis,2)
                cellP1 = FD_Dis(rr,cc); cellP2 = FD_Dur(rr,cc);
                FD_Dis_comp{ii}(rr,cc) = cell_list_op(cellP1,cell_list_op(cellP2,[],'not'),'and');
                FD_Dur_comp{ii}(rr,cc) = cell_list_op(cell_list_op(cellP1,[],'not'),cellP2,'and');
                FD_conj{ii}(rr,cc) = cell_list_op(cellP1,cellP2,'and');

                cellP1 = FT_Dis(rr,cc); cellP2 = FT_Dur(rr,cc);
                FT_Dis_comp{ii}(rr,cc) = cell_list_op(cellP1,cell_list_op(cellP2,[],'not'),'and');
                FT_Dur_comp{ii}(rr,cc) = cell_list_op(cell_list_op(cellP1,[],'not'),cellP2,'and');
                FT_conj{ii}(rr,cc) = cell_list_op(cellP1,cellP2,'and');
            end
        end
    end
break;
end
%% choosing 1 response fidelity and comparing across cell types, conditions, trials and intertrials
while 1
rfi = 4;

resp = [FD_Dur_comp{rfi} FD_Dis_comp{rfi} FD_conj{rfi} FT_Dur_comp{rfi} FT_Dis_comp{rfi} FT_conj{rfi}];

per_resp = 100*exec_fun_on_cell_mat(resp,'sum')./exec_fun_on_cell_mat(resp,'length');

[within,dvn,xlabels] = make_within_table({'TI','CT','Cond'},[2,3,3]);
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
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_CT_Cond','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3 3 3 3],[1 1.5]);
    hf = get_figure(5,[8 7 3.5 1]);
    tcolors = repmat(mData.colors(1:9),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[-0.04 0.01 0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_all.pdf'),600);
  %%
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
    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
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
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI','hsd'},[1.5 1 1]);
    xdata = make_xdata([2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(10:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',4,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 25]); format_axes(gca);
    xticks = xdata; xticklabels = {'T','I'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.5 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]})
    save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_TI.pdf'),600);
    
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1.5 1 1]);
    xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(7:end);
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
    break
end
%% %% choosing 1 condition and comparing across cell types, response fidelity, trials and intertrials
while 1
    %%
    cni = 1;
    rfi = 4;
    resp = [FD_Dur_comp{rfi}(:,cni) FD_Dis_comp{rfi}(:,cni) FD_conj{rfi}(:,cni) FT_Dur_comp{rfi}(:,cni) FT_Dis_comp{rfi}(:,cni) FT_conj{rfi}(:,cni)];
%     rfi = 2;
%     resp = [resp FD_Dur_comp{rfi}(:,cni) FD_Dis_comp{rfi}(:,cni) FD_conj{rfi}(:,cni) FT_Dur_comp{rfi}(:,cni) FT_Dis_comp{rfi}(:,cni) FT_conj{rfi}(:,cni)];

    per_resp = 100*exec_fun_on_cell_mat(resp,'sum')./exec_fun_on_cell_mat(resp,'length');

    [within,dvn,xlabels] = make_within_table({'RF','TI','CT'},[2,2,3]);
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
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'RF_TI_CT','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3 3],[1 1.5]);
    hf = get_figure(5,[8 7 3.5 1]);
    tcolors = repmat(mData.colors(1:9),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dur','Dis','Mix'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[-0.04 0.01 0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_all.pdf'),600);
    %%
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'RF_by_TI','hsd'},[1.5 1 1]);
     h(h==1) = 0;
%      p = ones(size(p)); 
%      ind = ismember(combs,[1 2],'rows');  p(ind) = ra.MC.hsd.RF_by_CT.pValue(1); 
%      ind = ismember(combs,[1 3],'rows');  p(ind) = ra.MC.hsd.CT_by_TI.pValue(2); 
%     ind = ismember(combs,[2 3],'rows');  p(ind) = ra.MC.hsd.CT_by_TI.pValue(4); 
%     ind = ismember(combs,[4 5],'rows'); p(ind) = ra.MC.hsd.CT_by_TI.pValue(7); 
%     ind = ismember(combs,[4 6],'rows'); p(ind) = ra.MC.hsd.CT_by_TI.pValue(8); 
%     ind = ismember(combs,[5 6],'rows'); p(ind) = ra.MC.hsd.CT_by_TI.pValue(10); 
%     h = p<0.05;
    xdata = make_xdata([2 2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%)',round(mra),round(semra)),'FontSize',6);
%     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'RF1','RF2'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.1 0.02]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_TI_CT.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_CT','hsd'},[1.5 1 1]);
     h(h==1) = 0;
%      p = ones(size(p)); 
%      ind = ismember(combs,[1 2],'rows');  p(ind) = ra.MC.hsd.RF_by_CT.pValue(1); 
%      ind = ismember(combs,[1 3],'rows');  p(ind) = ra.MC.hsd.CT_by_TI.pValue(2); 
%     ind = ismember(combs,[2 3],'rows');  p(ind) = ra.MC.hsd.CT_by_TI.pValue(4); 
%     ind = ismember(combs,[4 5],'rows'); p(ind) = ra.MC.hsd.CT_by_TI.pValue(7); 
%     ind = ismember(combs,[4 6],'rows'); p(ind) = ra.MC.hsd.CT_by_TI.pValue(8); 
%     ind = ismember(combs,[5 6],'rows'); p(ind) = ra.MC.hsd.CT_by_TI.pValue(10); 
%     h = p<0.05;
    xdata = make_xdata([3 3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%)',round(mra),round(semra)),'FontSize',6);
%     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dur','Dis','Mix'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.1 0.02]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    %%
    break;
end

%% %% choosing 1 cell type and comparing across response fidelity, conditions, trials and intertrials
while 1
    %%
    rfi = 1;
    resp = [FD_Dur_comp{rfi} FT_Dur_comp{rfi}];
    rfi = 2;
    resp = [resp FD_Dur_comp{rfi} FT_Dur_comp{rfi}];
    
%     rfi = 1;
%     resp = [FD_Dis_comp{rfi} FT_Dis_comp{rfi}];
%     rfi = 2;
%     resp = [resp FD_Dis_comp{rfi} FT_Dis_comp{rfi}];
% 
    rfi = 1;
    resp = [FD_conj{rfi} FT_conj{rfi}];
    rfi = 2;
    resp = [resp FD_conj{rfi} FT_conj{rfi}];
    
    per_resp = 100*exec_fun_on_cell_mat(resp,'sum')./exec_fun_on_cell_mat(resp,'length');

    [within,dvn,xlabels] = make_within_table({'RF','TI','Cond'},[2,2,3]);
    dataT = make_between_table({per_resp},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova

    
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'RF_TI_Cond','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3 3],[1 1.5]);
    hf = get_figure(7,[8 7 3.5 1]);
    tcolors = repmat(mData.colors(1:6),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[-0.04 0.01 0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_all.pdf'),600);
    %%
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'RF_by_CT','hsd'},[1.5 1 1]);
     h(h==1) = 0;
%      p = ones(size(p)); 
%      ind = ismember(combs,[1 2],'rows');  p(ind) = ra.MC.hsd.RF_by_CT.pValue(1); 
%      ind = ismember(combs,[1 3],'rows');  p(ind) = ra.MC.hsd.CT_by_TI.pValue(2); 
%     ind = ismember(combs,[2 3],'rows');  p(ind) = ra.MC.hsd.CT_by_TI.pValue(4); 
%     ind = ismember(combs,[4 5],'rows'); p(ind) = ra.MC.hsd.CT_by_TI.pValue(7); 
%     ind = ismember(combs,[4 6],'rows'); p(ind) = ra.MC.hsd.CT_by_TI.pValue(8); 
%     ind = ismember(combs,[5 6],'rows'); p(ind) = ra.MC.hsd.CT_by_TI.pValue(10); 
%     h = p<0.05;
    xdata = make_xdata([2 2 2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0.75,maxY-3,sprintf('Any Condition (%d\x00B1%d%%)',round(mra),round(semra)),'FontSize',6);
%     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'RF1','RF2'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.1 0.02]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_TI_CT.pdf'),600);
    %%
    break;
end

%%  4-within factors choosing all response fidelity and comparing across cell types, conditions, trials and intertrials
while 1
rfi = 1;
resp = [FD_Dur_comp{rfi} FD_Dis_comp{rfi} FD_conj{rfi} FT_Dur_comp{rfi} FT_Dis_comp{rfi} FT_conj{rfi}];

rfi = 2;
resp = [resp FD_Dur_comp{rfi} FD_Dis_comp{rfi} FD_conj{rfi} FT_Dur_comp{rfi} FT_Dis_comp{rfi} FT_conj{rfi}];

per_resp = 100*exec_fun_on_cell_mat(resp,'sum')./exec_fun_on_cell_mat(resp,'length');

[within,dvn,xlabels] = make_within_table({'RF','TI','CT','Cond'},[2,2,3,3]);
dataT = make_between_table({per_resp},dvn);
ra = RMA(dataT,within,{'hsd'});
ra.ranova

%%
[xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'RF_by_TI_by_CT','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3],[1 1.5]);
    hf = get_figure(5,[8 7 3.5 1]);
    tcolors = repmat(mData.colors(1:9),1,2);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 25]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[-0.04 0.01 0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_all.pdf'),600);
  %%
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
    tcolors = mData.dcolors;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
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
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'RF','hsd'},[1.5 1 1]);
    xdata = make_xdata([2],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(10:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',4,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 5;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 25]); format_axes(gca);
    xticks = xdata; xticklabels = {'T','I'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.5 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]})
    save_pdf(hf,mData.pdf_folder,sprintf('perc_cells_TI.pdf'),600);
    
    
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1.5 1 1]);
    xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(7:end);
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
    break
end

%% check the difference in zMI for Dis and Dur for the different cell types DurC, DisC, and DDM
while 1
    %%
    mean_dzMI = [];
    FD_Prop = dzMI_FD.diff_T_D; FT_Prop = dzMI_FT.diff_T_D; 
    FD_Prop = dzMI_FD.rs.diff_T_D; FT_Prop = dzMI_FT.rs.diff_T_D; 
%     FD_Prop = dzMI_FD.HaFD.diff_T_D; FT_Prop = dzMI_FT.HaFD.diff_T_D; 
%     FD_Prop = dzMI_FD.HiFD.diff_T_D; FT_Prop = dzMI_FT.HiFD.diff_T_D; 
    for rfi = 1:2
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
    
    [within,dvn,xlabels] = make_within_table({'RF','TI','CT','Cond'},[2,2,3,3]);
    dataT = make_between_table({mean_dzMI},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova
    %%
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1.5 1 1]);
    xdata = make_xdata([3],[1 1.5]);
    hf = get_figure(5,[8 7 1.5 1]);
    tcolors = mData.dcolors(7:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.015,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dur','Dis','Mix'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.4 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
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


%% for only one condition e.g., 3 check the difference in zMI for Dis and Dur for the different cell types DurC, DisC, and DDM
while 1
    %%
    cni = 1;
    mean_dzMI = [];
    FD_Prop = dzMI_FD.diff_T_D; FT_Prop = dzMI_FT.diff_T_D; 
%     FD_Prop = dzMI_FD.rs.diff_T_D; FT_Prop = dzMI_FT.rs.diff_T_D; 
%     FD_Prop = dzMI_FD.HaFD.diff_T_D; FT_Prop = dzMI_FT.HaFD.diff_T_D; 
%     FD_Prop = dzMI_FD.HiFD.diff_T_D; FT_Prop = dzMI_FT.HiFD.diff_T_D; 
    for rfi = 1:2
        TD = FD_Prop(:,cni);
        cell_resp = FD_Dur_comp{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp(:,cni))];
        cell_resp = FD_Dis_comp{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp(:,cni))];
        cell_resp = FD_conj{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp(:,cni))];

        TD = FT_Prop(:,cni);
        cell_resp = FT_Dur_comp{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp(:,cni)) ];
        cell_resp = FT_Dis_comp{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp(:,cni))];
        cell_resp = FT_conj{rfi};
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',cell_resp(:,cni))];
    end
    
%     [within,dvn,xlabels] = make_within_table({'TI','CT'},[2,3]);
%     dataT = make_between_table({mean_dzMI},dvn);
%     ra = RMA(dataT,within,{'hsd'});
%     ra.ranova
    
    [within,dvn,xlabels] = make_within_table({'RF','TI','CT'},[2,2,3]);
    dataT = make_between_table({mean_dzMI},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'RF_TI_CT','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3 3],[1 1.5]);
    hf = get_figure(5,[8 7 3.5 1]);
    tcolors = mData.dcolors(7:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.015,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dur','Dis','Mix'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('dzMI_CT.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'TI_by_CT','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3],[1 1.5]);
    hf = get_figure(5,[8 7 3.5 1]);
    tcolors = mData.dcolors(7:end);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.015,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY;
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Dur','Dis','Mix'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('dzMI_CT.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1.5 1 1]);
    xdata = make_xdata([3],[1 1.5]);
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
    rfi = 2;
    respAll = [FD_Dur_comp{rfi} FD_Dis_comp{rfi} FD_conj{rfi} FT_Dur_comp{rfi} FT_Dis_comp{rfi} FT_conj{rfi}];
    respAll = [FD_Dur_comp{rfi} FD_Dis_comp{rfi} FD_conj{rfi}];
    respAll = [FT_Dur_comp{rfi} FT_Dis_comp{rfi} FT_conj{rfi}];
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
    txl = {'3-T-Dur','4-T-Dur','5-T-Dur','3-T-Dis','4-T-Dis','5-T-Dis','3-T-Mix','4-T-Mix','5-T-Mix',...
        '3-I-Dur','4-I-Dur','5-I-Dur','3-I-Dis','4-I-Dis','5-I-Dis','3-I-Mix','4-I-Mix','5-I-Mix'};
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
