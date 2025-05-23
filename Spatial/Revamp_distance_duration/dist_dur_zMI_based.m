function dist_dur

%%
RsDt = o.Rs(:,[Ar_t_D ArL_t_D Ars_t_D]);  RsTt = o.Rs(:,[Ar_t_T ArL_t_T Ars_t_T]);
RsDi = o.Rs(:,[Ar_i_D ArL_i_D Ars_i_D]);  RsTi = o.Rs(:,[Ar_i_T ArL_i_T Ars_i_T]);

[dzMI_FD,dzMI_FT] = get_zMI_comp_dist_time(RsDt,RsTt,RsDi,RsTi);
%%
% respfids = {[30 40],[50 60],[70 80],[90 100]};
respfids = {[30 60],[70 100]};

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
            cell_list{ii,jj} = cell_list_op(props{ii,jj}.good_FR_and_tuned,props{ii,jj}.good_zMI,'and');
        end
    end

    clear FD_Dis_comp FD_Dur_comp FD_conj FT_Dis_comp FT_Dur_comp FT_conj
    for ii = 1:length(respfids)
        FD_Dur = cell_list_op(props{ii,1}.vals,props{ii,1}.good_FR,'and'); FD_Dis = cell_list_op(props{ii,2}.vals,props{ii,1}.good_FR,'and');
        FT_Dur = cell_list_op(props{ii,3}.vals,props{ii,3}.good_FR,'and'); FT_Dis = cell_list_op(props{ii,4}.vals,props{ii,3}.good_FR,'and');

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
rfi = 2;

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
