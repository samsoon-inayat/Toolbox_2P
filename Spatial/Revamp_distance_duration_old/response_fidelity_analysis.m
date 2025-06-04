function response_fidelity_analysis

%%
RsDt = o.Rs(:,[Ar_t_D ArL_t_D Ars_t_D]);  RsTt = o.Rs(:,[Ar_t_T ArL_t_T Ars_t_T]);
RsDi = o.Rs(:,[Ar_i_D ArL_i_D Ars_i_D]);  RsTi = o.Rs(:,[Ar_i_T ArL_i_T Ars_i_T]);

%%
respfids = {[30 40],[50 60],[70 80],[90 100]};
respfids = {[30 60],[70 100]};

raster_types = {'RsTt','RsDt','RsTi','RsDi'};
raster_types = {'RsTt','RsTi'};

clear props
for ii = 1:length(respfids)
    trialsR = respfids{ii};
    for jj = 1:length(raster_types)
        cmdTxt = sprintf('props{ii,jj} = get_props_Rs(%s,trialsR);',raster_types{jj});
        eval(cmdTxt);
    end
end

%%
ptcl = [];
for ii = 1:length(respfids)
    for jj = 1:length(raster_types)
        tcl = props{ii,jj}.good_FR;
        ptcl = [ptcl 100*exec_fun_on_cell_mat(tcl,'sum')./exec_fun_on_cell_mat(tcl,'length')];
    end
end

[within,dvn,xlabels] = make_within_table({'RF','TI','Cond'},[2,2,3]);
dataT = make_between_table({ptcl},dvn);
ra = RMA(dataT,within,{'hsd'});
ra.ranova


%%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'RF_TI_Cond','hsd'},[1.5 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3 3],[1 1.5]);
    hf = get_figure(5,[8 7 3.5 1]);
    tcolors = repmat(mData.colors(1:6),1,4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.15);
    maxY = maxY + 0;
    ylims = ylim;
    format_axes(gca);
%     htxt = text(0.75,maxY-3,sprintf('Any Condition or Type (%d\x00B1%d%%)',round(mra),round(semra)),'FontSize',6);
%     text(13,maxY-9,sprintf('%d\x00B1%d %%',round(mrall),round(semrall)),'FontSize',6);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'3','4','5'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[-0.04 0.01 0.1 0]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('resp_fid_2.pdf'),600);

    %%
     %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels] = get_vals_for_bar_graph_RMA(mData,ra,{'RF_by_TI','hsd'},[1.5 1 1]);
     h(h==1) = 0;
%      p = ones(size(p)); 
%      ind = ismember(combs,[1 2],'rows');  p(ind) = ra.MC.hsd.CT_by_TI.pValue(1); 
%      ind = ismember(combs,[1 3],'rows');  p(ind) = ra.MC.hsd.CT_by_TI.pValue(2); 
%     ind = ismember(combs,[2 3],'rows');  p(ind) = ra.MC.hsd.CT_by_TI.pValue(4); 
%     ind = ismember(combs,[4 5],'rows'); p(ind) = ra.MC.hsd.CT_by_TI.pValue(7); 
%     ind = ismember(combs,[4 6],'rows'); p(ind) = ra.MC.hsd.CT_by_TI.pValue(8); 
%     ind = ismember(combs,[5 6],'rows'); p(ind) = ra.MC.hsd.CT_by_TI.pValue(10); 
%     h = p<0.05;
    xdata = make_xdata([2 2 2 2],[1 1.5]);
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
    xticks = xdata; xticklabels = {'T','I'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[0.03 0.01 -0.1 0.02]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('resp_fid_2_ph.pdf'),600);

