function sensory_vs_spatial

ntrials = [70,100];
si_nb = [Ar_t_D ArL_t_D Ars_t_D];
% si_nb = [Ar_i_T ArL_i_T Ars_i_T];
props_nb = get_props_Rs(o.Rs(:,si_nb),ntrials);
resp_nb = props_nb.good_FR;
si = [Ab_T Abs_T];
Rs = o.Rs(:,si);mR = o.mR(:,si);
props1 = get_props_Rs(Rs,ntrials);
exc = props1.good_FR_and_exc;     sup = props1.good_FR_and_inh;    com = props1.good_FR_and_untuned;

respE_OR = cell_list_op(exc,[],'or');
respS_OR = cell_list_op(sup,[],'or');
respC_OR = cell_list_op(com,[],'or'); 

resp = [respE_OR(:,1) respS_OR(:,1) respC_OR(:,1) resp_nb];

disp('done');

clear valsE valsS valsC valsO
allzMIs = props_nb.zMI;
respE = repmat(exc(:,1),1,size(allzMIs,2)); 
respS = repmat(sup(:,1),1,size(allzMIs,2));
respC = repmat(com(:,1),1,size(allzMIs,2)); 
respO = repmat(cell_list_op(props1.good_FR(:,1),[],'not'),1,size(allzMIs,2));

for rr = 1:size(allzMIs,1)
    for cc = 1:size(allzMIs,2)
        valsE(rr,cc) = nanmean(allzMIs{rr,cc}(respE{rr,cc}));
        valsS(rr,cc) = nanmean(allzMIs{rr,cc}(respS{rr,cc}));
        valsC(rr,cc) = nanmean(allzMIs{rr,cc}(respC{rr,cc}));
        valsO(rr,cc) = nanmean(allzMIs{rr,cc}(respO{rr,cc}));
    end
end
disp('done');
%%
while 1
    %% For zMIs of cells over belt
    [within,dvn,xlabels] = make_within_table({'CT','Cond'},[4,3]);
    dataT = make_between_table({valsE,valsS,valsC,valsO},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT_Cond','hsd'},[1 1 1]);
    h(h==1) = 0;
    xdata = make_xdata([3 3 3 3],[1 2]); xdata(10:end) = xdata(10:end) + 1;
    hf = get_figure(5,[8 7 3 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = repmat(mData.colors(1:9),2,1);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) maxY]); format_axes(gca);
    xticks = xdata; xticklabels = {'Ar','ArL','Ar*'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[-0.041 0.0 0.12 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{{'zMI'},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_cells_E_S_C.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'CT','hsd'},[1 1 1]);
    xdata = make_xdata([4],[1 2]); 
    hf = get_figure(5,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(1:4);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.85,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
%     make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 6]); format_axes(gca);
    xticks = xdata; xticklabels = {'Exc','Sup','Com','Oth'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.5 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_cells_Sen_Sp.pdf'),600);
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Bin','hsd'},[1 1 1]);
    xdata = make_xdata([3],[1 2]); 
    hf = get_figure(6,[8 7 1.25 1]);
    % s = generate_shades(length(bins)-1);
    tcolors = mData.dcolors(3:5);
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',8,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     maxY = maxY;
    make_bars_hollow(hbs(10:end));
    ylims = ylim;
    format_axes(gca);
    set_axes_limits(gca,[0.35 xdata(end)+.65],[ylims(1) 6]); format_axes(gca);
    xticks = xdata; xticklabels = {'dB1','dB2','dB3'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(45)
    changePosition(gca,[0.05 0.0 -0.4 -0.09]); put_axes_labels(gca,{[],[0 0 0]},{{''},[0 0 0]});
%     ht = title('Across Lb and Lb*'); changePosition(ht,[-1 0 0]);
    save_pdf(hf,mData.pdf_folder,sprintf('zMI_cells_on_belt_pooled_Bin.pdf'),600);
    %%
    break;
end

