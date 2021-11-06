%% run conjunctive_representation file first to load data
selected_property = 'good_FR';
cmdTxt = sprintf('good_FR = props1.%s;',selected_property);

%% population vector and correlation sensory
while 1
    an = 4;
    titles = {'Ar','ArC','ArCB','ArB'};
    si = [Ar_t_D ArC_t_D ArCB_t_D ArB_t_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    eval(cmdTxt);
    ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[2 4],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.13],'widthHeightAdjustment',...
        [-50 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 3 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,good_FR,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
%     for ii = 1:length(ff.h_axes(1,:)) set_obj(ff.h_axes(2,ii),{'xtick',[1 18 38],'xticklabels',[-2 0 2.2]}); set_obj(ff.h_axes(2,ii),{'ytick',[1 18 38],'yticklabels',[-2 0 2.2]}); end 
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    changePosition(ff.h_axes(2,1).YLabel,[-2.5 0 0]); 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_light_%s.pdf',selected_property),600);

    % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 4],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.27],'widthHeightAdjustment',...
        [-50 -370]);    set(gcf,'color','w');    set(gcf,'Position',[10 8 2.5 0.85]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
%     for ii = 1:length(ff.h_axes(1,:)) set_obj(ff.h_axes(1,ii),{'xtick',[1 18 38],'xticklabels',[-2 0 2.2]}); set_obj(ff.h_axes(1,ii),{'ytick',[1 18 38],'yticklabels',[-2 0 2.2]}); end 
    changePosition(ff.h_axes(1,1).YLabel,[-2.5 0 0]); 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_light_%s.pdf',selected_property),600);
    %%
    break;
end

%% population vector and correlation sensory tuned
while 1
    an = 4;
    selected_property = 'tuned';
    titles = {'Lb','ArL-L','Lb*'};
    si = [Lb_T ArL_L_T Lbs_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    resp = props1.good_FR_and_untuned;
%     resp = props1.good_FR_and_inh;
    ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.13],'widthHeightAdjustment',...
        [-50 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.5 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) set_obj(ff.h_axes(2,ii),{'xtick',[1 18 38],'xticklabels',[-2 0 2.2]}); set_obj(ff.h_axes(2,ii),{'ytick',[1 18 38],'yticklabels',[-2 0 2.2]}); end 
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    changePosition(ff.h_axes(2,1).YLabel,[-2.5 0 0]); 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_light_%s.pdf',selected_property),600);

    % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.27],'widthHeightAdjustment',...
        [-50 -370]);    set(gcf,'color','w');    set(gcf,'Position',[10 8 2.5 0.85]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    for ii = 1:length(ff.h_axes(1,:)) set_obj(ff.h_axes(1,ii),{'xtick',[1 18 38],'xticklabels',[-2 0 2.2]}); set_obj(ff.h_axes(1,ii),{'ytick',[1 18 38],'yticklabels',[-2 0 2.2]}); end 
    changePosition(ff.h_axes(1,1).YLabel,[-2.5 0 0]); 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_light_%s.pdf',selected_property),600);
    %%
    break;
end

%% population vector and correlation sensory
while 1
    selected_property = 'tuned';
    an = 4;
    titles = {'Lb','Lb*'};
    si = [Lb_T Lbs_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    resp = props1.good_FR_and_tuned;
%     eval(cmdTxt);
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 2],...
        'spaceRowsCols',[0.01 0.05],'rightUpShifts',[0.16 0.11],'widthHeightAdjustment',...
        [-130 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 1.7 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_light_%s.pdf',selected_property),600);

    % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 2],...
        'spaceRowsCols',[0.01 0.05],'rightUpShifts',[0.16 0.25],'widthHeightAdjustment',...
        [-130 -320]);    set(gcf,'color','w');    set(gcf,'Position',[10 8 1.7 0.8]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_light_%s.pdf',selected_property),600);
    %%
    break;
end


%% percentage of cells sensory light different types
while 1
    an = 4;
    selected_property = 'tuned';
    titles = {'Lb','Lb*'};
    si = [Lb_T Lbs_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    resp1 = props1.good_FR_and_exc;     resp2 = props1.good_FR_and_inh;    resp3 = props1.good_FR_and_untuned;
    presp1 = 100 * exec_fun_on_cell_mat(resp1,'sum')./exec_fun_on_cell_mat(resp1,'length');
    presp2 = 100 * exec_fun_on_cell_mat(resp2,'sum')./exec_fun_on_cell_mat(resp2,'length');
    presp3 = 100 * exec_fun_on_cell_mat(resp3,'sum')./exec_fun_on_cell_mat(resp3,'length');
    
    [within,dvn,xlabels] = make_within_table({'Type','Cond'},[3,2]);
    dataT = make_between_table({presp1,presp2,presp3},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Type','hsd'},[1 1 1]);
    xdata = make_xdata([3],[1]);
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
    xticks = xdata; xticklabels = {'Exc','Sup','Com'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[0.065 0.0 -0.4 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    ht = title('Lb and Lb* (pooled)'); changePosition(ht,[-0.3 0 0])
    save_pdf(hf,mData.pdf_folder,sprintf('light_cell_types_percent.pdf'),600);
    %%
    break;
end

%% population vector and correlation sensory
while 1
    selected_property = 'untuned';
    an = 4;
    titles = {'Ab','Ab*'};
    si = [Ab_T Abs_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    resp = props1.good_FR_and_untuned;
%     eval(cmdTxt);
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 2],...
        'spaceRowsCols',[0.01 0.05],'rightUpShifts',[0.16 0.11],'widthHeightAdjustment',...
        [-130 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 1.7 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_air_%s.pdf',selected_property),600);

    % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 2],...
        'spaceRowsCols',[0.01 0.05],'rightUpShifts',[0.16 0.25],'widthHeightAdjustment',...
        [-130 -320]);    set(gcf,'color','w');    set(gcf,'Position',[10 8 1.7 0.8]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_air_%s.pdf',selected_property),600);
    %%
    break;
end


%% percentage of cells sensory air different types
while 1
    selected_property = 'tuned';
    an = 4;
    titles = {'Ab','Ab*'};
    si = [Ab_T Abs_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    resp1 = props1.good_FR_and_exc;     resp2 = props1.good_FR_and_inh;    resp3 = props1.good_FR_and_untuned;
    presp1 = 100 * exec_fun_on_cell_mat(resp1,'sum')./exec_fun_on_cell_mat(resp1,'length');
    presp2 = 100 * exec_fun_on_cell_mat(resp2,'sum')./exec_fun_on_cell_mat(resp2,'length');
    presp3 = 100 * exec_fun_on_cell_mat(resp3,'sum')./exec_fun_on_cell_mat(resp3,'length');
    
    [within,dvn,xlabels] = make_within_table({'Type','Cond'},[3,2]);
    dataT = make_between_table({presp1,presp2,presp3},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Type','hsd'},[1 1 1]);
    xdata = make_xdata([3],[1]);
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
    xticks = xdata; xticklabels = {'Exc','Sup','Com'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20 30]); xtickangle(45)
    changePosition(gca,[0.065 0.0 -0.4 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
    ht = title('Ab and Ab* (pooled)'); changePosition(ht,[-0.2 0 0])
    save_pdf(hf,mData.pdf_folder,sprintf('air_cell_types_percent_1.pdf'),600);

    %%
    break;
end


%% population vector and correlation temporal 
while 1
    titles = {'Ar-i-T','ArL-i-T','Ar*-i-T'};
    an = 4;
    si = [Ar_i_T ArL_i_T Ars_i_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    eval(cmdTxt);
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.11],'widthHeightAdjustment',...
        [-55 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.5 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,good_FR,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    changePosition(ff.h_axes(2,1).YLabel,[-1 0 0]); 
    for ii = 1:length(ff.h_axes(1,:)) set_obj(ff.h_axes(2,ii),{'xtick',[1 68 136],'xticklabels',[0 7.5 15]}); set_obj(ff.h_axes(2,ii),{'ytick',[1 68 136],'yticklabels',[0 7.5 15]}); end 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_temporal_%s.pdf',selected_property),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.25],'widthHeightAdjustment',...
        [-55 -350]);    set(gcf,'color','w');    set(gcf,'Position',[8 8 2.5 0.8]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    changePosition(ff.h_axes(1,1).YLabel,[-1 0 0]); 
    for ii = 1:length(ff.h_axes(1,:)) set_obj(ff.h_axes(1,ii),{'xtick',[1 68 136],'xticklabels',[0 7.5 15]}); set_obj(ff.h_axes(1,ii),{'ytick',[1 68 136],'yticklabels',[0 7.5 15]}); end 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_temporal_%s.pdf',selected_property),600);
    %%
    break;
end

%% population vector and correlation temporal 
while 1
    selected_property = 'tuned';
    titles = {'Ar-i-T','ArL-i-T','Ar*-i-T'};
    an = 4;
    si = [Ar_i_T ArL_i_T Ars_i_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    resp = props1.good_FR_and_Gauss_loose;
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.11],'widthHeightAdjustment',...
        [-55 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.5 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    changePosition(ff.h_axes(2,1).YLabel,[-3 0 0]); 
    for ii = 1:length(ff.h_axes(1,:)) set_obj(ff.h_axes(2,ii),{'xtick',[1 68 136],'xticklabels',[0 7.5 15]}); set_obj(ff.h_axes(2,ii),{'ytick',[1 68 136],'yticklabels',[0 7.5 15]}); end 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_temporal_%s.pdf',selected_property),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.25],'widthHeightAdjustment',...
        [-55 -350]);    set(gcf,'color','w');    set(gcf,'Position',[8 8 2.5 0.8]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    changePosition(ff.h_axes(1,1).YLabel,[-3 0 0]); 
    for ii = 1:length(ff.h_axes(1,:)) set_obj(ff.h_axes(1,ii),{'xtick',[1 68 136],'xticklabels',[0 7.5 15]}); set_obj(ff.h_axes(1,ii),{'ytick',[1 68 136],'yticklabels',[0 7.5 15]}); end 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_temporal_%s.pdf',selected_property),600);
    %%
    break;
end


%% population vector and correlation spatial 
while 1
    selected_property = 'tuned';
    titles = {'Ar-t-D','ArL-t-D','Ar*-t-D'};
    an = 4;
    si = [Ar_t_D ArL_t_D Ars_t_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    resp = props1.good_FR_and_Gauss_loose;
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.11],'widthHeightAdjustment',...
        [-55 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.5 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    changePosition(ff.h_axes(2,1).YLabel,[-3 0 0]); 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_spatial_%s.pdf',selected_property),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.25],'widthHeightAdjustment',...
        [-55 -350]);    set(gcf,'color','w');    set(gcf,'Position',[8 8 2.5 0.8]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    changePosition(ff.h_axes(1,1).YLabel,[-3 0 0]); 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_spatial_%s.pdf',selected_property),600);
    %%
    break;
end

%% percentage of gauss versus cells not gauss spatial
while 1
    an = 4;
    selected_property = 'tuned';
    si = [Ar_t_D ArL_t_D Ars_t_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    resp1 = props1.good_FR_and_Gauss_loose; resp2 = props1.good_FR_and_notGauss_loose;
    presp1 = 100 * exec_fun_on_cell_mat(resp1,'sum')./exec_fun_on_cell_mat(resp1,'length');
    presp2 = 100 * exec_fun_on_cell_mat(resp2,'sum')./exec_fun_on_cell_mat(resp2,'length');
   
    [within,dvn,xlabels] = make_within_table({'Type','Cond'},[2,3]);
    dataT = make_between_table({presp1,presp2},dvn);
    ra = RMA(dataT,within);
    ra.ranova
    %%
    [xdata,mVar,semVar,combs,p,h,colors,xlabels,extras] = get_vals_for_bar_graph_RMA(mData,ra,{'Type','hsd'},[1 1 1]);
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
    xticks = xdata; xticklabels = {'gT','gU'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0 10 20]); xtickangle(45)
    changePosition(gca,[0.065 0.0 -0.4 -0.05]); put_axes_labels(gca,{[],[0 0 0]},{'Cells (%)',[0 0 0]});
%     ht = title({'Pooled across','Conditions and Cell Types'}); changePosition(ht,[0.3 0 0])
    save_pdf(hf,mData.pdf_folder,sprintf('light_cell_types_percent.pdf'),600);
    %%
    break;
end


%% population vector and correlation spatial temporal
while 1
    titles = {'Ar-t-D','ArL-t-D','Ar*-t-D','C3it-Dist','C4it-Dist','C3''it-Dist'};
    an = 2;
    si1 = si_seq(setdiff(1:11,[1 11 9 2 10]));
    Rs1 = o.Rs(:,si1); mR1 = o.mR(:,si1);     ntrials = 50;     props1 = get_props_Rs(Rs1,ntrials);

    si = [si_air_dist_trials si_air_dist_itrials];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    eval(cmdTxt);
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 6],...
        'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.05 0.1],'widthHeightAdjustment',...
        [25 -60]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 6.99 2]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,good_FR,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
%     save_pdf(ff.hf,mData.pdf_folder,sprintf('population_vector_corr_time.pdf'),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 6],...
        'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.05 0.2],'widthHeightAdjustment',...
        [20 -240]);    set(gcf,'color','w');    set(gcf,'Position',[5 8 6.99 1]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
%     save_pdf(ff.hf,mData.pdf_folder,sprintf('average_population_vector_corr_time.pdf'),600);
    %%
    break;
end

%% population vector and correlation  temporal
while 1
    titles = {'Ar-i-T','ArL-i-T','Ar*-i-T'};
    an = 4;
    si = si_seq(setdiff(1:11,[1 11 9 2 10]));
    si = si([2 4 6]);
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    eval(cmdTxt);
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.1],'widthHeightAdjustment',...
        [-35 -60]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 3.45 2]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,props1.good_Gauss,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_temporal_%s.pdf',selected_property),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.2],'widthHeightAdjustment',...
        [-35 -240]);    set(gcf,'color','w');    set(gcf,'Position',[5 8 3.45 1]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_temporal_%s.pdf',selected_property),600);
    %%
    break;
end


%% population vector and correlation spatial with zMID larger
while 1
    titles = {'Ar-t-D','ArL-t-D','Ar*-t-D'};
    an = 4;
    si = si_seq(setdiff(1:11,[1 11 9 2 10]));
    si = si([1 3 5]);
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    eval(cmdTxt);
    good_FR = dzMI_T.resp_D_g_T;
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.1],'widthHeightAdjustment',...
        [-35 -60]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 3.45 2]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,good_FR,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_spatial_%s.pdf',selected_property),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.2],'widthHeightAdjustment',...
        [-35 -240]);    set(gcf,'color','w');    set(gcf,'Position',[5 8 3.45 1]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_spatial_%s.pdf',selected_property),600);
    %%
    break;
end

%% population vector and correlation spatial first three trials
while 1
    titles = {'Ar-t-D','ArL-t-D','Ar*-t-D'};
    an = 4;
    sel_ot = out3;
    Rs = sel_ot.Rs; mR = sel_ot.mR;
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.1],'widthHeightAdjustment',...
        [-35 -60]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 3.45 2]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,gFR_OR,1);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_spatial_trials.pdf'),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.2],'widthHeightAdjustment',...
        [-35 -240]);    set(gcf,'color','w');    set(gcf,'Position',[5 8 3.45 1]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_spatial_trials.pdf'),600);
    %%
    break;
end
