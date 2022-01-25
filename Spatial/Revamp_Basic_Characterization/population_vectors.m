%% get cell populations
ntrials = 50;
si = [Lb_T ArL_L_T Lbs_T Ab_T Abs_T Ar_t_D ArL_t_D Ars_t_D Ar_i_T ArL_i_T Ars_i_T];
props1 = get_props_Rs(o.Rs(:,si),ntrials);
resp = props1.good_FR;
only_cells = [];
for ii = 1:length(si)
    cond_mat = -(1:length(si));
    cond_mat(ii) = -cond_mat(ii);
    only_cells = [only_cells get_cell_list(resp,cond_mat,0)];
    
end
perc_only_cells = 100*exec_fun_on_cell_mat(only_cells,'sum')./exec_fun_on_cell_mat(only_cells,'length');

% only light
cond_mat = [1;3];
temp_cell_list = get_cell_list(resp,cond_mat,0);
int_cell_list = [temp_cell_list resp(:,setdiff(1:11,cond_mat))];
only_ones = get_cell_list(int_cell_list,[1 -(2:size(int_cell_list,2))]);
perc_only_ones = 100*exec_fun_on_cell_mat(only_ones,'sum')./exec_fun_on_cell_mat(only_ones,'length');
only_light = only_ones;
perc_only_light = perc_only_ones;


% only air
cond_mat = [4;5];
temp_cell_list = get_cell_list(resp,cond_mat,0);
int_cell_list = [temp_cell_list resp(:,setdiff(1:11,cond_mat))];
only_ones = get_cell_list(int_cell_list,[1 -(2:size(int_cell_list,2))]);
perc_only_ones = 100*exec_fun_on_cell_mat(only_ones,'sum')./exec_fun_on_cell_mat(only_ones,'length');
only_air = only_ones;
perc_only_air = perc_only_ones;

% only spatial
cond_mat = [6;7;8];
temp_cell_list = get_cell_list(resp,cond_mat,0);
int_cell_list = [temp_cell_list resp(:,setdiff(1:11,cond_mat))];
only_ones = get_cell_list(int_cell_list,[1 -(2:size(int_cell_list,2))]);
perc_only_ones = 100*exec_fun_on_cell_mat(only_ones,'sum')./exec_fun_on_cell_mat(only_ones,'length');
only_spatial = only_ones;
perc_only_spatial = perc_only_ones;

% only temporal
cond_mat = [9;10;11];
temp_cell_list = get_cell_list(resp,cond_mat,0);
int_cell_list = [temp_cell_list resp(:,setdiff(1:11,cond_mat))];
only_ones = get_cell_list(int_cell_list,[1 -(2:size(int_cell_list,2))]);
perc_only_ones = 100*exec_fun_on_cell_mat(only_ones,'sum')./exec_fun_on_cell_mat(only_ones,'length');
only_temporal = only_ones;
perc_only_temporal = perc_only_ones;

%% population vector and correlation sensory
while 1
    an = 4;
    titles = {'Lb','ArL-L','Lb*'};
    si = [Lb_T ArL_L_T Lbs_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    eval(cmdTxt);
    ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.13],'widthHeightAdjustment',...
        [-50 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.5 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,good_FR,0);
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

%% population vector and correlation sensory tuned
while 1
    an = 4;
    selected_property = 'tuned';
    titles = {'Lb','ArL-L','Lb*'};
    si = [Lb_T ArL_L_T Lbs_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    resp = props1.good_FR_and_tuned;
%     resp = props1.good_FR_and_inh;
    ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.13],'widthHeightAdjustment',...
        [-50 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.5 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) set_obj(ff.h_axes(2,ii),{'xtick',[1 18 38],'xticklabels',[-2 0 2.2]}); set_obj(ff.h_axes(2,ii),{'ytick',[1 18 38],'yticklabels',[-2 0 2.2]}); end 
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    changePosition(ff.h_axes(2,1).YLabel,[-2.5 0 0]); 
    colormap parula
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_light_%s.pdf',selected_property),600);

    % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.27],'widthHeightAdjustment',...
        [-50 -370]);    set(gcf,'color','w');    set(gcf,'Position',[10 8 2.5 0.85]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    for ii = 1:length(ff.h_axes(1,:)) set_obj(ff.h_axes(1,ii),{'xtick',[1 18 38],'xticklabels',[-2 0 2.2]}); set_obj(ff.h_axes(1,ii),{'ytick',[1 18 38],'yticklabels',[-2 0 2.2]}); end 
    changePosition(ff.h_axes(1,1).YLabel,[-2.5 0 0]); 
    colormap parula
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

%% population vector and correlation sensory air trials and intertrials separate
while 1
    selected_property = 'untuned';
    an = 4;
    titles = {'Ab-t-T','Ab*-t-T','Ab-i-T','Ab*-i-T'};
    si = [Ab_t_T Abs_t_T Ab_i_T Abs_i_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    respA = props1.good_FR_and_tuned;
    resp1 = cell_list_op(respA(:,1:2),[],'or');
    resp2 = cell_list_op(respA(:,3:4),[],'or');
    resp = respA;%[resp1(:,1:2),resp2(:,1:2)];
%     eval(cmdTxt);
    ffM = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 4],...
        'spaceRowsCols',[0.01 0.03],'rightUpShifts',[0.08 0.11],'widthHeightAdjustment',...
        [-40 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 3 1.5]);
    sInds = 1:2;
    [CRc,aCRc1,mRR] = find_population_vector_corr(Rs(:,sInds),mR(:,sInds),resp(:,sInds),0);
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),mRR(an,:),CRc(an,:),[],[],0);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    
    sInds = 3:4;
    [CRc,aCRc2,mRR] = find_population_vector_corr(Rs(:,sInds),mR(:,sInds),resp(:,sInds),0);
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    
    
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_air_%s.pdf',selected_property),600);

    % average correlation of all animals
    ffM = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 4],...
        'spaceRowsCols',[0.01 0.03],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
        [-60 -320]);    set(gcf,'color','w');    set(gcf,'Position',[10 8 3 0.8]);
    sInds = 1:2;
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),[],aCRc1,[],[],0);
    sInds = 3:4;
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),[],aCRc1,[],[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_air_%s.pdf',selected_property),600);
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
    ntrials = [30,100];
    props1 = get_props_Rs(Rs,ntrials);
    resp = props1.good_FR_and_tuned;
%     resp = only_air(:,[1 2]);
%     eval(cmdTxt);
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 2],...
        'spaceRowsCols',[0.01 0.05],'rightUpShifts',[0.16 0.11],'widthHeightAdjustment',...
        [-130 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 1.7 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    colormap parula
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_air_%s.pdf',selected_property),600);

    % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 2],...
        'spaceRowsCols',[0.01 0.05],'rightUpShifts',[0.16 0.25],'widthHeightAdjustment',...
        [-130 -320]);    set(gcf,'color','w');    set(gcf,'Position',[10 8 1.7 0.8]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    colormap parula
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
    Rs = o.Rs(:,si); mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
%     resp = cell_list_op(props1.good_FR,[],'and');
    resp = props1.good_FR;
%     eval(cmdTxt);
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.11],'widthHeightAdjustment',...
        [-55 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.5 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,1);
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
    resp = props1.good_FR_and_tuned;
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


%% population vector and correlation spatial for presentation
while 1
    selected_property = 'tuned';
    titles = {'Ar-t-D','Ar-t-D','Ar-t-D'};
    an = 4;
    si = [Ar_t_D Ar_t_D Ar_t_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    q_type = 'All';
    resp = [props1.good_FR(:,1) props1.good_zMI(:,1) props1.good_Gauss(:,1)]; 
%     resp = cell_list_op(props1.good_FR_and_Gauss_loose,[],'or');
%     resp = cell_list_op(props1.good_FR,[],'or');
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.11],'widthHeightAdjustment',...
        [-55 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.5 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    changePosition(ff.h_axes(2,1).YLabel,[-3 0 0]); 
    colormap parula
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_spatial_%s_%s.pdf',selected_property,q_type),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.25],'widthHeightAdjustment',...
        [-55 -350]);    set(gcf,'color','w');    set(gcf,'Position',[8 8 2.5 0.8]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    changePosition(ff.h_axes(1,1).YLabel,[-3 0 0]); 
    colormap parula
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_spatial_%s_%s.pdf',selected_property,q_type),600);
    %%
    break;
end


%% population vector and correlation spatial 
while 1
    selected_property = 'tuned';
    titles = {'Ar-t-D','ArL-t-D','Ar*-t-D'};
    an = 3;
    si = [Ar_t_D ArL_t_D Ars_t_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    props1 = get_props_Rs(Rs,[50,100]);
    q_type = '1040';
%     resp = cell_list_op(props1.good_FR,[],'and');
    resp = props1.good_FR;
%     resp = repmat(respSen(:,2),1,3);
%     resp = cell_list_op(props1.good_FR_and_Gauss_loose,[],'or');
%     resp = cell_list_op(props1.good_FR,[],'or');
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.11],'widthHeightAdjustment',...
        [-55 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.5 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    changePosition(ff.h_axes(2,1).YLabel,[-3 0 0]); 
    colormap parula
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_spatial_%s_%s.pdf',selected_property,q_type),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.25],'widthHeightAdjustment',...
        [-55 -350]);    set(gcf,'color','w');    set(gcf,'Position',[8 8 2.5 0.8]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    changePosition(ff.h_axes(1,1).YLabel,[-3 0 0]); 
    colormap parula
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_spatial_%s_%s.pdf',selected_property,q_type),600);
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
    titles = {'Ar-t-D','ArL-t-D','Ar*-t-D','Ar-i-T','ArL-i-T','Ar*-i-T'};
    an = 4;
%     resp = [propsD.good_FR(:,1:3) propsT.good_FR(:,4:6)];
    respDo = propsD.good_Gauss_loose(:,1:3); respTo = propsT.good_Gauss_loose(:,4:6);
    respDi = cell_list_op(respDo,[],'or'); respTi = cell_list_op(respTo,[],'or');
    respD = cell_list_op(respDi,respTi,'sep');
    respT = cell_list_op(respTi,respDi,'sep');
    orderD = 1; orderT = 1;
%     resp = [respD(:,1:3),respT(:,1:3)];
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 6],...
        'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.055 0.1],'widthHeightAdjustment',...
        [25 -60]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 5 1.5]);
    ffL = ff; ffL.axesPos = ffL.axesPos(:,1:3); ffL.h_axes = ffL.h_axes(:,1:3);
    [CRc,aCRcD,mRR] = find_population_vector_corr(RsD(:,1:3),mRD(:,1:3),respD,orderD);
    ffL = show_population_vector_and_corr(mData,ffL,RsD(an,1:3),mRR(an,:),CRc(an,:),[],[],0);
    for ii = 1:length(ffL.h_axes(1,:)) ht = get_obj(ffL.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    ffL = ff; ffL.axesPos = ffL.axesPos(:,4:6); ffL.h_axes = ffL.h_axes(:,4:6);
    [CRc,aCRcT,mRR] = find_population_vector_corr(RsT(:,4:6),mRT(:,4:6),respT,orderT);
    ffL = show_population_vector_and_corr(mData,ffL,RsT(an,4:6),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ffL.h_axes(1,:)) ht = get_obj(ffL.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('population_vector_corr_dist_time.pdf'),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 6],...
        'spaceRowsCols',[0 -0.024],'rightUpShifts',[0.055 0.25],'widthHeightAdjustment',...
        [20 -310]);    set(gcf,'color','w');    set(gcf,'Position',[10 8 5 0.85]);
    ffL = ff; ffL.axesPos = ffL.axesPos(:,1:3); ffL.h_axes = ffL.h_axes(:,1:3);
    ffL = show_population_vector_and_corr(mData,ffL,RsD(an,1:3),[],aCRcD,[],[],0);
    ffL = ff; ffL.axesPos = ffL.axesPos(:,4:6); ffL.h_axes = ffL.h_axes(:,4:6);
    ffL = show_population_vector_and_corr(mData,ffL,RsT(an,4:6),[],aCRcT,[],[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('average_population_vector_corr_dist_time.pdf'),600);
    %%
    break;
end

%% population vector and correlation  temporal
while 1
    selected_property = 'tuned';
    titles = {'Ar-i-T','ArL-i-T','Ar*-i-T'};
    an = 4;
    si = [Ar_i_T ArL_i_T Ars_i_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    q_type = 'Resp';
    resp = props1.good_FR_and_tuned;
%     resp = cell_list_op(props1.good_FR,[],'or');
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.01],'rightUpShifts',[0.11 0.1],'widthHeightAdjustment',...
        [-35 -60]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 3.45 2]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
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


%% population vector and correlation temporal 
while 1
    selected_property = 'tuned';
    titles = {'Ar-i-T','ArL-i-T','Ar*-i-T'};
    an = 4;
    si = [Ar_i_T ArL_i_T Ars_i_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    q_type = 'Resp_Gauss_Loose';
    props1 = get_props_Rs(Rs,ntrials);
    resp = props1.good_FR_IT;
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.11],'widthHeightAdjustment',...
        [-55 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.5 1.5]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[],[]);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set_obj(ht,{'String',titles{ii}}); end
    changePosition(ff.h_axes(2,1).YLabel,[-3 0 0]); 
    for ii = 1:length(ff.h_axes(1,:)) set_obj(ff.h_axes(2,ii),{'xtick',[1 68 136],'xticklabels',[0 7.5 15]}); set_obj(ff.h_axes(2,ii),{'ytick',[1 68 136],'yticklabels',[0 7.5 15]}); end 
%     colormap parula
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_temporal_%s_%s.pdf',selected_property,q_type),600);

   % average correlation of all animals
    ff = makeFigureRowsCols(108,[1 0.5 4 0.5],'RowsCols',[1 3],...
        'spaceRowsCols',[0 0.03],'rightUpShifts',[0.11 0.25],'widthHeightAdjustment',...
        [-55 -350]);    set(gcf,'color','w');    set(gcf,'Position',[8 8 2.5 0.8]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[],[]);
    changePosition(ff.h_axes(1,1).YLabel,[-3 0 0]); 
    for ii = 1:length(ff.h_axes(1,:)) set_obj(ff.h_axes(1,ii),{'xtick',[1 68 136],'xticklabels',[0 7.5 15]}); set_obj(ff.h_axes(1,ii),{'ytick',[1 68 136],'yticklabels',[0 7.5 15]}); end 
%     colormap parula
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_temporal_%s_%s.pdf',selected_property,q_type),600);
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



%% population vector and correlation spatial for presentation
while 1
    selected_property = 'tuned';
    titles = {'Ar-t-D','ArL-t-D','Ar*-t-D'};
    an = 5;
    si = [Ar_t_D ArL_t_D Ars_t_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    resp = props1.good_FR_and_Gauss_loose;
    resp = cell_list_op(props1.good_FR_and_Gauss_loose,[],'or');
    resp = cell_list_op(props1.good_FR,[],'or');
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


%% population vector grand
while 1
    selected_property = 'tuned';
%     titles = {'Ar-t-D','ArL-t-D','Ar*-t-D'};
    titles = {'Lb','ArL-L','Lb*','Ab-t','Ab-i','Ab*-t','Ab*-i','Ar-t','ArL-t','Ar*-t','Ar-i','ArL-i','Ar*-i'};
    an = 4;
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = [50,100];
    props1 = get_props_Rs(Rs,ntrials);
    resp = props1.good_FR;
    resp = get_cell_list(resp,[1;3]);
    resp_perc = get_cell_list(resp,[1;3],1);
    resp_perc(an,:)
%     resp = cell_list_op(props1.good_FR,[],'or');
    szC = [];
    for cc = 1:size(resp,2)
        szC(1,cc) = size(mR{an,cc},2);
    end
    numCells = sum(resp{an,1},1); numConds = size(resp,2);
    giantmR = NaN(numCells,sum(szC));
    for cc = 1:size(resp,2)
        if cc == 1
            colS = 1; colE = szC(1);
        else
            colS = colE+1; colE = colE+szC(cc);
        end
        giantmR(1:numCells,colS:colE) = mR{an,cc}(resp{an,cc},:);
        stC = colS+1;
    end
    figure(10000);clf;imagesc(corr(giantmR));
    hold on;
    for cc = 1:size(resp,2)
        if cc == 1
            colS = 1; colE = szC(1);
        else
            colS = colE+1; colE = colE+szC(cc);
        end
        plot([0 sum(szC)],[colE colE],'r');
        plot([colE colE],[0 sum(szC)],'r');
    end
    
    %%
    break;
end
