function air_figure
n = 0;
while 1
    Rs = o.Rs(:,si_air_rest);
    mR = o.mR(:,si_air_rest);
    ntrials = 50; props1 = get_props_Rs(Rs,ntrials);
    untuned = cell_list_op(props1.vals,[],'not');
    tuned = props1.vals;    good_FR = props1.good_FR;
    good_FR_tuned = cell_list_op(good_FR,tuned,'and');
    good_FR_untuned = cell_list_op(good_FR,untuned,'and');
    resp_fractionC = exec_fun_on_cell_mat(good_FR,'sum')./exec_fun_on_cell_mat(good_FR,'length')
    good_FR_OR = cell_list_op(good_FR,[],'or');
    resp_OR_fractionC = exec_fun_on_cell_mat(good_FR_OR(:,1),'sum')./exec_fun_on_cell_mat(good_FR_OR(:,1),'length');
    exc = o.resp.exc(:,si_air_rest); inh = o.resp.inh(:,si_air_rest);
%     polarity = find_resp_polarity(Rs,good_FR);
    break;
end
%%
respE_OR = cell_list_op(exc,[],'or');
respI_OR = cell_list_op(inh,[],'or');

outE = find_population_vector_corr_remap(Rs,mR,respE_OR);
outI = find_population_vector_corr_remap(Rs,mR,respI_OR);
n = 0;
%%
trials = mat2cell([1:10]',ones(size([1:10]')));
resp = get_cell_list(resp_valsC,[]);
outE1 = find_population_vector_corr_remap_trials(Rs(:,1),respE_OR,trials);
outE2 = find_population_vector_corr_remap_trials(Rs(:,2),respE_OR,trials);
outI1 = find_population_vector_corr_remap_trials(Rs(:,1),respI_OR,trials);
outI2 = find_population_vector_corr_remap_trials(Rs(:,2),respI_OR,trials);

n = 0;
%% Show stimulus train
if 1
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 0.5],'color','w');
    hold on;
    b = ei{1}.b;
    ss = b.air_puff_r(1) - round(1e6 * 3/b.si);
    ee = b.air_puff_f(10) + round(1e6 * 3/b.si);
    air_rest_sig = b.air_puff_raw(ss:ee);
    ts = b.ts(ss:ee)-b.ts(ss);
    plot(ts,air_rest_sig);
    axis off;
    text(1,1.5,{'5sec air pulse', '10sec inter-trial interval'},'FontSize',5);
    changePosition(gca,[-0.07 -0.13 -0.15 -0.1]);
    save_pdf(hf,mData.pdf_folder,sprintf('air_signal_train'),600);
end

%% Show sample rasters for presentation
% an = 5; cn = 2;
% plotRasters_simplest(Rs{an,cn})
if 1
    titles = {'Ar','Ar*'};
    si = [Ab_T Abs_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    propsA = get_props_Rs(Rs,50);
    an = 1; cn = 1;
    % plotRasters_simplest(Rs{an,cn})
    % find(resp_valsC{an}(:,cn));
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.15 0.07],'rightUpShifts',[0.07 0.25],'widthHeightAdjustment',...
        [-100 -475]);
    set(gcf,'color','w'); set(gcf,'Position',[10 4 4 1]);
%     ff = sample_rasters(Rs{an,cn},[140 17 319 168],ff);
    ff = sample_rasters(Rs{an,cn},[45 8 25 181],ff);
    for ii = 1:length(ff.h_axes)
        set_obj(ff.h_axes(ii),{'xticklabels',{'-5','0','5','10'}})
    end
    colormap parula
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rasters'),600);
end

%% Show sample rasters
% an = 5; cn = 2;
% plotRasters_simplest(Rs{an,cn})
if 1
    titles = {'Ar','Ar*'};
    si = [Ab_T Abs_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    an = 1; cn = 1;
    % plotRasters_simplest(Rs{an,cn})
    % find(resp_valsC{an}(:,cn));
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 2],...
        'spaceRowsCols',[0.15 0.1],'rightUpShifts',[0.15 0.25],'widthHeightAdjustment',...
        [-175 -475]);
    set(gcf,'color','w'); set(gcf,'Position',[10 4 2 1]);
    ff = sample_rasters(Rs{an,cn},[140 181],ff);
%     ff = sample_rasters(Rs{an,cn},[140 66 162 181],ff);
    for ii = 1:length(ff.h_axes)
        set_obj(ff.h_axes(ii),{'xticklabels',{'-5','0','5','10'}})
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rasters'),600);
end

%% population vector and correlation single animal
while 1
    an = 1;
    ff = makeFigureRowsCols(106,[1 0.5 4 1],'RowsCols',[2 2],...
        'spaceRowsCols',[0 -0.02],'rightUpShifts',[0.15 0.1],'widthHeightAdjustment',...
        [-80 -70]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 5 2.2 2]);
%     [resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC,resp_exc_inh] = get_responsive_fraction(Rs);
%     resp = get_cell_list(resp_valsC,[1;2]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,good_FR_untuned,0);
    % ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[-0.1 1],[]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[-0.1 1],[]);
    set_obj(ff,{'FontWeight','Normal','FontSize',6,'LineWidth',0.25});
    ht = get_obj(ff,'title'); hyl = get_obj(ff,'ylabel'); changePosition(hyl(1,1),[7 0 0]);
    set_obj(ht,'String',{'Pop. Activity','Pop. Activity';'Pop. Correlation','Pop.Correlation'});
    set_obj(ht,{'FontSize',5,'FontWeight','Normal'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_population_vector_corr.pdf'),600);
    break;
end
%% average population correlation (from all animals)
while 1
    ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 2],...
        'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.15 0.2],'widthHeightAdjustment',...
        [-70 -240]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 3 2.2 1]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[-0.1 1],[]);
    set_obj(ff,{'FontWeight','Normal','FontSize',6,'LineWidth',0.25});
    ht = get_obj(ff,'title'); hyl = get_obj(ff,'ylabel'); changePosition(hyl(1,1),[7 0 0]);
    set_obj(ht,'String',{'Avg. Pop. Correlation','Avg. Pop.Correlation'});
    set_obj(ht,{'FontSize',5,'FontWeight','Normal'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_average_population_vector_corr.pdf'),600);
    break;
end
%% population correlation clustering and showing that mean of two clusters are significantly different across animals
if 1
    [resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC,resp_exc_inh] = get_responsive_fraction(Rs);
    resp = get_cell_list(resp_valsC,[1;2]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,resp,0);
    mask = logical(triu(ones(size(CRc{1,1})),1));
    for rr = 1:size(CRc,1)
        for cc = 1:size(CRc,2)
            tempC = CRc{rr,cc};
            vals = tempC;%(mask);
            [clusi,clusc] = kmeans(vals,2);
%             rng('default');  % For reproducibility
%             eva = evalclusters(vals,'kmeans','gap','KList',[1:4])
%             optimalK(rr,cc) = eva.OptimalK;
            if clusc(1) < clusc(2)
                clusS(rr,cc) = clusc(1);
                clusB(rr,cc) = clusc(2);
                clusi1p(rr,cc) = 100*(sum(clusi == 1)/length(clusi));
                clusi2p(rr,cc) = 100*(sum(clusi == 2)/length(clusi));
            else
                clusS(rr,cc) = clusc(2);
                clusB(rr,cc) = clusc(1);
                clusi1p(rr,cc) = 100*(sum(clusi == 2)/length(clusi));
                clusi2p(rr,cc) = 100*(sum(clusi == 1)/length(clusi));
            end
        end
    end
    [within,dvn,xlabels] = make_within_table({'Cond','clusInd'},2,2);
    dataT = make_between_table({clusi1p(:,1),clusi2p(:,1),clusi1p(:,2),clusi2p(:,2)},dvn);
    ra = repeatedMeasuresAnova(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',colors,'sigColor','k',...
        'ySpacing',0.02,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.05);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = [xdata(1:end)]; xticklabels = xlabels;
    set(gca,'xtick',xticks,'xticklabels',xticklabels);  xtickangle(30)
    changePosition(gca,[0.2 0.03 -0.2 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Air Responsive','Cells (%)'},[0 0 0]});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('air_population_vector_corr_clus_ind_bar_garph.pdf'),600);
end
%% Percentage of Responsive Cells
if 1
    within = make_within_table({'Cond'},2);
    dataT = make_between_table({resp_fractionC*100},{'C21','C22'});
    ra = RMA(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph_RMA(mData,ra,{'Cond','hsd'},[1 1 1]);
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',colors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.05);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = [xdata(1:end)]; xticklabels = {'C2','C2'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30)
    any_mean = mean(100*resp_OR_fractionC);    any_sem = std(100*resp_OR_fractionC)/sqrt(5);
    pmchar=char(177); any_text = sprintf('%.0f%c%.0f%%',any_mean,pmchar,any_sem); text(0.75,45,any_text,'FontSize',6);
    changePosition(gca,[0.2 0.03 -0.5 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Air Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('percentage_air_responsive'),600);
end
%% Percentage of excitatory inhibitory responsive cells
if 1
    [respE1,respE_OR1,respE_AND1,respE_fraction1] = get_cell_list_exc_inh(resp_exc_inh,1,1);
    [respE2,respE_OR2,respE_AND2,respE_fraction2] = get_cell_list_exc_inh(resp_exc_inh,2,1);
    
    [respI1,respI_OR1,respI_AND1,respI_fraction1] = get_cell_list_exc_inh(resp_exc_inh,1,0);
    [respI2,respI_OR2,respI_AND2,respI_fraction2] = get_cell_list_exc_inh(resp_exc_inh,2,0);
    [within,dvn,xlabels] = make_within_table({'Cond','EI'},[2,2]);
    between = make_between_table({100*respE_fraction1',100*respI_fraction1',100*respE_fraction2',100*respI_fraction2'},dvn);
    ra = repeatedMeasuresAnova(between,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0 1]);
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');hold on;
%     colors = mData.colors;
%     tcolors = {colors{3};colors{3};colors{4};colors{4};};
    tcolors = {colors{1},colors{1}/3,colors{2},colors{2}/3};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     for ii = [2 4]
%         set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
%     end
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 30],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C2-Ex','C2-Su','C2''-Ex','C2''-Su'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.2 0.03 -0.25 -0.05]);
%     changePosition(gca,[0.11 0.03 -0.2 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{[],[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'air_responsive_exc_inh',600);
end

%% Spatial correlation between adjacent trails and considering exc and inh as two separate groups
if 1
    [within,dvn,xlabels] = make_within_table({'Condition','TrialPairs'},[2,9]);
    varE1 = arrayfun(@(x) mean(x{1}),outE1.adj_SP_corr_diag); varE2 = arrayfun(@(x) mean(x{1}),outE2.adj_SP_corr_diag);
    varI1 = arrayfun(@(x) mean(x{1}),outI1.adj_SP_corr_diag); varI2 = arrayfun(@(x) mean(x{1}),outI2.adj_SP_corr_diag);

%     varE1 = arrayfun(@(x) mean(x{1}),outE1.adj_PV_corr_diag); varE2 = arrayfun(@(x) mean(x{1}),outE2.adj_PV_corr_diag);
%     varI1 = arrayfun(@(x) mean(x{1}),outI1.adj_PV_corr_diag); varI2 = arrayfun(@(x) mean(x{1}),outI2.adj_PV_corr_diag);between = make_between_table({varE1,varE2;varI1,varI2},dvn);

%     varE1 = arrayfun(@(x) nanmean(x{1}),outE1.adj_RR_SP); varE2 = arrayfun(@(x) nanmean(x{1}),outE2.adj_RR_SP);
%     varI1 = arrayfun(@(x) nanmean(x{1}),outI1.adj_RR_SP); varI2 = arrayfun(@(x) nanmean(x{1}),outI2.adj_RR_SP);
    
    between = make_between_table({varE1,varI1;varE2,varI2},dvn);
    ra = repeatedMeasuresAnova(between,within);
    %%
    hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 2 1],'color','w');
    hold on;
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_line_graph(mData,ra,0);
    tcolors = {colors{1},colors{1}/3,colors{2},colors{2}/3};
    ii = 1; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 2; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linestyle','-.','linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 3; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.25,'linestyle','none','capsize',1);
    ii = 4; plot(xdata,mVar(ii,:),'color',tcolors{ii},'linestyle','-.','linewidth',0.5); errorbar(xdata,mVar(ii,:),semVar(ii,:),'color',tcolors{ii},'linewidth',0.5,'linestyle','none','capsize',1);
    legs = {'C2-Ex','C2-Su','C2''-Ex','C2''-Su'};
    legs{end+1} = [0.5 0.2 0.07 0.25];
    putLegendH(gca,legs,tcolors)
    ylims = ylim;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[ylims(1) 0.08],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'T1-T2','T2-T3','T3-T4','T4-T5','T5-T6','T6-T7','T7-T8','T8-T9','T9-T10'};
    xticklabels = repmat(xticklabels,1,2);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.03 0.11 0.05 -0.07]);
    put_axes_labels(gca,{'Trial-Pairs',[0 0 0]},{{'Correlation'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'spatial_correlation_trials',600);
end
%% Spatial correlation between adjacent trails (taking mean) and considering exc and inh as within factors (not groups)
if 1
    [within,dvn,xlabels] = make_within_table({'Condition','EI'},[2 2]);
    varE1 = mean(arrayfun(@(x) mean(x{1}),outE1.adj_SP_corr_diag),2); varE2 = mean(arrayfun(@(x) mean(x{1}),outE2.adj_SP_corr_diag),2);
    varI1 = mean(arrayfun(@(x) mean(x{1}),outI1.adj_SP_corr_diag),2); varI2 = mean(arrayfun(@(x) mean(x{1}),outI2.adj_SP_corr_diag),2);
    between = make_between_table({varE1,varI1,varE2,varI2},dvn);
    ra = repeatedMeasuresAnova(between,within);

    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0.5 1]);
    hollowsep = 19;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors = {colors{1},colors{1}/3,colors{2},colors{2}/3};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
            'ySpacing',0.01,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
            'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = hollowsep:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',colors{ii});
    end
    ylims = ylim;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[ylims(1) maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C2-Ex','C2-Su','C2''-Ex','C2''-Su'};
    xticklabels = repmat(xticklabels,1,2);
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
%     changePosition(gca,[0.1 0.11 -0.2 -0.05]);
    changePosition(gca,[0.2 0.03 -0.2 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Correlation'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'spatial_correlation_trials_bar_graph',600);
end

%% Rate change score between adjacent trails (taking mean) and considering exc and inh as within factors (not groups)
if 1
    [within,dvn,xlabels] = make_within_table({'Condition','EI'},[2 2]);
    varE1 = mean(arrayfun(@(x) nanmean(x{1}),outE1.adj_RR_SP),2); varE2 = mean(arrayfun(@(x) nanmean(x{1}),outE2.adj_RR_SP),2);
    varI1 = mean(arrayfun(@(x) nanmean(x{1}),outI1.adj_RR_SP),2); varI2 = mean(arrayfun(@(x) nanmean(x{1}),outI2.adj_RR_SP),2);
    between = make_between_table({varE1,varI1,varE2,varI2},dvn);
    ra = repeatedMeasuresAnova(between,within);

    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0.5 1]);
    hollowsep = 19;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors = {colors{1},colors{1}/3,colors{2},colors{2}/3};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
            'ySpacing',0.03,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
            'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.01);
    for ii = hollowsep:length(hbs)
        set(hbs(ii),'facecolor','none','edgecolor',colors{ii});
    end
    ylims = ylim;
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0.92 1.1],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C2-Ex','C2-Su','C2''-Ex','C2''-Su'};
    xticklabels = repmat(xticklabels,1,2);
    set(gca,'xtick',xticks,'xticklabels',xticklabels,'ytick',[0.92 1 1.08]);
    xtickangle(30);
%     changePosition(gca,[0.1 0.11 -0.2 -0.05]);
    changePosition(gca,[0.2 0.03 -0.2 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'\Delta FR Score'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'Rate_change_score_trials_bar_graph',600);
end


%% Spatial Correlation Between Conditions of both excitatory and inhibitory responses
if 1
    within = make_within_table({'EI'},2);
    var_CE = arrayfun(@(x) mean(x{1}),outE.adj_SP_corr_diag);var_CI = arrayfun(@(x) mean(x{1}),outI.adj_SP_corr_diag);
%     var_CE = arrayfun(@(x) mean(x{1}),outE.adj_PV_corr_diag);var_CI = arrayfun(@(x) mean(x{1}),outI.adj_PV_corr_diag);
%     var_CE = arrayfun(@(x) mean(x{1}),outE.adj_RR_SP);var_CI = arrayfun(@(x) mean(x{1}),outI.adj_RR_SP);
    dataT = make_between_table({var_CE,var_CI},{'E','I'});
    ra = repeatedMeasuresAnova(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0 1]);
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors ={colors{7};colors{8};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.02,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C2-C2''-Ex','C2-C2''-Su'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     set(hbs(2),'facecolor','none','edgecolor',tcolors{2});
    xtickangle(30);
%     changePosition(gca,[0.2 0.03 -0.55 -0.05]);
    changePosition(gca,[0.1 0.03 -0.4 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Correlation'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'air_rest_cell corr remap',600);
return;
end

%% Delta FR Score Between Conditions of both excitatory and inhibitory responses
if 1
    within = make_within_table({'EI'},2);
%     var_CE = arrayfun(@(x) mean(x{1}),outE.adj_PV_corr_diag);var_CI = arrayfun(@(x) mean(x{1}),outI.adj_PV_corr_diag);
    var_CE = arrayfun(@(x) mean(x{1}),outE.adj_RR_SP);var_CI = arrayfun(@(x) mean(x{1}),outI.adj_RR_SP);
    dataT = make_between_table({var_CE,var_CI},{'E','I'});
    ra = repeatedMeasuresAnova(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 0 1]);
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors ={colors{7};colors{8};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.04,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.01);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0.9 1.03],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C2-C2''-Ex','C2-C2''-Su'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     set(hbs(2),'facecolor','none','edgecolor',tcolors{2});
    xtickangle(30);
%     changePosition(gca,[0.2 0.03 -0.55 -0.05]);
    changePosition(gca,[0.1 0.03 -0.4 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'\Delta FR Score'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'air_rest_DeltaFR_score remap',600);
return;
end

%% average correlation of all animals
ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[2 2],...
    'spaceRowsCols',[0.13 0.13],'rightUpShifts',[0.13 0.2],'widthHeightAdjustment',...
    [-220 -220]);
set(gcf,'color','w');
set(gcf,'Position',[5 3 2 2]);
ff = show_remapping_corr_plots(mData,ff,mean_corr_CA,xs_CA,[-0.1 1]);
save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rest_remap_corr.pdf'),600);

%% overlap RM bar graph
if 1
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    
    OI_E = get_overlap_index(resp_exc_inh(:,1));
    OI_I = get_overlap_index(resp_exc_inh(:,2));
    for ii = 1:5
        C12(ii,1) = OI_E{ii}(1,2);
        C12(ii,2) = OI_I{ii}(1,2);
    end

    dataT = array2table([C12]);
    dataT.Properties.VariableNames = {'E','I'};
    colVar1 = [1 2];    
    within = table(colVar1');
    within.Properties.VariableNames = {'EI'};
    within.EI = categorical(within.EI);
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    
    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
    xdata = [1 2];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    tcolors ={colors{7};colors{8};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.02,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'Ex','Su'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     set(hbs(2),'facecolor','none','edgecolor',tcolors{2});
    xtickangle(30)
    xticks = xdata(1:end)+0; xticklabels = {'C2-C2''-Ex','C2-C2''-Su'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     set(hbs(2),'facecolor','none','edgecolor',tcolors{2});
    xtickangle(30);
%     changePosition(gca,[0.2 0.03 -0.55 -0.05]);
    changePosition(gca,[0.1 0.03 -0.4 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Overlap Index'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('air_overlap_stats'),600);
end
%% pca

tmR = mR{1,1};
[coeff,score,latent,tsquared,explained,mu] = pca(tmR');
cexplained = cumsum(explained);
ind = find(cexplained > 95,1,'first')
