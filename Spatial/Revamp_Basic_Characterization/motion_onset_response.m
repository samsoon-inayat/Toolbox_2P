function firing_rate_motion_vs_rest

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_C = evalin('base','ei'); 

selContexts = [0,0];
rasterNames = {'motionOnsets','motionOffsets'};

Rs = get_rasters_data(ei_C,selContexts,rasterNames);
Rs = find_responsive_rasters(Rs,[]);
mR = calc_mean_rasters(Rs,[]);
[resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC,resp_exc_inh] = get_responsive_fraction(Rs);

[respE,respE_OR,respE_AND,respE_fraction] = get_cell_list_exc_inh(resp_exc_inh,1,1);
[CRcE,aCRcE,mRRE] = find_population_vector_corr(Rs,mR,respE,0);

[respI,respI_OR,respI_AND,respI_fraction] = get_cell_list_exc_inh(resp_exc_inh,1,0);
[CRcI,aCRcI,mRRI] = find_population_vector_corr(Rs,mR,respI,0);

n = 0;
%% Show sample rasters
an = 3; cn = 1;
plotRasters_simplest(Rs{an,cn})
%%
if 1
    an = 3; cn = 1;
    % plotRasters_simplest(Rs{an,cn})
    % find(resp_valsC{an}(:,cn));
    ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
        'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
        [-75 -475]);
    set(gcf,'color','w'); set(gcf,'Position',[10 4 3.25 1]);
    ff = sample_rasters(Rs{an,cn},[3 36 38 226],ff);
    axes(ff.h_axes(1,1));
%     text(0,13.5,{'Representative rasters - Condition C2'},'FontSize',7,'FontWeight','Normal');
    save_pdf(ff.hf,mData.pdf_folder,sprintf('motionOnset_rasters'),600);
end
%% population vector and correlation single animal
if 1
    an = 1;
    ff = makeFigureRowsCols(106,[1 0.5 4 1],'RowsCols',[2 2],...
        'spaceRowsCols',[0 -0.02],'rightUpShifts',[0.15 0.1],'widthHeightAdjustment',...
        [-80 -70]);
    set(gcf,'color','w');
    set(gcf,'Position',[5 5 2.2 2]);
    [resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC,resp_exc_inh] = get_responsive_fraction(Rs);
%     resp = get_cell_list(resp_valsC,[1;2]);
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,1,0);
    % ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[-0.1 1],[]);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[-0.1 1],[]);
    set_obj(ff,{'FontWeight','Normal','FontSize',6,'LineWidth',0.25});
    ht = get_obj(ff,'title'); hyl = get_obj(ff,'ylabel'); changePosition(hyl(1,1),[0 0 0]);
    set_obj(ht,'String',{'Pop. Activity','Pop. Activity';'Pop. Correlation','Pop.Correlation'});
    set_obj(ht,{'FontSize',5,'FontWeight','Normal'});
    save_pdf(ff.hf,mData.pdf_folder,sprintf('motion_population_vector_corr.pdf'),600);
end


%% average population correlation (from all animals)
if 1
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
    save_pdf(ff.hf,mData.pdf_folder,sprintf('motion_average_population_vector_corr.pdf'),600);
end
%% Percentage of Responsive Cells
if 1
    within = make_within_table({'Cond'},2);
    dataT = make_between_table({resp_fractionC*100},{'M_On','M_Off'});
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w'); hold on;
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',colors,'sigColor','k',...
        'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.05);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = [xdata(1:end)]; xticklabels = {'C2','C2'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30)
    any_mean = mean(100*resp_OR_fractionC);    any_sem = std(100*resp_OR_fractionC)/sqrt(5);
    pmchar=char(177); any_text = sprintf('%.0f%c%.0f%%',any_mean,pmchar,any_sem); %text(0.75,31,any_text,'FontSize',6);
    changePosition(gca,[0.2 0.03 -0.5 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Air Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('percentage_motion_responsive'),600);
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
        'ySpacing',5,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.1);
%     for ii = [2 4]
%         set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
%     end
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 30],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'MOn-Ex','MOn-Su','MOff-Ex','MOff-Su'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.2 0.03 -0.25 -0.05]);
%     changePosition(gca,[0.11 0.03 -0.2 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{[],[0 0 0]});
    save_pdf(hf,mData.pdf_folder,'motion_responsive_exc_inh',600);
end
%%
if 1
    [within,dvn,xlabels] = make_within_table({'Conditions'},7);
    dataT = make_between_table({out_CC.mean_mcr*100},dvn);
    ra = repeatedMeasuresAnova(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
    
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 9 1.35 1],'color','w');

    hold on;
%     tcolors ={colors{1};colors{2};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',colors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.0001);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
    xticks = xdata; xticklabels = {'C1','C2','C3','C4','C3''','C1''','C2'''};
    xtickangle(30)
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.1 0.05 -0.08 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{sprintf('Avg. %cF/Fo (%%)',916)},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Ca_Signal_conditions'),600);
return;
end
%%
if 1
    [within,dvn,xlabels] = make_within_table({'Conditions'},7);
    dataT = make_between_table({out_CC.mean_msr},dvn);
    ra = repeatedMeasuresAnova(dataT,within);
    [xdata,mVar,semVar,combs,p,h,colors,hollowsep] = get_vals_for_bar_graph(mData,ra,0,[1 1 1]);
    
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 9 1.5 1],'color','w');

    hold on;
%     tcolors ={colors{1};colors{2};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',colors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.0001);
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
    xticks = xdata; xticklabels = {'C1','C2','C3','C4','C3''','C1''','C2'''};;
    xtickangle(30)
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.1 0.05 -0.08 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Avg. FR (AU)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Firing_Rate_conditions'),600);
return;
end

%%
if 1
    
    s = generate_shades(4);
    tcolors = s.g;
    [within,dvn,xlabels] = make_within_table({'Type'},2);
    dataT = make_between_table({out_C.m_sp_animal_level_rest',out_C.m_sp_animal_level_motion'},dvn);
    ra = repeatedMeasuresAnova(dataT,within);

    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
     xdata = [1 1.75]; 
    colors = mData.colors;
        hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 9 1.25 1],'color','w');

    hold on;
%     tcolors ={colors{1};colors{2};colors{1};colors{2};colors{3};colors{4}};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.05,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.4,'sigLinesStartYFactor',0.0001);
    set(gca,'xlim',[0.5 xdata(end)+0.5],'ylim',[0 maxY],'FontSize',6,'FontWeight','Normal','TickDir','out','xcolor','k','ycolor','k');
    xticks = xdata; xticklabels = {'Rest','Motion'};
    hatch(hbs(2),30,'k','-',4,0.1); %angle,color,style,step,width
%     hatch(hbs(1),150,'k','-',2,0.1); %angle,color,style,step,width
    xtickangle(30)
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    changePosition(gca,[0.17 0.05 -0.4 -0.1]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Avg. FR (AU)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Firing_Rate_Motion_vs_Rest'),600);
return;
end

%%
if 1
    tcolors = {'k','r'};
   distD(:,1) = out_C.m_sp_animal_motion;
   distD(:,2) = out_C.m_sp_animal_rest;
   [distDo,allVals,allValsG] = plotDistributions(distD);
   minBin = min(allVals);
   maxBin = max(allVals);
   [h,p,ks2stat] = kstest2(allValsG{1},allValsG{2});
   [h,p,cd,ks2stat] = ttest2(allValsG{1},allValsG{2});
    tcolors = mData.colors;
   incr = 0.001; %maxBin =
   hf = figure(8);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.85 1],'color','w');
   hold on;
%    [ha,hb,hca] = plotDistributions(distD,'colors',tcolors,'maxY',maxBin,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
   [ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin);
   set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
   changePosition(gca,[0.09 0.13 -0.05 -0.13]);
    put_axes_labels(gca,{'Average Firing Rate (AU)',[0 0 0]},{{'Percentage','of Neurons'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution_firing_rate'),600);
end
