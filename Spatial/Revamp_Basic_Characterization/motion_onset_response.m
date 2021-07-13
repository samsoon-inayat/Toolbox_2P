function firing_rate_motion_vs_rest

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_C = evalin('base','ei'); 

selContexts = [0,0];
rasterNames = {'motionOnsets','motionOffsets'};

Rs = get_rasters_data(ei_C,selContexts,rasterNames);
Rs = find_responsive_rasters(Rs,[]);
[resp_fractionCT,resp_valsCT,OICT,mean_OICT,resp_ORCT,resp_OR_fractionCT,resp_ANDCT,resp_AND_fractionCT] = get_responsive_fraction(Rs);
mR = calc_mean_rasters(Rs,[]);

n = 0;
%%
if 1
    an = 5; 
    plotRasters_simplest(Rs{an,1},[])
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
    ht = get_obj(ff,'title'); hyl = get_obj(ff,'ylabel'); changePosition(hyl(1,1),[7 0 0]);
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
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 35],'FontSize',6,'FontWeight','Normal','TickDir','out');
    xticks = [xdata(1:end)]; xticklabels = {'C2','C2'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels); xtickangle(30)
    any_mean = mean(100*resp_OR_fractionC);    any_sem = std(100*resp_OR_fractionC)/sqrt(5);
    pmchar=char(177); any_text = sprintf('%.0f%c%.0f%%',any_mean,pmchar,any_sem); %text(0.75,31,any_text,'FontSize',6);
    changePosition(gca,[0.2 0.03 -0.5 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Air Responsive','Cells (%)'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('percentage_motion_responsive'),600);
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
