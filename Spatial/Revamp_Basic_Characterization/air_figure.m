function light_figure

mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei_11_15 = evalin('base','ei'); 

selContexts = [2 7];
rasterNames = {'air55T','air55T'};
Rs = get_rasters_data(ei_11_15,selContexts,rasterNames);
% Rs = get_rasters_data(ei_2_3,selContexts,rasterNames);
Rs = find_responsive_rasters(Rs,1:10);
mR = calc_mean_rasters(Rs,1:10);
[resp_fractionC,resp_valsC,OIC,mean_OIC,resp_ORC,resp_OR_fractionC,resp_ANDC,resp_AND_fractionC,resp_exc_inh] = get_responsive_fraction(Rs);
[CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,1,0);

[respE,respE_OR,respE_AND,respE_fraction] = get_cell_list_exc_inh(resp_exc_inh,1,1);
[CRcE,aCRcE,mRRE] = find_population_vector_corr(Rs,mR,respE,0);

[respI,respI_OR,respI_AND,respI_fraction] = get_cell_list_exc_inh(resp_exc_inh,1,0);
[CRcI,aCRcI,mRRI] = find_population_vector_corr(Rs,mR,respI,0);

[all_corr_CA,all_corr_cell_CA,mean_corr_CA,mean_cell_corr_CA,xs_CA,paramA] = find_population_vector_corr_remap(Rs,mR,resp_ORC);

n = 0;

%%
% an = 5; cn = 2;
% plotRasters_simplest(Rs{an,cn})
%%
an = 1; cn = 1;
% plotRasters_simplest(Rs{an,cn})
% find(resp_valsC{an}(:,cn));
ff = makeFigureRowsCols(2020,[0.5 0.5 4 1],'RowsCols',[1 4],...
    'spaceRowsCols',[0.15 0.06],'rightUpShifts',[0.08 0.25],'widthHeightAdjustment',...
    [-75 -475]);
set(gcf,'color','w'); set(gcf,'Position',[10 4 4 1]);
ff = sample_rasters(Rs{an,cn},[140 66 162 181],ff);
save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rasters'),600);
%% population vector and correlation single animal
an = 1;
ff = makeFigureRowsCols(106,[1 0.5 4 1],'RowsCols',[2 2],...
    'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.15 0.1],'widthHeightAdjustment',...
    [-70 -60]);
set(gcf,'color','w');
set(gcf,'Position',[5 5 2.2 2]);
% ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRR(an,:),CRc(an,:),[-0.1 1],[]);
ff = show_population_vector_and_corr(mData,ff,Rs(an,:),mRRE(an,:),CRcE(an,:),[-0.1 1],[]);
save_pdf(ff.hf,mData.pdf_folder,sprintf('air_population_vector_corr.pdf'),600);

%% average correlation of all animals
ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 2],...
    'spaceRowsCols',[0 -0.03],'rightUpShifts',[0.15 0.2],'widthHeightAdjustment',...
    [-70 -220]);
set(gcf,'color','w');
set(gcf,'Position',[5 3 2.2 1]);
ff = show_population_vector_and_corr(mData,ff,Rs(an,:),[],aCRc,[-0.1 1],[]);
save_pdf(ff.hf,mData.pdf_folder,sprintf('air_average_population_vector_corr.pdf'),600);


%% average correlation of all animals
ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[2 2],...
    'spaceRowsCols',[0.13 0.13],'rightUpShifts',[0.13 0.2],'widthHeightAdjustment',...
    [-220 -220]);
set(gcf,'color','w');
set(gcf,'Position',[5 3 2 2]);
ff = show_remapping_corr_plots(mData,ff,mean_corr_CA,xs_CA,[-0.1 1]);
save_pdf(ff.hf,mData.pdf_folder,sprintf('air_rest_remap_corr.pdf'),600);


%% excitatory inhibitory responsive percentages
if 1
    [respE1,respE_OR1,respE_AND1,respE_fraction1] = get_cell_list_exc_inh(resp_exc_inh,1,1);
    [respE2,respE_OR2,respE_AND2,respE_fraction2] = get_cell_list_exc_inh(resp_exc_inh,2,1);
    
    [respI1,respI_OR1,respI_AND1,respI_fraction1] = get_cell_list_exc_inh(resp_exc_inh,1,0);
    [respI2,respI_OR2,respI_AND2,respI_fraction2] = get_cell_list_exc_inh(resp_exc_inh,2,0);
    
    dataT = array2table(100*[respE_fraction1' respI_fraction1' respE_fraction2' respI_fraction2']);
    dataT.Properties.VariableNames = {'C1E','C1I','C2E','C2I'};
    within = array2table([[1 1 2 2]' [1 2 1 2]']);
    within.Properties.VariableNames = {'Cond','EI'}; within.Cond = categorical(within.Cond); within.EI = categorical(within.EI);
    
    ra = repeatedMeasuresAnova(dataT,within,0.05);
    
    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
    xdata = [1 2 3 4];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    tcolors ={colors{1};colors{1};colors{2};colors{2};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    for ii = [2 4]
        set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
    end
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C2-E','C2-I','C2''-E','C2''-I'};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.06 0.03 0.02 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Air responsive Cells (%)'},[0 0 0]});

    save_pdf(hf,mData.pdf_folder,'air_responsive_exc_inh',600);
    
end

%%
if 0
% out_CC = all_out_C.all_corr_an{3};
% out_AA = all_out_A.all_corr_an{1};
    temp = mean_cell_corr_CA(:,:,1);
    mask = ones(size(temp)); mask = triu(mask,1) & ~triu(mask,2);
    for ii = 1:size(mean_cell_corr_CA,3)
        temp = mean_cell_corr_CA(:,:,ii);
        var_C(ii,:) = temp(mask)';
    end
    mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
    combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
    xdata = [1 2];
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.5 1],'color','w');
    hold on;
    tcolors ={colors{1};colors{2};colors{3};colors{1};colors{2};colors{3};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',0.2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.1);
    % plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end)+0; xticklabels = {'C3-C4','C4-C3'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30);
    changePosition(gca,[0.06 0.03 0.02 -0.11]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Correlation'},[0 0 0]});

    save_pdf(hf,mData.pdf_folder,'air_rest_cell corr remap',600);
return;
end

%%
dataT = array2table(resp_fractionC*100);
dataT.Properties.VariableNames = {'A1','A2'};
within = array2table([1 2]');
within.Properties.VariableNames = {'Cond'};
within.Cond = categorical(within.Cond);
ra = repeatedMeasuresAnova(dataT,within,0.05);

%%
mVar = ra.est_marginal_means.Mean;semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p < 0.05;
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
% row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
% row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

xdata = [1 2]; 
colors = mData.colors;
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
hold on;
%     tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4};};
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'ySpacing',3,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0.05);
% hbe = bar(3,mean(100*resp_OR_fractionC));
% hbee = errorbar(3,mean(100*resp_OR_fractionC),std(100*resp_OR_fractionC)/sqrt(5), 'k', 'linestyle', 'none','CapSize',3);
% set(hbe,'FaceColor',colors{3},'EdgeColor',colors{3});
% plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
%     maxY = maxY + 5 + 2;
set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 30],'FontSize',6,'FontWeight','Bold','TickDir','out');
xticks = [xdata(1:end)]; xticklabels = {'C2','C2'''};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30)
changePosition(gca,[0.2 0.03 -0.5 -0.05]);
put_axes_labels(gca,{[],[0 0 0]},{{'Air Responsive','Cells (%)'},[0 0 0]});
save_pdf(hf,mData.pdf_folder,sprintf('percentage_air_responsive'),600);

hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[3 7 1.25 1],'color','w');
hold on;

hbe = bar(1,mean(100*resp_OR_fractionC),'barwidth',0.6);
hbee = errorbar(1,mean(100*resp_OR_fractionC),std(100*resp_OR_fractionC)/sqrt(5), 'k', 'linestyle', 'none','CapSize',3);
set(hbe,'FaceColor',colors{3},'EdgeColor',colors{3});
set(gca,'xlim',[0.25 1.75],'ylim',[0 30],'FontSize',6,'FontWeight','Bold','TickDir','out');
xticks = [xdata(1:end)]; xticklabels = {'C2 or C2'''};
set(gca,'xtick',xticks,'xticklabels',xticklabels);
%     xtickangle(30)
changePosition(gca,[0.2 0.03 -0.55 -0.05]);
put_axes_labels(gca,{[],[0 0 0]},{{'Air Responsive Cells','in any Condition(%)'},[0 0 0]});
save_pdf(hf,mData.pdf_folder,sprintf('percentage_air_responsive_any'),600);
%%
%% overlap RM bar graph
hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
    for ii = 1:5
        C12(ii) = OIC{ii}(1,2);
    end
%%
    [mVar,semVar] = findMeanAndStandardError(C12);
    combs = []; p = 1; h = p < 0.05;
%     row = [1 2]; ii = ismember(combs,row,'rows'); p(ii) = mcGroup{1,6}; h(ii) = 1; 
    % row = [5 6]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{2,6}; h(ii) = 1; 
    % row = [3 4]; ii = ismember(combs,row,'rows'); p(ii) = mcTI{3,6}; h(ii) = 1; 

    xdata = [1]; 
    colors = mData.colors;
    hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 1.25 1],'color','w');
    hold on;
%     tcolors = {colors{1};colors{1};colors{2};colors{2};colors{3};colors{3};colors{4};colors{4}};
    tcolors = {colors{1};colors{2};colors{3};colors{4};colors{1};colors{2};colors{3};colors{4};};
    [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
        'ySpacing',2,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.001,...
        'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.5,'sigLinesStartYFactor',0.05);

    set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[0 maxY],'FontSize',6,'FontWeight','Bold','TickDir','out');
    xticks = xdata(1:end); xticklabels = {'C2-C2'''};
    set(gca,'xtick',xticks,'xticklabels',xticklabels);
    xtickangle(30)
    changePosition(gca,[0.3 0.03 -0.5 -0.05]);
    put_axes_labels(gca,{[],[0 0 0]},{{'Overlap Index'},[0 0 0]});
    save_pdf(hf,mData.pdf_folder,sprintf('air_overlap_stats'),600);