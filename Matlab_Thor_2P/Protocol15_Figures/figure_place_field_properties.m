function figure_place_field_properties(fn,allRs,ccs)
adata = evalin('base','data15');
mData = evalin('base','mData15');
colors = mData.colors;
sigColor = mData.sigColor;
% selAnimals = 1:11;
selAnimals = [1 2 3 4 5];
% selAnimals = 10;
mData.belt_length = adata{selAnimals(1)}{1}{1}.belt_length;
n = 0;

%%
trials = 3:10;
trials10 = 3:9;
for ii = 1:3%length(data)
    distDi = []; pwsi = []; pcsi = []; rsi = []; sii = [];
    for jj = 1:length(selAnimals)
        [ii jj selAnimals(jj)]
        [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'placeCells5',ii);
        distDi = [distDi tempD];
        [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'pws',ii);
        pwsi = [pwsi tempD];
        [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'centers',ii);
        pcsi = [pcsi tempD];
        [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'rs',ii);
        rsi = [rsi tempD];
        [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'SI',ii);
        sii = [sii tempD];
    end
    distD{ii} = distDi; pw{ii} = pwsi; pc{ii} = pcsi; Rsq{ii} = rsi; SI{ii} = sii;
end
pcs5 = distD; distD = [];
% varNames = {'pw','pc','Rsq','SI'};
% for ii = 1:length(varNames)
%     for jj = 1:4
%         cmdTxt = sprintf('tempVals = %s{jj};',varNames{ii}); eval(cmdTxt);
%         tempVals = tempVals(find(distD{jj}));
%         cmdTxt = sprintf('%s{jj} = tempVals;',varNames{ii}); eval(cmdTxt);
%     end
% end

%% Percentage of place cells vs place field centers
runThis = 1;
if runThis
distD = [];
for ii = 1:length(pc)
    tempV = pc{ii}(find(pcs5{ii}));
    distD{ii} = tempV(~isnan(tempV));
end    
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 7.5 6 2.5],'color','w');hold on;
[ha,hb,hca,sigR] = plotDistributions(distD,'colors',colors,'maxY',30,'cumPos',[0.6 0.4 0.23 0.35],'min',0,'incr',10,'max',mData.belt_length,'BaseValue',0.15);
set(gca,'TickDir','out','FontSize',mData.axes_font_size,'FontWeight','Bold');
hx = xlabel('Place Field Center (cm)');changePosition(hx,[0 0.0 0]);
hy = ylabel('Percentage (%)');changePosition(hy,[1 -1 0]);
legs = {'Context 1','Context 2','Context 3',[60 5 30 3]};
% sigR = significanceTesting(distD);
putLegend(gca,legs,'colors',colors,'sigR',{sigR,'ks',sigColor,10});
changePosition(gca,[-0.04 0.09 0.11 -0.06]);
axes(hca);set(gca,'FontSize',10);
save_pdf(hf,mData.pdf_folder,'figure_percent_place_cells_vs_centers_15.pdf',600);
return;
end

%%

%% find difference place field widths 1st half beld vs 2nd half belt
runThis = 1;
if runThis
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 6 2.5],'color','w');hold on;
ya = []; xa = []; gr = [];
ind = 1;
for ii = 1:length(pw)
    pwidths = pw{ii}(find(pcs5{ii}))';
    pcenters = pc{ii}(find(pcs5{ii}))';
    indsL = pcenters < (mData.belt_length/2);
    indsG = pcenters > (mData.belt_length/2);
    pwL = pcenters(indsL);
    pwG = pcenters(indsG);
    pwData{ind} = pwL;ind = ind + 1;
    pwData{ind} = pwG;ind = ind + 1;
end
sigR = significanceTesting(pwData);
hs = sigR.anova.multcompare.h; ps = sigR.anova.multcompare.p;
plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[hs ps],'colors',{colors{1},colors{1},colors{2},colors{2},colors{3},colors{3}},'sigColor',sigColor,...
    'maxY',250,'ySpacing',12,'sigTestName','ANOVA','sigLineWidth',0.25,'BaseValue',0.001,...
    'sigFontSize',9,'sigAsteriskFontSize',15);
xlabel('Contexts - First Half (FH), Second Half (SH)'); hy = ylabel('Mean Place Field Width (cm)'); changePosition(hy,[-0.1 -30 0]);
% legs = {'Context 1','Context 2','Context 3',[11 0.25 0.27 0.025]};
% putLegend(gca,legs,'colors',colors,'sigR',{[],'ks',sigColor,10});
set(gca,'XTick',[1:6],'XTickLabel',{'FH','SH','FH','SH','FH','SH'});
set(gca,'FontSize',mData.axes_font_size,'FontWeight','Bold','TickDir','out');changePosition(gca,[-0.01 0.01 0.06 -0.01]);
% text(1,0.29,'Light Responsive and Place Cells','FontSize',mData.axes_font_size,'FontWeight','Bold');
save_pdf(hf,mData.pdf_folder,'place_field_widths_1st_2nd_half_15',600);
return;
end


%% scatter place field widths vs centers
runThis = 0;
if runThis
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 6 2.5],'color','w');hold on;
ya = []; xa = []; gr = [];
for ii = 1:length(pw)
    pwidths = pw{ii}(find(pcs5{ii}))';
    pcenters = pc{ii}(find(pcs5{ii}))';
    scatter(pcenters,pwidths,20,'.','MarkerEdgeColor',colors{ii},'MarkerFaceColor','none');
    ft = fittype('poly1');
    [ftc,gof,output] = fit(pcenters,pwidths,ft);    rsq(ii) = gof.rsquare;
    co = coeffvalues(ftc);
    pwf = ft(co(1),co(2),pcenters);
    hold on;plot(pcenters,pwf,'color',colors{ii},'LineWidth',1)
    xa = [xa;pcenters]; ya = [ya;pwidths]; gr = [gr;ones(size(pcenters))*ii];
end
sigR = do_ancova(xa,ya,gr);
ylim([0 100]);
set(gca,'TickDir','out','FontSize',mData.axes_font_size,'FontWeight','Bold');
hx = xlabel('Place Field Center (cm)');changePosition(hx,[0 0.0 0]);
hy = ylabel('Place Field Width (cm)');changePosition(hy,[2.2 -3.0 0]);
for ii = 1:length(pw)
    legs{ii} = sprintf('C %d',ii);%{'Context 1','Context 2','Context 3','Context 4',[30 3 35 3]};
end
legs{ii+1} = [35 5 93 7];putLegend(gca,legs,'colors',colors,'sigR',{sigR,'ancovas',sigColor,9}); text(30,100,'Slope','FontWeight','Bold');

changePosition(gca,[-0.03 0.09 0.09 -0.06]);
save_pdf(hf,mData.pdf_folder,'figure_scatter_centers_widths_15.pdf',600);
return;
end

%% scatter zMI vs centers
runThis = 0;
if runThis
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 6 2.5],'color','w');hold on;
ya = []; xa = []; gr = [];
for ii = 1:length(pw)
    SIs = SI{ii}(find(pcs5{ii}))';
    pcenters = pc{ii}(find(pcs5{ii}))';
    scatter(pcenters,SIs,20,'.','MarkerEdgeColor',colors{ii},'MarkerFaceColor','none');
    ft = fittype('poly1');
    [ftc,gof,output] = fit(pcenters,SIs,ft);  rsq(ii) = gof.rsquare;
    co = coeffvalues(ftc);
    pwf = ft(co(1),co(2),pcenters);
    hold on;plot(pcenters,pwf,'color',colors{ii},'LineWidth',1)
    xa = [xa;pcenters]; ya = [ya;SIs]; gr = [gr;ones(size(pcenters))*ii];
end
sigR = do_ancova(xa,ya,gr);
ylim([0 35]);
set(gca,'TickDir','out','FontSize',mData.axes_font_size,'FontWeight','Bold');
hx = xlabel('Place Field Center (cm)');changePosition(hx,[0 0.0 0]);
hy = ylabel('Mutual Info. (Z Score)');changePosition(hy,[0 0.0 0]);
for ii = 1:length(pw)
    legs{ii} = sprintf('C %d',ii);%{'Context 1','Context 2','Context 3','Context 4',[30 3 35 3]};
end
legs{ii+1} = [130 5 33 2];putLegend(gca,legs,'colors',colors,'sigR',{sigR,'ancovai',sigColor,9}); text(125,36,'Intercept','FontWeight','Bold'); 
legs{end} = [70 5 33 2];putLegend(gca,legs,'colors',colors,'sigR',{sigR,'ancovas',sigColor,9}); text(65,36,'Slope','FontWeight','Bold');

changePosition(gca,[-0.03 0.09 0.09 -0.06]);
save_pdf(hf,mData.pdf_folder,'figure_scatter_centers_zMI_15.pdf',600);
return;
end

%% place width distribution vs place field centers
runThis = 0;
if runThis
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[13 4 6 2.5],'color','w');hold on;
bins = 0:10:mData.belt_length;
xs = 5:10:mData.belt_length;
sigD = [];
for ii = 1:length(pw)
    pwidths = pw{ii}(find(pcs5{ii}));
    pcenters = pc{ii}(find(pcs5{ii}));
    mpws = []; sempws = [];
    for jj = 1:(length(bins)-1)
        st = bins(jj);
        se = bins(jj+1);
        inds = find(pcenters>st & pcenters<=se);
        [mpws(jj) sempws(jj)] = findMeanAndStandardError(pwidths(inds));
        sigD{ii,jj} = pwidths(inds);
    end
    inds = isnan(mpws);
    mpws(inds) = []; sempws(inds) = [];
    errorbar(xs(1:length(mpws)),mpws,sempws,'color',colors{ii},'linewidth',1,'capsize',3)
end
xlim([0 140]);
set(gca,'TickDir','out','FontSize',mData.axes_font_size,'FontWeight','Bold');
hx = xlabel('Place Field Center (cm)');changePosition(hx,[0 0.0 0]);
hy = ylabel('Place Field Width (cm)');changePosition(hy,[1 -2 0]);
legs = {'Context 1','Context 2','Context 3','Context 4',[30 5 80 7]};
sigR = significanceTesting(sigD,8);
putLegend(gca,legs,'colors',colors,'sigR',{sigR,'ancova',sigColor,10});
changePosition(gca,[-0.04 0.09 0.11 -0.06]);
save_pdf(hf,mData.pdf_folder,'figure_widths_vs_centers.pdf',600);
return;
end
%% MI distribution vs place field centers
runThis = 0;
if runThis
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 6 2.5],'color','w');hold on;
bins = 0:10:145;
xs = 5:10:145;
sigD = [];
for ii = 1:length(pw)
    pwidths = pw{ii}(find(pcs5{ii}));
    pcenters = pc{ii}(find(pcs5{ii}));
    SIs = SI{ii}(find(pcs5{ii}));
    mpws = []; sempws = [];
    for jj = 1:(length(bins)-1)
        st = bins(jj);
        se = bins(jj+1);
        inds = find(pcenters>st & pcenters<=se);
        [mpws(jj) sempws(jj)] = findMeanAndStandardError(SIs(inds));
        sigD{ii,jj} = SIs(inds);
    end
    inds = isnan(mpws);
    mpws(inds) = []; sempws(inds) = [];
    errorbar(xs(1:length(mpws)),mpws,sempws,'color',colors{ii},'linewidth',1,'capsize',3)
end
xlim([0 140]); ylim([1 30]);
set(gca,'TickDir','out','FontSize',mData.axes_font_size,'FontWeight','Bold');
hx = xlabel('Place Field Center (cm)');changePosition(hx,[0 0.0 0]);
hy = ylabel('Mutual Information (z-score)');changePosition(hy,[2 -3.5 0]);
legs = {'Context 1','Context 2','Context 3','Context 4',[60 5 30 3]};
sigR = significanceTesting(sigD,8);
putLegend(gca,legs,'colors',colors(1:end),'sigR',{sigR,'anova',sigColor,10});
changePosition(gca,[-0.04 0.09 0.11 -0.06]);
save_pdf(hf,mData.pdf_folder,'figure_MIs_vs_centers.pdf',600);
return
end
%% Rsq distribution vs place field centers
runThis = 0;
if runThis
ff = makeFigureWindow__one_axes_only(106,[6 4 5 2.5],[0.15 0.25 0.8 0.7]);
set(gcf,'color','w');
set(gcf,'Position',[22 7 2.5 1.5]);
axes(ff.ha);hold on;

bins = 0:10:145;
xs = 5:10:145;
sigD = [];
for ii = 1:length(pw)
    pwidths = pw{ii};
    Rsqs = Rsq{ii};
    pcenters = pc{ii};
    mpws = []; sempws = [];
    for jj = 1:(length(bins)-1)
        st = bins(jj);
        se = bins(jj+1);
        inds = find(pcenters>st & pcenters<=se);
        [mpws(jj) sempws(jj)] = findMeanAndStandardError(Rsqs(inds));
        sigD{ii,jj} = Rsqs(inds);
    end
    inds = isnan(mpws);
%     mpws(inds) = []; sempws(inds) = [];
    errorbar(xs(1:length(mpws)),mpws,sempws,'color',colors{selCon(ii)},'linewidth',0.5,'capsize',3)
end
xlim([0 140]);
ylim([0.2 1.2]);
set(gca,'TickDir','out','FontSize',7,'FontWeight','Normal');
hx = xlabel('Place Field Center (cm)');changePosition(hx,[0 0.0 0]);
hy = ylabel('Goodness of Fit (R^2)');changePosition(hy,[2 0.0 0]);
legs = {'Context 1','Context 2','Context 3','Context 4',[100 5 1.05 0.05]};
% legs = {'Context 1','Context 2',[100 5 1.05 0.05]};
sigR = significanceTesting(sigD,9);
putLegend(gca,legs,'colors',colors(1:end),'sigR',{sigR,'annova',sigColor,4});
save2pdf('figure_Rsqs_vs_centers.pdf',ff.hf,600);
return;
end

%% difference in centers  vs place field centers
% 
% ff = makeFigureWindow__one_axes_only(106,[6 4 5 2.5],[0.15 0.25 0.8 0.7]);
% set(gcf,'color','w');
% set(gcf,'Position',[22 7 2.5 1.5]);
% axes(ff.ha);hold on;
% 
% bins = 0:10:145;
% xs = 5:10:145;
% sigD = [];
% for ii = 1:length(dpc)
%     dpcs = dpc{ii};
%     pcenters = pc{ii};
%     mpws = []; sempws = [];
%     for jj = 1:(length(bins)-1)
%         st = bins(jj);
%         se = bins(jj+1);
%         inds = find(pcenters>st & pcenters<=se);
%         [mpws(jj) sempws(jj)] = findMeanAndStandardError(dpcs(inds));
%         sigD{ii,jj} = dpcs(inds);
%     end
%     inds = isnan(mpws);
%     mpws(inds) = []; sempws(inds) = [];
%     errorbar(xs(1:length(mpws)),mpws,sempws,'color',colors{selCon(ii+1)},'linewidth',0.5,'capsize',3)
% end
% % xlim([0 140]);
% ylim([-20 30]);
% set(gca,'TickDir','out','FontSize',7,'FontWeight','Normal');
% hx = xlabel('Place Field Center (cm)');changePosition(hx,[0 0.0 0]);
% hy = ylabel('Difference centers (cm)');changePosition(hy,[2 0.0 0]);
% legs = {'Context 1','Context 2','Context 3','Context 4',[100 5 1.05 0.05]};
% legs = {'(Context 1) - (Context 2)','(Context 2) - (Context 3)','(Context 3) - (Context 4)',[20 5 30 4]};
% sigR = significanceTesting(sigD,10);
% putLegend(gca,legs,'colors',colors(2:end),'sigR',{sigR,'annova',sigColor,4});
% save2pdf('figure_dpc_vs_centers.pdf',ff.hf,600);

%% difference place field centers
runThis = 0;
if runThis
distD = [];
for ii = 1:length(dpc)
    distD{ii} = dpc{ii}';
end

ff = makeFigureWindow__one_axes_only(107,[10 7 2.5 1.5],[0.14 0.23 0.85 0.73]);
set(gcf,'color','w');
axes(ff.ha);hold on;
[ha,hb,hca,sigR] = plotDistributions(distD,'colors',colors,'maxY',80,'cumPos',[0.6 0.32 0.23 0.35],'min',-25,'incr',5,'max',25,'BaseValue',0.15);
legs = {'(Context 1) - (Context 2)','(Context 2) - (Context 3)','(Context 3) - (Context 4)',[-15 2 80 5]};
putLegend(ff.ha,legs,'colors',colors,'sigR',{sigR,'ks',sigColor,4});
h = xlabel('Difference in Place Field Centers (cm)');
changePosition(h,[0 0.7 0]);
h = ylabel('Percentage');
changePosition(h,[0.5 0 0]);
axes(hca);
save2pdf('figure_diff_centers_dist.pdf',ff.hf,600);

% bar graph from anova
% ff = makeFigureWindow__one_axes_only(2,[1 4 1.75 2],[0.2 0.35 0.7 0.6]);
ff = makeFigureWindow__one_axes_only(3,[2 4 2 2],[0.3 0.24 0.68 0.72]);
set(gcf,'color','w');
set(gcf,'Position',[22 4 1.2 1.5]);
axes(ff.ha); hs = sigR.annova.multcompare.h; ps = sigR.annova.multcompare.p;
plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[hs ps],'colors',colors,'sigColor',sigColor,'maxY',15,'ySpacing',2,'sigTestName','ANOVA');
xlim([0.4 0.6+length(sigR.means)]);
%     ylim([0 max(mVals)]);
hyl = ylabel('Avg. Diff. PF Center (cm)');
pos = get(hyl,'Position');pos = pos + [+0.1 0 0];set(hyl,'Position',pos);
% set(ff.ha,'linewidth',1);
set(ff.ha,'TickDir','out','FontSize',7,'FontWeight','Normal');
set(ff.ha,'XTick',[1 2 3 4],'XTickLabel',{'Context 1-2','Context 2-3','Context 3-4'});
xtickangle(25);
save2pdf('temp.pdf',ff.hf,600);
pause(0.1);
append_pdfs('figure_diff_centers_dist.pdf','temp.pdf');
return;
end