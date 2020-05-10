function figure_place_field_dynamics_contexts(fn,allRs,ccs)

data = evalin('base','data');
mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;
n = 0;
%% cell selection

% cells that remained disrupted after context 1
selCells = getCellPop('C_3_4_O');

% selCells = getCellPop('C_1_2_P_3_P_4_P');

% selCells = getCellPop('C_2_3_P_4_P');
% 
% selCells = getCellPop('C_2_3_P_4_P');

% selCells = getCellPop('all')

% selCells = selectCells(data,mData,'Common',[1 2]);


%%
trials = 3:10;
ii = 3;

dataForAlignment = data{ii};
[P,C,cellNums] = findPopulationVectorPlot(dataForAlignment,selCells,trials);
allPeakPos = []; allP = []; allPR = []; allC = [];
for sii = 1:4
%     dataForAlignment = data{sii};
%     [P,C,cellNums] = findPopulationVectorPlot(dataForAlignment,selCells,trials);

    allRs = data(sii);
    Rs = allRs{1};
    [P,C] = findPopulationVectorPlot(Rs,selCells,trials,cellNums);
    allP{sii} = P;
    allC{sii} = C;
    minC(sii) = min(C(:));
    [~,temp] = max(P,[],2);
    allPeakPos(:,sii) = Rs.dist(temp);
    temp = findMeanRasters(Rs,trials);
    temp = temp(selCells,:);
    allPR{sii} = temp(cellNums,:);
end
n = 0;
%%
ff = makeFigureRowsCols(105,[25 0.5 4 1],'RowsCols',[2 4],...
    'spaceRowsCols',[-0.06 -0.065],'rightUpShifts',[0.1 0.04],'widthHeightAdjustment',...
    [30 -15]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[22 6 3.5 2]);
FS = 5;
for sii = 1:4
    axes(ff.h_axes(1,sii));changePosition(gca,[0 0.05 -0.091 -0.1]);
    imagesc(allP{sii});
    axis off
    if sii == 1
        text(-5,1,'1','FontSize',FS,'FontWeight','Normal');
        if size(P,1) < 100
            text(-7,size(P,1),sprintf('%d',size(P,1)),'FontSize',FS,'FontWeight','Normal');
        else
            text(-10,size(P,1),sprintf('%d',size(P,1)),'FontSize',FS,'FontWeight','Normal');
        end
        text(-21,25,sprintf('Cells'),'FontSize',FS+3,'FontWeight','Bold','rotation',90);
    end
    text(3,size(P,1)+round(size(P,1)/10),sprintf('Context %d',sii),'FontSize',FS,'FontWeight','Normal');
    set(gca,'Ydir','Normal','linewidth',1.5,'FontSize',FS,'FontWeight','Bold','YTick',[1 length(selCells)]);
    cols = size(Rs.rasters(:,:,1),2);
    colsHalf = ceil(cols/2);
    ts = round(Rs.dist);
    set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    set(gca,'XTick',[]);
    axes(ff.h_axes(2,sii));
    dec = -0.09;
    changePosition(gca,[0.0 0.05 dec dec]);
    imagesc(allC{sii},[-1 1]);
    box off;
    axis equal
    axis off
    if sii == 1        
        text(-8,3,'0','FontSize',FS,'FontWeight','Normal');
        text(-10,50,'142','FontSize',FS,'FontWeight','Normal');
    end
    text(-1,-3,'0','FontSize',FS,'FontWeight','Normal');
    text(44,-3,'142','FontSize',FS,'FontWeight','Normal');
    if sii == 2
        text(35,-13,sprintf('Position (cm)'),'FontSize',FS+3,'FontWeight','Bold','rotation',0);
    end
    if sii == 1
        text(-21,-3,sprintf('Position (cm)'),'FontSize',FS+3,'FontWeight','Bold','rotation',90);
    end
    set(gca,'Ydir','Normal','linewidth',1,'FontSize',FS,'FontWeight','Bold');
    h = xlabel('Position (cm)');
    changePosition(h,[0 0 0]);
    h = ylabel('Position (cm)');
    changePosition(h,[1 0 0]);
    cols = size(Rs.rasters(:,:,1),2);
    colsHalf = ceil(cols/2);
    ts = round(Rs.dist);
    set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1) ts(colsHalf) ts(cols)],'Ydir','Normal');
end

colormap parula
mI = min(minC);
for ii = 1:4
    axes(ff.h_axes(2,ii));
    caxis([mI 1]);
end

axes(ff.h_axes(2,4));
pos = get(gca,'Position');
h = axes('Position',pos);
changePosition(h,[0.14 0.06 0.0 -0.1]);
hc = colorbar; axis off
set(hc,'YTick',[],'linewidth',0.01);
changePosition(hc,[0 0 0 -0.01]);
ylims = ylim;
xlims = xlim;
text(xlims(1)+0.33,ylims(1)-0.05,sprintf('%.2f',mI),'FontSize',4.5);
text(xlims(1)+0.39,ylims(2)+0.045,'1','FontSize',4.5);

axes(ff.h_axes(1,4));
pos = get(gca,'Position');
h = axes('Position',pos);
changePosition(h,[0.14 0.06 0.0 -0.1]);
hc = colorbar; axis off
set(hc,'YTick',[],'linewidth',0.01);
changePosition(hc,[0 0 0 -0.01]);
ylims = ylim;
xlims = xlim;
text(xlims(1)+0.4,ylims(1)-0.05,sprintf('0'),'FontSize',4.5);
text(xlims(1)+0.39,ylims(2)+0.045,'1','FontSize',4.5);

save2pdf('figure_dynamics_contexts.pdf',ff.hf,600);
sigR = significanceTesting(allP);
return;
%% plot anova
figure(101);clf;
set(gcf,'units','inches');
set(gcf,'Position',[22 4 1.25 2]);hs = sigR.annova.multcompare.h; ps = sigR.annova.multcompare.p;
plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[hs ps],'colors',colors,'sigColor',sigColor,'maxY',0.33,'ySpacing',0.015,'sigTestName','ANOVA','BaseValue',0.001);
xlim([0.4 0.6+length(sigR.means)]);
% hx = xlabel('Trials');changePosition(hx,[0 0.0 0]);
hy = ylabel('Mean Normalized FR');
% ylim([0 0.27]);
xlim([0.5 4.5]);
set(gca,'FontSize',7,'XTick',[1 2 3 4],'XTickLabel',{'Context 1','Context 2','Context 3','Context 4'},'TickDir','out');
xtickangle(30);
save2pdf('anova_dynamics.pdf',gcf,600);

%% peaks shifts

ff = makeFigureWindow__one_axes_only(106,[6 4 5 2.5],[0.19 0.2 0.76 0.75]);
set(gcf,'color','w');
set(gcf,'Position',[15 6 2 2]);
axes(ff.ha);hold on;
distD = [];
distD{1} = allPeakPos(:,2)'-allPeakPos(:,1)';
distD{2} = allPeakPos(:,3)'-allPeakPos(:,2)';
distD{3} = allPeakPos(:,4)'-allPeakPos(:,3)';
legs = {'Context 2 - Context 1','Context 3 - Context 2','Context 4 - Context 3',[-100 0.03 80 5]};
[hb,hca] = plotDistributions(distD,'min',-140,'max',140,'incr',20,'colors',colors,'maxY',80,'cumPos',[0.5 0.35 0.3 0.3],'legend',legs);
axes(ff.ha);
h = xlabel('Difference Peak Position (cm)');
changePosition(h,[0 0.5 0]);
h = ylabel('Percentage');
changePosition(h,[0.1 0 0]);
axes(hca);
save2pdf('figure_3_peak_shifts_dist_contexts.pdf',ff.hf,600);

