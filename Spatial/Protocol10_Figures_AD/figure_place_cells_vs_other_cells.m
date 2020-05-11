function figure_place_cells_vs_other_cells(fn,allRs,ccs)

adata = evalin('base','dataa');
mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;
selAnimals = 1:3;
selAnimals = 1:4;
mData.belt_length = adata{selAnimals(1)}{1}{1}.belt_length;
n = 0;

%%
trials = 3:10;
trials10 = 3:9;
for ii = 1:4%length(data)
    distDi = [];
    mRsi = [];
    for jj = 1:length(selAnimals)
        [ii jj selAnimals(jj)]
         [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'placeCells5',ii);
         distDi = [distDi tempD];
         [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'rasters',ii);
         try
             mR = findMeanRasters(tempD,trials);
         catch
             mR = findMeanRasters(tempD,trials10);
         end
         mRsi = [mRsi;mR];
    end
    distD{ii} = distDi;
    allRs{ii} = mRsi;
    [allP{ii},allC{ii}] = findPopulationVectorPlot(mRsi,find(distDi));
end

% allP = []; allPSig = []; allC = [];
% for sii = 1:length(data)
%     allRs = data(sii);
%     Rs = allRs{1};
%     [P,C] = findPopulationVectorPlot(Rs,selCells{sii},trials,cellNums{sii});
% %     [P,C] = findPopulationVectorPlot(Rs,oSelCells{sii},trials,oCellNums{sii});
%     allP{sii} = P;
%     P(P == 0) = [];
%     
%     allC{sii} = C;
%     minC(sii) = min(C(:));
% %     [~,temp] = max(P,[],2);
% %     allPeakPos(:,sii) = Rs.dist(temp);
% %     temp = findMeanRasters(Rs,trials);
% %     temp = temp(selCells,:);
% %     allPR(:,:,sii) = temp(cellNums,:);
% end
n = 0;
%%
ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 4],...
    'spaceRowsCols',[-0.06 -0.05],'rightUpShifts',[0.1 0.04],'widthHeightAdjustment',...
    [15 -15]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[1 6 3.5 2]);
FS = 5;
for sii = 1:4
    P = allP{sii};
    axes(ff.h_axes(1,sii));changePosition(gca,[0 0.05 -0.091 -0.1]);
    imagesc(P);
    axis off
%     if sii == 1
        text(-5,1,'1','FontSize',FS,'FontWeight','Normal');
        if size(P,1) < 100
            text(-7,size(P,1),sprintf('%d',size(P,1)),'FontSize',FS,'FontWeight','Normal');
        else
            text(-10,size(P,1),sprintf('%d',size(P,1)),'FontSize',FS,'FontWeight','Normal');
        end
        if sii == 1
            text(-21,25,sprintf('Cells'),'FontSize',FS+3,'FontWeight','Bold','rotation',90);
        end
%     end
    text(3,size(P,1)+round(size(P,1)/10),sprintf('Context %d',sii),'FontSize',FS,'FontWeight','Normal');
    set(gca,'Ydir','Normal','linewidth',1.5,'FontSize',FS,'FontWeight','Bold','YTick',[1 sum(distD{sii})]);
    cols = 50;
    colsHalf = ceil(cols/2);
    ts = round(adata{selAnimals(1)}{1}{1}.dist);
    set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    set(gca,'XTick',[]);
    axes(ff.h_axes(2,sii));
    dec = -0.09;
    changePosition(gca,[0.0 0.05 dec dec]);
    imagesc(allC{sii},[-1 1]);
    minC(sii) = min(allC{sii}(:));
    box off;
    axis equal
    axis off
    if sii == 1        
        text(-8,3,'0','FontSize',FS,'FontWeight','Normal');
        text(-10,50,num2str(mData.belt_length),'FontSize',FS,'FontWeight','Normal');
    end
    text(-1,-3,'0','FontSize',FS,'FontWeight','Normal');
    text(44,-3,num2str(mData.belt_length),'FontSize',FS,'FontWeight','Normal');
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
    cols = 50;%size(Rs.rasters(:,:,1),2);
    colsHalf = ceil(cols/2);
%     ts = round(Rs.dist);
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

save_pdf(ff.hf,mData.pdf_folder,'figure_place_cells_py.pdf',600);

return;
%%
% bar graph from anova
sigR = significanceTesting(allP);
ff = makeFigureWindow__one_axes_only(3,[2 4 2 2],[0.3 0.22 0.68 0.72]);
set(gcf,'color','w');
set(gcf,'Position',[3 4 1.2 1.5]);
axes(ff.ha); hs = sigR.anova.multcompare.h; ps = sigR.anova.multcompare.p;
plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[hs ps],'colors',colors,'sigColor',sigColor,'maxY',0.3,'ySpacing',0.03,'sigTestName','ANOVA');
% plotBarsWithSigLines(sigR.means,sigR.sems,sigR.combs,[hs ps],'colors',colors,'sigColor',sigColor,'maxY',0.57,'ySpacing',0.04,'sigTestName','ANOVA');
xlim([0.4 0.6+length(sigR.means)]);
%     ylim([0 max(mVals)]);
hyl = ylabel('Average Normalized FR');
pos = get(hyl,'Position');pos = pos + [+0.4 0 0];set(hyl,'Position',pos);
% set(ff.ha,'linewidth',1);
set(ff.ha,'TickDir','out','FontSize',7,'FontWeight','Normal');
set(ff.ha,'XTick',[1 2 3 4],'XTickLabel',{'Context 1','Context 2','Context 3','Context4'});
xtickangle(25);
save2pdf('figure_place_cells_bars.pdf',ff.hf,600);