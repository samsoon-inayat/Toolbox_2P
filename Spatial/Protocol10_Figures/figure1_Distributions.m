function figure1_Distributions

adata = evalin('base','datab');
mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;
axes_font_size = mData.axes_font_size;
% allCells = mData.allCells;
selAnimals = 1:11;
selAnimals = 1:4;
n = 0;
%%
runThis = 0;
if runThis
for ii = 1:4%length(data)
    distDi = [];
    for jj = 1:length(selAnimals)
         distDi = [distDi getVariableValues(adata{selAnimals(jj)},'rs',ii)];
    end
    distD{ii} = distDi;
end
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 6.5 2.5],'color','w');
hold on;
[ha,hb,hca,sigR] = plotDistributions(distD,'colors',colors,'maxY',30,'cumPos',[0.5 0.36 0.25 0.4],'min',0,'incr',0.05,'max',1);
hold on;
legs = [];
for ii = 1:length(distD)
    legs{ii} = sprintf('C%d (N = %d)',ii,length(distD{ii}));
end
legs{ii+1} = [0.3 0.03 30 3.5];
putLegend(ha,legs,'colors',colors,'sigR',{sigR,'anova',sigColor,10});
h = xlabel('Goodness of fit R^2');changePosition(h,[0 -2 0]);h = ylabel('Percentage');changePosition(h,[0.0 0 0]);
set(gca,'FontSize',axes_font_size,'FontWeight','Bold');changePosition(ha,[-0.04 0.1 0.12 -0.06]);
axes(hca);set(gca,'FontSize',9);ylim([20 100]);
save_pdf(hf,mData.pdf_folder,'Distribution Of Rs2',600);
return;
end
%%
runThis = 1;
if runThis
for ii = 1:4%length(data)
    distDi = [];
    for jj = 1:length(selAnimals)
         distDi = [distDi getVariableValues(adata{selAnimals(jj)},'SI',ii)];
    end
    distD{ii} = distDi;
end
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 6.9 2.5],'color','w');
hold on;
[ha,hb,hca,sigR] = plotDistributions(distD,'colors',colors,'maxY',50,'cumPos',[0.5 0.36 0.25 0.3],'min',0,'incr',1,'max',21);
hold on;
legs = [];
for ii = 1:length(distD)
    legs{ii} = sprintf('C%d (N = %d)',ii,length(distD{ii}));
end
legs{ii+1} = [5 0.3 45 4.3];
putLegend(ha,legs,'colors',colors,'sigR',{sigR,'anova',sigColor,8});
axes(ha);h = xlabel('Mutual Information (Z Score)');changePosition(h,[0 -3 0]);
set(gca,'FontSize',axes_font_size,'FontWeight','Bold');
changePosition(h,[0 0.7 0]);h = ylabel('Percentage');changePosition(h,[0.0 0 0]);
axes(hca);set(gca,'FontSize',6);
set(gca,'FontSize',axes_font_size,'FontWeight','Bold');changePosition(ha,[-0.04 0.09 0.13 -0.05]);
save_pdf(hf,mData.pdf_folder,'Distribution Of zMI',600);
return;
end
%%
annovaVar = [];
for ii = 1:length(data)
    annovaVar = [annovaVar data{ii}.SI(selCells)'];
    names{ii} = data{ii}.name;
    [mVals(ii) semVals(ii)] = findMeanAndStandardError(data{ii}.SI(selCells));
end
[p,tbl,stats] = anova1(annovaVar);%,names);
% [p,tbl,stats] = kruskalwallis(y2wa,subst,'on');
figure(2001);
[c,~,~,gnames] = multcompare(stats,'CType','hsd');
pdf = c(:,6);
hdf = pdf<0.05;
selGroups = 1:length(data);
combs = nchoosek(1:length(selGroups),2);

for ii = 1:size(combs,1)
    inds = ismember(c(:,[1 2]),[selGroups(combs(ii,1)) selGroups(combs(ii,2))],'rows');
    prT(ii,1) = c(inds,6);
end
hrT = prT<0.05;
%%
% bar graph from anova
% ff = makeFigureWindow__one_axes_only(2,[1 4 1.75 2],[0.2 0.35 0.7 0.6]);
ff = makeFigureWindow__one_axes_only(3,[2 4 2 2],[0.25 0.17 0.7 0.79]);
set(gcf,'color','w');
set(gcf,'Position',[25 4 1.5 2]);
axes(ff.ha);
plotBarsWithSigLines(mVals,semVals,combs,[hrT prT],'colors',thisCols_all(selGroups),'ySpacingFactor',10);
xlim([0.4 0.6+length(selGroups)]);
%     ylim([0 max(mVals)]);
hyl = ylabel('Average MI (z-score)');
pos = get(hyl,'Position');pos = pos + [+0.1 0 0];set(hyl,'Position',pos);
set(ff.ha,'linewidth',1);
set(ff.ha,'TickDir','out','FontSize',8,'FontWeight','bold');
set(ff.ha,'XTick',[1 2 3 4],'XTickLabel',{'Context 1','Context 2','Context 3','Context4'});
xtickangle(25);
save2pdf('MI_anova_bars.pdf',ff.hf,600);







% %%
% allRs = data(1);
% ccs = allCells;
% [~,cellNums] = sort(allRs{1}.SI);
% % if exist('sorting','var')
% %     ptc = findMeanRasters(allRs{sorting});
% %     ptc = ptc(ccs,:);
% %     ptc = normalizeSignal(ptc,2);
% %     [~,peakPos] = max(ptc,[],2);
% %     [~,cellNums] = sort(peakPos);
% % end
% 
% numberOfRows = 2;
% numberOfCols = 1;
% graphsOnOneFigure = numberOfRows * numberOfCols;
% % numberOfData = length(ccsi);
% % numberOfGroups = ceil(numberOfData/graphsOnOneFigure);
% % numberOfFullGroups = floor(numberOfData/graphsOnOneFigure);
% % indices = NaN(1,(numberOfGroups*graphsOnOneFigure));
% % indices(1:numberOfData) = 1:numberOfData;
% % groupIndices = reshape(indices,numberOfRows,numberOfCols,numberOfGroups);
% % % for gg = 1:numberOfGroups
% % %     groupIndices(:,:,gg) = groupIndices(:,:,gg)';
% % % end
% 
% ff = makeFigureRowsCols(105,[25 0.5 4 1],'RowsCols',[numberOfRows numberOfCols],...
%     'spaceRowsCols',[0.05 0.05],'rightUpShifts',[0.1 0.17],'widthHeightAdjustment',...
%     [-70 -170]);
% gg = 1;
% set(gcf,'color','w');
% set(gcf,'Position',[25 4 1 2]);
% 
% FS = 8;
% 
% rows = numberOfRows;%length(allRs);
% 
% indices = 1:(2*rows);
% indices = reshape(indices,2,rows);
% indices = indices';
% % figure(fn);clf;
% for ii = 1
%     Rs = allRs{ii};
%     if exist('sorting','var')
%         [P,C] = findPopulationVectorPlot(Rs,ccs,cellNums);
%     else
%         [P,C] = findPopulationVectorPlot(Rs,ccs);
%     end
%     axes(ff.h_axes(1,1));
%     changePosition(gca,[0 0.05 -0.1 -0.1]);
%     imagesc(P);
% %     colorbar
% 
%     colormap parula;
%         % caxis([minC maxC]);
%     set(gca,'Ydir','Normal','linewidth',1.5,'FontSize',FS,'FontWeight','Bold');
% 
%     xlabel('Position (cm)');
%     ylabel('Cell Number');
%     cols = size(Rs.rasters(:,:,1),2);
%     colsHalf = ceil(cols/2);
%     ts = round(Rs.dist);
%     set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
%     h = xlabel('Position (cm)');
%     changePosition(h,[0 0 0]);
% 
%     axes(ff.h_axes(2,1));
%     dec = -0.09;
%     changePosition(gca,[0.0 0.05 dec dec]);
%     imagesc(C,[-1 1]);
% %     xlim([0 142]);
% %     ylim([0 142]);
%     box off;
% %     hc = colorbar
% %     changePosition(hc,[0.0 0.05 0 0.33]);
%     colormap parula;
%     % caxis([minC maxC]);
%     axis equal
%     set(gca,'Ydir','Normal','linewidth',1,'FontSize',FS,'FontWeight','Bold');
%     h = xlabel('Position (cm)');
%     changePosition(h,[0 0 0]);
%     h = ylabel('Position (cm)');
%     changePosition(h,[1 0 0]);
%     cols = size(Rs.rasters(:,:,1),2);
%     colsHalf = ceil(cols/2);
%     ts = round(Rs.dist);
%     set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
%     set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1) ts(colsHalf) ts(cols)],'Ydir','Normal');
% 
% end
% save2pdf('figure_1_pop_vec.pdf',ff.hf,600);
% 
% function [ptc,CRc] = findPopulationVectorPlot(Rs,ccs,cellNums)
% ptc = findMeanRasters(Rs);
% ptc = ptc(ccs,:);
% ptc = normalizeSignal(ptc,2);
% if ~exist('cellNums','var')
%     [~,peakPos] = max(ptc,[],2);
%     [~,cellNums] = sort(peakPos);
%     ptc = ptc(cellNums,:);
% else
%     ptc = ptc(cellNums,:);
% end
% CRc = corrcoef(ptc);
% 
% 
