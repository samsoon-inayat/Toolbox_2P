function motion_cells
% close all
thisCols_all = {[0 0 0],'b','r',[0 0.7 0.3],[0 0.7 0.3],'r',[0 0.7 0.3],[0 0 0],'b','m','r'};

ei = evalin('base','ei');
owsi = 0;
tei = ei([1 2]);

fileName = makeName('trialsRastersMetaData.mat',tei{1}.folders.thispFolder);
trialsRastersMD = load(fileName);
trials = trialsRastersMD.trials;
% trials = {[],[1:15],[16:31],[32:47]};
trl1 = trials{1}; trl2 = trials{2}; trl3 = trials{3}; trl4 = trials{4};
% trl1 = 1:10; trl3 = 11:26; trl4 = 27:42;
% nbb_motion = getRasters_motion_onsets(tei,trl1);
% temp = find_mutual_information(tei,nbb_air,owsi); nbb_air.SI = temp.zMI;
% nbb_belt = getRasters_belt(tei,trl1);
% temp = find_mutual_information(tei,nbb_belt,owsi); nbb_belt.SI = temp.zMI;
nbc_motion = getRasters_motion_onsets(tei,trl2);
temp = find_mutual_information(tei,nbc_motion,owsi); nbc_motion.SI = temp.zMI;

selCells{1} = 1:size(nbc_motion.rasters,3);
% zScoreTh = 5;
% selCells{1} = find(nbc_motion.SI > zScoreTh);% & nbb_belt.SI < zScoreTh & bc_air_belt.SI > zScoreTh);
n = 0;
% plotRastersMulti3(nbb_belt,bc_air_belt,bb_air_belt,nbb_air,1:size(nbb_air.rasters,3),0,0);
plotRastersMulti({nbc_motion},selCells{1},0,0);
% % plotRastersMeanMulti(nbb_air,nbb_belt,bc_air_belt,bb_air_belt,selCells,0,0);
% return;



% ff = makeFigureRowsCols(10,[1 1 6 2.5],'RowsCols',[1 4],'spaceRowsCols',[-0.009 0.0009],'rightUpShifts',[0.04 0.01],'widthHeightAdjustment',[0 0]);
rows = 4;cols = 2;
totalPlots = rows * cols;
plotNums = 1:totalPlots;
plotInds = reshape(plotNums,cols,rows)';
varNames = {'nbb_air','nbb_belt','bc_air_belt','bb_air_belt'};
cmdText = sprintf('ptc = findMeanRasters(%s);',varNames{3});
eval(cmdText);
ccs = find(bc_air_belt.SI>7);
% ccs = find(bb_air_belt.SI>5);
[~,peakPos] = max(ptc(ccs,:),[],2);[~,cellNums] = sort(peakPos);
hf = figure(20);clf;
for ii = 1:4
    cmdText = sprintf('ptc = findMeanRasters(%s);',varNames{ii});
    eval(cmdText);
    ptc = ptc(ccs,:);
%     normSel = 1;
    [~,peakPos] = max(ptc,[],2);[~,cellNums] = sort(peakPos);
    ptc = ptc(cellNums,:);
    ptc = normalizeSignal(ptc,2);
    CRc = corrcoef(ptc);
    subplot(4,2,plotInds(ii,1));
    imagesc(ptc);colorbar
    colormap jet;
    % caxis([minC maxC]);
    set(gca,'Ydir','Normal','linewidth',1.5,'FontSize',12,'FontWeight','Bold');
    xlabel('Position (cm)');ylabel('Cells');
    set(gca,'XTick',[1 25 50],'XTickLabel',[0 71 142]);
%     text(73,100,'Normalized Spike Rate','Rotation',90,'FontWeight','Bold');

    subplot(4,2,plotInds(ii,2));
    imagesc(CRc);colorbar
    colormap jet;
    % caxis([minC maxC]);
    set(gca,'Ydir','Normal','linewidth',1.5,'FontSize',12,'FontWeight','Bold');
    axis equal
    set(gca,'XTick',[1 25 50],'YTick',[1 25 50],'XTickLabel',[0 71 142],'YTickLabel',[0 71 142]);
    xlabel('Position (cm)');ylabel('Position (cm)');
%     text(76,-5,'Pearson Correlation Coefficient','Rotation',90,'FontWeight','Bold');
    n = 0;
end
n=0;


figure(300);clf;
SIs = [nbb_air.SI;nbb_belt.SI;bc_air_belt.SI;bb_air_belt.SI];
imagesc(SIs);
colorbar

% for distribution of place field centers
maxSI = 50;
minSI = 1;
bins = minSI:1:maxSI;
[bar1 xs] = hist(nbb_air.pcc(~isnan(nbb_air.pcc)),bins); bar1 = 100*bar1/sum(bar1);
[bar2 xs] = hist(nbb_belt.pcc(~isnan(nbb_belt.pcc)),bins); bar2 = 100*bar2/sum(bar2);
[bar3 xs] = hist(bc_air_belt.pcc(~isnan(bc_air_belt.pcc)),bins); bar3 = 100*bar3/sum(bar3);
[bar4 xs] = hist(bb_air_belt.pcc(~isnan(bb_air_belt.pcc)),bins); bar4 = 100*bar4/sum(bar4);
% allBars = [100*bar1/sum(bar1);100*bar2/sum(bar2);100*bar3/sum(bar3);100*bar4/sum(bar4)];
allBars = [bar1;bar2;bar3;bar4];

ff = makeFigureWindow__one_axes_only(6,[6 4 7 2.5],[0.13 0.27 0.85 0.7]);
axes(ff.ha);hold on;
hb = bar(xs,allBars');
for ii = 1:length(hb)
    set(hb(ii),'facecolor',thisCols_all{ii});
end
xlim([bins(1) bins(end)]);
ylim([0 15]);
set(gca,'TickDir','out','FontSize',14,'FontWeight','Bold');
xlabel('Mutual Information (z-score)');
ylabel('Percentage');

x1 = 5; x2 = x1+3; y1 = (15:-1.5:0); y1 = y1(1:4); y2 = y1;
legendFontSize = 11;
legs = {'nbb_air','nbb_belt','bc_air_belt','bb_air_belt'};
for ii = 1:length(legs)
    plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols_all{ii},'linewidth',2);
    text(x2+0.5,y1(ii),sprintf('%s',legs{ii}),'Color',thisCols_all{ii},'FontSize',legendFontSize);
end

cs1 = cumsum(bar1);cs2 = cumsum(bar2);
cs3 = cumsum(bar3);cs4 = cumsum(bar4);
axesPos = ff.pos + [0.6 0.35 0 0];
axesPos(3:4) = [0.25 0.3];
axes('Position',axesPos);hold on;
plot(xs,cs1,'color',thisCols_all{1},'linewidth',1.5);
plot(xs,cs2,'color',thisCols_all{2},'linewidth',1.5);
plot(xs,cs3,'color',thisCols_all{3},'linewidth',1.5);
plot(xs,cs4,'color',thisCols_all{4},'linewidth',1.5);
xlim([bins(1) bins(end)]);
ylim([0 100]);
set(gca,'TickDir','out','FontSize',8,'FontWeight','Bold');

