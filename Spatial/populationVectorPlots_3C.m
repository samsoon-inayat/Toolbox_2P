function populationVectorPlots(tei)
% close all
thisCols_all = {[0 0 0],'b','r',[0 0.7 0.3],[0 0.7 0.3],'r',[0 0.7 0.3],[0 0 0],'b','m','r'};
if ~exist('tei','var')
    ei = evalin('base','ei');
    tei = ei([1 2]);
end
owsi = 0;
fileName = makeName('trialsRastersMetaData.mat',tei{1}.folders.thispFolder);
trialsRastersMD = load(fileName);
trials = trialsRastersMD.trials;
% trials = {[],[1:15],[16:31],[32:47]};
trl1 = trials{1}; trl2 = trials{2}; trl3 = trials{3}; trl4 = trials{4};
% trl1 = 1:10; trl3 = 11:26; trl4 = 27:42;
if ~isempty(trl1)
    nbb_air = getRasters_air(tei,trl1);
    temp = find_mutual_information(tei,nbb_air,owsi); nbb_air.SI = temp.zMI;
    nbb_belt = getRasters_belt(tei,trl1);
    temp = find_mutual_information(tei,nbb_belt,owsi); nbb_belt.SI = temp.zMI;
end
nbc_air = getRasters_air(tei,trl2);
temp = find_mutual_information(tei,nbc_air,owsi); nbc_air.SI = temp.zMI;
nbc_belt = getRasters_belt(tei,trl2);
temp = find_mutual_information(tei,nbc_belt,owsi); nbc_belt.SI = temp.zMI;
bc_air_belt = getRasters_air_belt(tei,trl3);
temp = find_mutual_information(tei,bc_air_belt,owsi); bc_air_belt.SI = temp.zMI;
bb_air_belt = getRasters_air_belt(tei,trl4);
temp = find_mutual_information(tei,bb_air_belt,owsi); bb_air_belt.SI = temp.zMI;

% ff = makeFigureRowsCols(10,[1 1 6 2.5],'RowsCols',[1 4],'spaceRowsCols',[-0.009 0.0009],'rightUpShifts',[0.04 0.01],'widthHeightAdjustment',[0 0]);
rows = 4;cols = 2;
totalPlots = rows * cols;
plotNums = 1:totalPlots;
plotInds = reshape(plotNums,cols,rows)';
varNames = {'nbc_air','nbc_belt','bc_air_belt','bb_air_belt'};
cmdText = sprintf('ptc = findMeanRasters(%s);',varNames{3});
eval(cmdText);
ccs = find(bc_air_belt.SI>0);
% ccs = find(bb_air_belt.SI>5);
[~,peakPos] = max(ptc(ccs,:),[],2);[~,cellNums] = sort(peakPos);
hf = figure(20);clf;
for ii = 1:4
    cmdText = sprintf('ptc = findMeanRasters(%s);',varNames{ii});
    eval(cmdText);
    ptc = ptc(ccs,:);
%     normSel = 1;
%     [~,peakPos] = max(ptc,[],2);[~,cellNums] = sort(peakPos);
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

