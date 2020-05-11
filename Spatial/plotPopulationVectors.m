function plotPopulationVectors(fn,allRs,ccs,sorting)
FS = 9;
if exist('sorting','var')
    ptc = findMeanRasters(allRs{sorting});
    ptc = ptc(ccs,:);
    ptc = normalizeSignal(ptc,2);
    [~,peakPos] = max(ptc,[],2);
    [~,cellNums] = sort(peakPos);
end

rows = length(allRs);

indices = 1:(2*rows);
indices = reshape(indices,2,rows);
indices = indices';
figure(fn);clf;
for ii = 1:rows
    Rs = allRs{ii};
    if exist('sorting','var')
        [P,C] = findPopulationVectorPlot(Rs,ccs,cellNums);
    else
        [P,C] = findPopulationVectorPlot(Rs,ccs);
    end
    subplot(rows,2,indices(ii,1));imagesc(P);colorbar

    colormap parula;
        % caxis([minC maxC]);
    set(gca,'Ydir','Normal','linewidth',1.5,'FontSize',FS,'FontWeight','Bold');
    if isfield(Rs,'dist')
        xlabel('Position (cm)');
%         ylabel({sprintf('%s',Rs.name),'Cells'});
        ylabel('Cell Number');
        cols = size(Rs.rasters(:,:,1),2);
        colsHalf = ceil(cols/2);
        ts = round(Rs.dist);
        set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
        xlabel('Position (cm)');
    else
        xlabel('Time (secs)');ylabel({sprintf('%s',Rs.name),'Cells'});
        cols = size(Rs.rasters(:,:,1),2);
        colsHalf = ceil(cols/2);
        ts = round(Rs.cells(1).times - Rs.cells(1).times(colsHalf));
        set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    end

    subplot(rows,2,indices(ii,2));imagesc(C,[-1 1]);
    colorbar
    colormap parula;
    % caxis([minC maxC]);
    axis equal
    set(gca,'Ydir','Normal','linewidth',1.5,'FontSize',FS,'FontWeight','Bold');
    if isfield(Rs,'dist')
        xlabel('Position (cm)');ylabel('Position (cm)');
        cols = size(Rs.rasters(:,:,1),2);
        colsHalf = ceil(cols/2);
        ts = round(Rs.dist);
        set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
        set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1) ts(colsHalf) ts(cols)],'Ydir','Normal');
    else
        xlabel('Time (secs)');ylabel('Time (secs)');
        cols = size(Rs.rasters(:,:,1),2);
        colsHalf = ceil(cols/2);
        ts = round(Rs.cells(1).times - Rs.cells(1).times(colsHalf));
        set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
        set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1) ts(colsHalf) ts(cols)],'Ydir','Normal');
    end
end


function [ptc,CRc] = findPopulationVectorPlot(Rs,ccs,cellNums)
ptc = findMeanRasters(Rs);
ptc = ptc(ccs,:);
ptc = normalizeSignal(ptc,2);
if ~exist('cellNums','var')
    [~,peakPos] = max(ptc,[],2);
    [~,cellNums] = sort(peakPos);
    ptc = ptc(cellNums,:);
else
    ptc = ptc(cellNums,:);
end
CRc = corrcoef(ptc);


