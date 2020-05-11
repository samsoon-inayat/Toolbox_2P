function plotPopulationVectors(fn,allRs,allRs1,ccs,sorting)
if iscell(ccs)
    ccs1 = ccs{1};
    ccs2 = ccs{2};
else
    ccs1 = ccs;
    ccs2 = ccs;
end
FS = 9;
if exist('sorting','var')
    ptc = findMeanRasters(allRs{sorting(1)});
    ptc = ptc(ccs1,:);
    ptc = normalizeSignal(ptc,2);
    [~,peakPos] = max(ptc,[],2);
    [~,cellNums1] = sort(peakPos);
    if length(sorting)>1
        ptc = findMeanRasters(allRs1{sorting(2)});
        ptc = ptc(ccs2,:);
        ptc = normalizeSignal(ptc,2);
        [~,peakPos] = max(ptc,[],2);
        [~,cellNums2] = sort(peakPos);
    end
end

rows = length(allRs);
scols = 4;
indices = 1:(scols*rows);
indices = reshape(indices,4,rows);
indices = indices';
figure(fn);clf;
for ii = 1:rows
    Rs = allRs{ii};
    if exist('sorting','var')
        [P,C] = findPopulationVectorPlot(Rs,ccs1,cellNums1);
    else
        [P,C,cellNumsC{ii}] = findPopulationVectorPlot(Rs,ccs1);
    end
    subplot(rows,scols,indices(ii,1));imagesc(P);colorbar

    colormap jet;
        % caxis([minC maxC]);
    set(gca,'Ydir','Normal','linewidth',1.5,'FontSize',FS,'FontWeight','Bold');
    if isfield(Rs,'dist')
        xlabel('Position (cm)');ylabel({sprintf('%s',Rs.name),'Cells'});
        cols = size(Rs.rasters(:,:,1),2);
        colsHalf = ceil(cols/2);
        ts = round(Rs.dist);
        set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
        xlabel('Position (cm)');
    else
        xlabel('Time (secs)');ylabel({sprintf('%s',Rs.name),'Cells'});
        cols = size(Rs.rasters(:,:,1),2);
        colsHalf = ceil(cols/2);
        ts = round(Rs.cells(1).times);% - Rs.cells(1).times(colsHalf));
        set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    end

    subplot(rows,scols,indices(ii,2));imagesc(C,[-1 1]);
    colorbar
    colormap jet;
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
        ts = round(Rs.cells(1).times);% - Rs.cells(1).times(colsHalf));
        set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
        set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1) ts(colsHalf) ts(cols)],'Ydir','Normal');
    end
end

for ii = 1:rows
    Rs = allRs1{ii};
    if exist('sorting','var')
        [P,C] = findPopulationVectorPlot(Rs,ccs2,cellNums2);
    else
        [P,C] = findPopulationVectorPlot(Rs,ccs2,cellNumsC{ii});
    end
    subplot(rows,scols,indices(ii,3));imagesc(P);colorbar

    colormap jet;
        % caxis([minC maxC]);
    set(gca,'Ydir','Normal','linewidth',1.5,'FontSize',FS,'FontWeight','Bold');
    if isfield(Rs,'dist')
        xlabel('Position (cm)');ylabel({sprintf('%s',Rs.name),'Cells'});
        cols = size(Rs.rasters(:,:,1),2);
        colsHalf = ceil(cols/2);
        ts = round(Rs.dist);
        set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
        xlabel('Position (cm)');
    else
        xlabel('Time (secs)');ylabel({sprintf('%s',Rs.name),'Cells'});
        cols = size(Rs.rasters(:,:,1),2);
        colsHalf = ceil(cols/2);
        ts = round(Rs.cells(1).times);% - Rs.cells(1).times(colsHalf));
        set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
    end

    subplot(rows,scols,indices(ii,4));imagesc(C,[-1 1]);
    colorbar
    colormap jet;
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
        ts = round(Rs.cells(1).times);% - Rs.cells(1).times(colsHalf));
        set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
        set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1) ts(colsHalf) ts(cols)],'Ydir','Normal');
    end
end




function [ptc,CRc,cellNums] = findPopulationVectorPlot(Rs,ccs,cellNums)
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


