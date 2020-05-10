function dists = getDistanceBetweenCells (ei,selCells,overwrite)

if ~exist('overwrite','var')
    overwrite = 0;
end

fileName = makeName(sprintf('cellDists.mat'),ei.folders.thispFolder);
if exist(fileName,'file') && overwrite == 0
    load(fileName);
%     rasters = findOccupancyNormalizedRasters(rasters,ei.deconv.spSigAll);
    if exist('selCells','var')
        dists = dists(selCells,selCells);
    end
    return;
end

FS = 12;
xrange = ei.ops1{1}.xrange;
yrange = ei.ops1{1}.yrange;
% mimg = ei.ops1{1}.mimg;
% figure(2000);clf;
% imagesc(mimg,[min(mimg(:)) max(mimg(:))-1000]);
% colormap gray;
% axis equal;
% axis off;
% hold on;
if ~exist('selCells','var')
    selCells = 1:length(ei.areCells);
end
ccs = ei.areCells(selCells);
combs = nchoosek(ccs,2);
dists = zeros(length(selCells),length(selCells));

display('Finding distance matrix');
for ii = 1:size(combs,1)
    ccs1 = combs(ii,1);
    ccs2 = combs(ii,2);
    if iscell(ei.tP.stat)
        stat1 = ei.tP.stat{ccs1};
        stat2 = ei.tP.stat{ccs2};
    else
        stat1 = ei.tP.stat(ccs1);
        stat2 = ei.tP.stat(ccs2);
    end
    cellX = double(stat1.xpix) + double(min(xrange));
    cellY = double(stat1.ypix) + double(min(yrange));
    temp = zeros(max(cellY-min(cellY))+1,max(cellX-min(cellX))+1);
    for jj = 1:length(cellX)
        temp(cellY(jj)-min(cellY)+1,cellX(jj)-min(cellX)+1) = 1;
    end
    cent = regionprops(temp,'Centroid');
    C1 = cent.Centroid + [min(cellX)-1 min(cellY)-1];
    cellX = double(stat2.xpix) + double(min(xrange));
    cellY = double(stat2.ypix) + double(min(yrange));
    temp = zeros(max(cellY-min(cellY))+1,max(cellX-min(cellX))+1);
    for jj = 1:length(cellX)
        temp(cellY(jj)-min(cellY)+1,cellX(jj)-min(cellX)+1) = 1;
    end
    cent = regionprops(temp,'Centroid');
    C2 = cent.Centroid + [min(cellX)-1 min(cellY)-1];
    iC1 = find(selCells == find(ei.areCells == ccs1));
    iC2 = find(selCells == find(ei.areCells == ccs2));
    dists(iC1,iC2) = sqrt(sum((C2-C1).^2));
% %     figure(2001);clf;imagesc(temp);axis equal;hold on;
%     plot(cellX,cellY,'.','color','r');
%     plot(C1(1),C1(2),'*','color','m');
%     plot(C2(1),C2(2),'*','color','c');
%     n = 0;
end
distsT = dists';
dists = dists + distsT;
save(fileName,'dists');
display('Done!!!');
