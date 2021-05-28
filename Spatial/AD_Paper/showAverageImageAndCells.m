function showAverageImageAndCells(ei,pl)


mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
% cellsOrNot = 1; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = 1.96; fwids = [1 150]; fcens = [0 150]; rs_th = 0.3; FR = NaN;
conditionsAndRasterTypes = [11;21;31;41];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN,FR);
out = read_data_from_base_workspace_AD(selC)

ei_C = out.eis{1}; ei_A = out.eis{2};
pMs_C = out.pMs{1}; pMs_A = out.pMs{2};
paramMs_C = out.paramMs{1}; paramMs_A = out.paramMs{2};
cpMs_C = out.cpMs{1}; cpMs_A = out.cpMs{2};
selAnimals_C = out.selAnimals{1}; selAnimals_A = out.selAnimals{2};
perc_cells_C = out.perc_cells{1}; perc_cells_A = out.perc_cells{2};

writeToFile = 1;

an_C = 3; an_A = 3;

h = figure(102);clf;
showAverageImage(ei_C{an_C},1);
set(h,'units','inches');
set(h,'Position',[5 4 1 1]);
changePosition(gca,[-0.195 -0.147 0.267 0.27]);
if writeToFile
save_pdf(h,mData.pdf_folder,sprintf('cellsAvgImg_C.pdf'),600);
end

h = figure(103);clf;
showROIs(ei_C{an_C},1);
set(h,'units','inches');
set(h,'Position',[6.5 4 1 1]);
changePosition(gca,[-0.195 -0.147 0.267 0.27]);
if writeToFile
save_pdf(h,mData.pdf_folder,sprintf('cellsAvgImg_C_ROIs.pdf'),600);
end


h = figure(104);clf;
showAverageImage(ei_A{an_A},1);
set(h,'units','inches');
set(h,'Position',[9 4 1 1]);
changePosition(gca,[-0.195 -0.147 0.267 0.27]);
if writeToFile
save_pdf(h,mData.pdf_folder,sprintf('cellsAvgImg_A.pdf'),600);
end

h = figure(105);clf;
showROIs(ei_A{an_A},1);
set(h,'units','inches');
set(h,'Position',[10.5 4 1 1]);
changePosition(gca,[-0.195 -0.147 0.267 0.27]);
if writeToFile
save_pdf(h,mData.pdf_folder,sprintf('cellsAvgImg_A_ROIs.pdf'),600);
end


function showAverageImage(ei,pl)
FS = 8;
micronsPerPixel = ei.thorExp.widthUM/ei.thorExp.pixelX;
xrange = ei.plane{pl}.tP.ops.xrange;
yrange = ei.plane{pl}.tP.ops.yrange;
mimg = ei.plane{pl}.tP.ops.meanImgE;
% mimg = ei.plane{pl}.tP.ops.refImg;
imagesc(mimg,[min(mimg(:)) max(mimg(:))]);
colormap gray;
% axis equal;
axis off;
hold on;
nPixels = round(50/micronsPerPixel);
scaleBarS = xrange(2)- 20 - nPixels; scaleBarE = scaleBarS + nPixels;
scaleBarY = yrange(2) - 30;
plot([scaleBarS,scaleBarE],[scaleBarY,scaleBarY],'r','linewidth',1.5);
xlim(xrange); ylim(yrange);

function showROIs(ei,pl)
FS = 8;
micronsPerPixel = ei.thorExp.widthUM/ei.thorExp.pixelX;
xrange = ei.plane{pl}.tP.ops.xrange;
yrange = ei.plane{pl}.tP.ops.yrange;
mimg = double(ei.plane{pl}.tP.ops.meanImgE);
maskZ = zeros(size(mimg));
colormap gray;
% axis equal;
axis off;
hold on;
allmask = maskZ;
ccs = find(ei.plane{pl}.tP.iscell(:,1));
for cc = 1:length(ccs)
    stat =  ei.plane{pl}.tP.stat(ccs(cc));
    mask = maskZ;
    xpix = stat{1}.xpix; ypix = stat{1}.ypix;
    ipix = sub2ind(size(mask),ypix,xpix);
    mask(ipix) = 1;
%     mask = expandOrCompressMask(mask,stat{1}.footprint);
%     mask = mask';
    allmask(mask==1) = 1;
end
% mimg = ei.plane{pl}.tP.ops.refImg;
imagesc((allmask));
nPixels = round(50/micronsPerPixel);
scaleBarS = xrange(2)- 20 - nPixels; scaleBarE = scaleBarS + nPixels;
scaleBarY = yrange(2) - 30;
plot([scaleBarS,scaleBarE],[scaleBarY,scaleBarY],'r','linewidth',1.5);
xlim(xrange); ylim(yrange);
set(gca,'Ydir','reverse');
