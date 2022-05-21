function checking_ca_signal_from_within_roi
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei = evalin('base','d15'); Rs = evalin('base','raster_data'); 
sC = evalin('base','selContexts'); rN = evalin('base','rasterNames');
 n = 0;
 %%
 
an = 1; pp = 1;
rec = ei{an};
pla = rec.plane{pp}.tP;
ops = pla.ops;
stat = pla.stat;
xrange = ops.xrange;
yrange = ops.yrange;
mImg = double(ops.meanImgE);
% mImg = double(ops.max_proj);
% mImg = mImg(yrange(1):yrange(2),xrange(1):xrange(2));
ccs = find(pla.iscell(:,1));
roiN = ccs(17);
roi = stat{roiN};
xpix = double(roi.xpix);
ypix = double(roi.ypix);
xlims = [min(xpix) max(xpix)] + [-50 20];
ylims = [min(ypix) max(ypix)] + [-50 20];

% xlims = [1 size(mImg,2)];
% ylims = [1 size(mImg,1)];

figure(100);clf;
imagesc(mImg);
xlim(xlims);
ylim(ylims);
% xlim(xrange);
% ylim(yrange);

maskZ = zeros(size(mImg));
ipix = sub2ind(size(mImg),ypix,xpix);
maskZ(ipix) = 1;
maskZ = expandOrCompressMask(maskZ,0.25);
% maskZ = maskZ';
figure(101);clf;
imagesc(0.75*max(mImg(:))*maskZ+mImg);
xlim(xlims);
ylim(ylims);
