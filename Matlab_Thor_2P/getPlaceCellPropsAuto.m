function pfs = getPlaceCellPropsAuto(allRs,ids,selCells)
plotFlag = 1;
Rs = allRs{ids(1)};
n = 0;

raster = Rs.rasters(:,:,1);
mSig = nanmean(raster);
xs = 1:length(mSig);
uxs = linspace(1,length(mSig),length(mSig)*10);

for ii = 1:length(selCells)
    thisCell = selCells(ii);
    raster = Rs.rasters(:,:,thisCell);
    mSig = nanmean(raster);
    umSig = interp1(xs,mSig,uxs);
    umSig = applyGaussFilt(umSig,10);
    
    [pks,locs,wpk,pkp] = findpeaks(umSig,'MinPeakDistance',50,'MinPeakProminence',50);
    thresh = mean(umSig) + 0.5*std(umSig);
    inds = find(pks > thresh);
    pks = pks(inds);locs = locs(inds);wpk = wpk(inds);pkp = pkp(inds);
    nPeaks(ii) = length(inds);
    plotTheCurves(xs,uxs,umSig,raster,Rs.dist(end)/size(raster,2),locs);
    n=0;
end


function plotTheCurves(xs,uxs,umSig,ccSignalA,cm_per_bin,locs)
figure(10002);clf
subplot 221;imagesc(ccSignalA);set(gca,'XTick',(1:5:50),'XTickLabel',round(((1:5:50))*cm_per_bin));
subplot 222;imagesc(corrcoef(ccSignalA));
colorbar;
set(gca,'Ydir','Normal');
subplot(2,2,3);
plot(uxs*cm_per_bin,umSig,'k');hold on;
plot([uxs(1) uxs(end)]*cm_per_bin,[mean(umSig) mean(umSig)],'m');
plot(uxs(locs)*cm_per_bin,umSig(locs),'*r');
set(gca,'XTick',round((uxs(1:50:500))*cm_per_bin));
xlim([uxs(1) uxs(end)]*cm_per_bin);
subplot(2,2,4);
plot(xs*cm_per_bin,nanmean(ccSignalA).*nanstd(ccSignalA));
% 
% E = imbinarize(ccSignalA,nanmean(ccSignalA(:))+1*nanstd(ccSignalA(:)));
% imagesc(E);colorbar;
% title(sprintf('%.3f - %.3f',skm,k
