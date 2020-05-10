function cellCorrelation(ei,fn)

thisCols_all = {[0 0 0],'b','r',[0 0.7 0.3],[0 0.7 0.3],'r',[0 0.7 0.3],[0 0 0],'b','m','r'};

aei = evalin('base','ei');
if ~exist('ei','var')
    ei = aei([1]);
    fn = 10002;
end
cei = combineRecordings(ei);
% cei.cellsN = size(ei{1}.signals,1);
% cei.cellsInds = 1:cei.cellsN;

mData.cellsN = cei.cellsN;
mData.cellsInds = cei.cellsInds;
mData.sith = [10 10 10 10 5 5 5 5];

tei{1} = cei;

[a_ddata,a_tdata] = getDataContexts(cei,1:9,'air');
[b_ddata,b_tdata] = getDataContexts(cei,1:9,'belt');
data = {a_ddata{2} b_ddata{2:4}};
[mon_ddata,mon_tdata] = getDataContexts(cei,1:9,'motionOnsets22');
selCells = selectCells({b_ddata{3}},mData,1,[1],[5]);
showCells(100,tei{1},selCells)
% selCells = selectCells({mon_tdata{8}},mData,1,[1],[3]);

frames2 = find(tei{1}.b.frames_f >  tei{1}.contexts(2).markers.air_onsets(1) & tei{1}.b.frames_f <  tei{1}.contexts(2).markers.air_offsets(end));
frames3 = find(tei{1}.b.frames_f >  tei{1}.contexts(3).markers.air_onsets(1) & tei{1}.b.frames_f <  tei{1}.contexts(3).markers.air_offsets(end));
frames4 = find(tei{1}.b.frames_f >  tei{1}.contexts(4).markers.air_onsets(1) & tei{1}.b.frames_f <  tei{1}.contexts(4).markers.air_offsets(end));


for ii = 1:length(tei{1}.deconv.caSigAll)
    cellSigs(:,ii) = tei{1}.deconv.caSigAll{ii};
end

% for ii = 1:length(tei{2}.deconv.spSigAll)
%     cellSigs(:,ii+length(tei{1}.deconv.spSigAll)) = tei{2}.deconv.spSigAll{ii};
% end
dists = getDistanceBetweenCells(tei{1},selCells);
CRc2 = corrcoef(cellSigs(frames2,selCells));
CRc3 = corrcoef(cellSigs(frames3,selCells));
CRc4 = corrcoef(cellSigs(frames4,selCells));
ccth = 0.5;
CRc2(CRc2<ccth) = 0;
CRc3(CRc3<ccth) = 0;
CRc4(CRc4<ccth) = 0;
[rows,cols] = find(triu(CRc3) > ccth);
rowsC = selCells(rows);
colsC = selCells(cols);
for cc = 1:length(cols)
    if rows(cc) == cols(cc)
        continue;
    end
    plotRastersMulti(data,[rowsC(cc) colsC(cc)],0,0,0);
    figure(2000);clf;
    sig2 = cellSigs(frames2,colsC(cc));
    sig1 = cellSigs(frames2,rowsC(cc));
    subplot 311
    plot(normalizeSignal(sig1));hold on;
    plot(normalizeSignal(sig2));
    title(sprintf('%d - %d - [%.2f %.2f %.2f]',rowsC(cc),colsC(cc),CRc2(rows(cc),cols(cc)),CRc3(rows(cc),cols(cc)),CRc4(rows(cc),cols(cc))));
    sig2 = cellSigs(frames3,colsC(cc));
    sig1 = cellSigs(frames3,rowsC(cc));
    subplot 312
    plot(normalizeSignal(sig1));hold on;
    plot(normalizeSignal(sig2));
    cchere = corrcoef(sig1,sig2);
    title(sprintf('%.2f',cchere(1,2)));
    sig2 = cellSigs(frames4,colsC(cc));
    sig1 = cellSigs(frames4,rowsC(cc));
    subplot 313
    plot(normalizeSignal(sig1));hold on;
    plot(normalizeSignal(sig2));
    cchere = corrcoef(sig1,sig2);
    title(sprintf('%.2f',cchere(1,2)));
    showCells(100,tei{1},[rowsC(cc) colsC(cc)])
    n =0 ;
end

% CRC = xcorr(cellSigs,'maxlag',0)
figure(2000);clf
subplot(1,3,1);imagesc(CRc2,[-1 1]);colorbar
subplot(1,3,2);imagesc(CRc3,[-1 1]);colorbar
subplot(1,3,3);imagesc(CRc4,[-1 1]);colorbar
colormap jet
% % subplot(1,2,1);
% % imagesc(cellSigs);colorbar
% % colormap jet;
% % caxis([minC maxC]);
% % set(gca,'Ydir','Normal','linewidth',1.5,'FontSize',12,'FontWeight','Bold');
% % xlabel('Position (cm)');ylabel('Cells');
% % set(gca,'XTick',[1 25 50],'XTickLabel',[0 71 142]);
% %     text(73,100,'Normalized Spike Rate','Rotation',90,'FontWeight','Bold');
% 
% % subplot(1,2,2);
% imagesc(CRc,[0 0.2]);colorbar
% colormap jet;
% caxis([minC maxC]);
% set(gca,'Ydir','Normal','linewidth',1.5,'FontSize',12,'FontWeight','Bold');
% axis equal
% set(gca,'XTick',[1 25 50],'YTick',[1 25 50],'XTickLabel',[0 71 142],'YTickLabel',[0 71 142]);
% xlabel('Position (cm)');ylabel('Position (cm)');
%     text(76,-5,'Pearson Correlation Coefficient','Rotation',90,'FontWeight','Bold');
allInds = 1:size(CRc2,1);
for ii = 1:size(CRc2,1)
    n = 0;
    thisCs1 = CRc2(ii,:);%.*dists(ii,:);
    thisCs2 = CRc3(ii,:);%.*dists(ii,:);
    thisCs3 = CRc4(ii,:);%.*dists(ii,:);
    inds = allInds; inds(ii) = [];
    allCs(ii,:) = [mean(thisCs1(inds)) mean(thisCs2(inds)) mean(thisCs3(inds))];
    dAllCs(ii,:) = abs([sum(thisCs2(inds) - thisCs1(inds)) sum(thisCs3(inds) - thisCs1(inds)) sum(thisCs3(inds) - thisCs2(inds))]);
end
data{1}.vals = dAllCs(:,1)'; data{2}.vals = dAllCs(:,2)'; data{3}.vals = dAllCs(:,3)';
data{1}.name = '1'; data{2}.name = '2'; data{3}.name = '3'; 
data{1}.vals = allCs(:,1)'; data{2}.vals = allCs(:,2)'; data{3}.vals = allCs(:,3)';
plotDistributions(fn,data,20);
return;
for ii = 1:size(CRc2,1)
    n = 0;
    thisCs1 = CRc2(ii,:);
    thisCs2 = CRc3(ii,:);
    thisCs3 = CRc4(ii,:);
    inds = allInds; inds(ii) = [];
    allCs = [thisCs1(inds) thisCs2(inds) thisCs3(inds)];
    if max(allCs) == 0
        continue;
    end
    scale = [min(allCs) max(allCs)]
    [ccimg1,bound] = showCs(10,tei{1},thisCs1,ii,scale,selCells);
    ccimg2 = showCs(11,tei{1},thisCs2,ii,scale,selCells);
    ccimg3 = showCs(12,tei{1},thisCs3,ii,scale,selCells);
    figure(10002);clf;
    subplot 121;
    imagesc(ccimg2-ccimg1);axis equal; hold on;
    plot(bound(:,2),bound(:,1),'w');
    subplot 122;
    imagesc(ccimg3-ccimg2);axis equal; hold on;
    plot(bound(:,2),bound(:,1),'w');
    n = 0;
end

n = 0;


function [ccimg,bound] = showCs(ah,ei,thisCs,ii,scale,selCells)
FS = 12;
micronsPerPixel = ei.widthUM/ei.pixelX;
xrange = ei.ops1{1}.xrange;
yrange = ei.ops1{1}.yrange;
mimg = ei.ops1{1}.mimg;
ccimg = zeros(size(mimg));
tcimg = ccimg;
% axes(ah);
% imagesc(mimg);hold on;
% colormap gray
ccimg = zeros(size(mimg));
ccs = ei.areCells(selCells);
for cc = 1:length(ccs)
    cellX = ei.tP.stat(ccs(cc)).xpix + min(xrange);
    cellY = ei.tP.stat(ccs(cc)).ypix + min(yrange);
    if cc == ii
        mx = min(cellX); Mx = max(cellX); dx = (Mx-mx)/2;
        my = max(cellY); My = max(cellY); dy = (My-my)/2;
        for jj = 1:length(cellX)
            tcimg(cellY(jj),cellX(jj)) = 1;
        end
        B = bwboundaries(tcimg);
        continue;
    end
%     plot(cellX,cellY,'.','color','r');
    for jj = 1:length(cellX)
        ccimg(cellY(jj),cellX(jj)) = thisCs(cc);
    end
end
figure(ah);clf;
imagesc(ccimg,scale);colorbar
hold on;
bound = B{1};
plot(B{1}(:,2),B{1}(:,1),'w');
text(mx+dx,my+dy,num2str(ccs(ii)),'color','w','fontsize',FS);
colormap jet;
axis equal
