function figure_tracePlots_2 (ei)
%%
adata = evalin('base','data');
aei = evalin('base','sei');
mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;
% allCells = mData.allCells;
selAnimals = 1:3;
n = 0;

% allCells = 1:sum(mData.cellsN);
data_1 = adata(1:3);
data = {data_1{1:2}};
% data = populate_mean_raster_fitting(data);

n = 0;


%%
ei = aei{1};
pp = 1;
b = ei.b; b.frames_f = ei.plane{pp}.b.frames_f;
signals = ei.plane{pp}.tP.signals;
spSigAll = ei.plane{pp}.tP.deconv.spSigAll;
caSigAll = ei.plane{pp}.tP.deconv.caSigAll;
coiSI = find(ei.plane{pp}.contexts(1).rasters.airD.SI>5);
rasters = ei.plane{pp}.contexts(1).rasters.airD.rasters;
allSIs = ei.plane{pp}.contexts(1).rasters.airD.SI(coiSI);
% labelCells(ei,find(coiSI),allSIs);
ccsi = coiSI(randi([1 length(coiSI)],1,6));
% ccsi = [488 136 319 39 295 77    78]; %137, 39, 77, 5, 136, 429, 319
% ccsi = [39,48,390,429,248,136];
ccsi = [21 24 7 18 40 105];
% ccsi = sort(ccsi);
n = 0;
%%
A = data{1}{1}{1};

numberOfRows = 2;
numberOfCols = 3;
graphsOnOneFigure = numberOfRows * numberOfCols;
numberOfData = length(ccsi);
numberOfGroups = ceil(numberOfData/graphsOnOneFigure);
numberOfFullGroups = floor(numberOfData/graphsOnOneFigure);
indices = NaN(1,(numberOfGroups*graphsOnOneFigure));
indices(1:numberOfData) = 1:numberOfData;
groupIndices = reshape(indices,numberOfRows,numberOfCols,numberOfGroups);
% for gg = 1:numberOfGroups
%     groupIndices(:,:,gg) = groupIndices(:,:,gg)';
% end

ff = makeFigureRowsCols(105,[0.5 0.5 4 1],'RowsCols',[numberOfRows numberOfCols],...
    'spaceRowsCols',[0.07 0.06],'rightUpShifts',[0.1 0.1],'widthHeightAdjustment',...
    [-90 -180]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[10 4 2.75 2]);
while gg<=numberOfGroups
    for rr = 1:numberOfRows
        for cc = 1:numberOfCols
            cn = ccsi(groupIndices(rr,cc));
            axes(ff.h_axes(rr,cc));
            thisRaster = rasters(:,:,cn);
            
            mSig = nanmean(thisRaster);
            ft = fittype(A.formula);
            xs = 1:size(A.rasters,2);
            coeff = A.coeff(:,cn);
            fitplot = ft(coeff(1),coeff(2),coeff(3),xs);
            
            imagesc(thisRaster);hold on;
            plot(size(thisRaster,1)*fitplot/max(fitplot),'linewidth',0.5,'color','r');
            plot(size(thisRaster,1)*nanmean(thisRaster)/max(nanmean(thisRaster)),'linewidth',0.25,'color','w');
%             hc = colorbar;
            box off;
%             colorbarPos = get(hc,'Position');
%             colorbarPos = colorbarPos + [0.11 0 -0.01 0];
%             set(hc,'Position',colorbarPos);
            if rr < 2
                set(gca,'XTick',[]);
            else
                set(gca,'XTick',[1 24 48],'XTickLabel',[0 70 140]);
                if cc == 2
                    hx = xlabel('Position (cm)');
                    changePosition(hx,[0 -18 0]);                    
                end
            end
            if cc >1
                set(gca,'YTick',[]);
            end
            if cc == 1
                h = ylabel('Trials');
                changePosition(h,[8 0 0]);
                text
            end
%                 text(65,15,'Occupancy Normalized Spike Rate (Hz/sec)','Rotation',90,'FontWeight','Normal','FontSize',5);
            text(1,12.2,sprintf('Cell %d - Max FR %d',cn,round(max(thisRaster(:)))),'FontSize',5,'color','k');
            text(53,4,sprintf('MI = %.2f',allSIs(coiSI==cn)),'FontSize',5,'color','k','rotation',90);
            set(gca,'FontSize',7,'FontWeight','Normal','linewidth',0.75,'Ydir','normal','TickDir','out');
            if rr == 1 & cc == 1
                text(-11,15,'Rasters - Occupancy normalized firing rate (FR)','FontSize',6.5);
            end
            if cc == 3 && rr == 1
                pos = get(gca,'Position');
                h = axes('Position',pos);
                changePosition(h,[0.03 0.33 0.02 -0.1]);
                hc = colorbar('north'); axis off
                set(hc,'YTick',[],'linewidth',0.01);
                changePosition(hc,[0 0 0 -0.01]);
                ylims = ylim;
                xlims = xlim;
                text(xlims(1)+0.15,ylims(1)+0.75,'0','FontSize',4.5);
                text(xlims(2)-0.4,ylims(1)+0.75,'Max FR','FontSize',4.5);
            end
        end
    end
    gg = gg + 1;
end

figure(105);
save_pdf(gcf,mData.pdf_folder,'Sample Rasters',600);
return;

%%
ff = makeFigureRowsCols(100,[1 1 5.4 1.5],'RowsCols',[(length(ccsi))+1 1],'spaceRowsCols',[-0.009 0.0009],'rightUpShifts',[0.04 0.02],'widthHeightAdjustment',[-70 7]);
set(gcf,'Position',[10 4 4.5 1.6]);
set(gcf,'color','w');
ccs = find(ei.plane{pp}.tP.areCells);%(ccsi);
lengthSigs = size(signals,2);
traceTime = ei.b.ts(b.frames_f(1:lengthSigs));
for cc = 1:length(ccsi)
    tsp = spSigAll{ccsi(cc)}';
    alltsp(cc,:) = tsp;
    caSig = signals(ccs(ccsi(cc)),:)';
    allcaSig(cc,:) = caSig;
    sCaSig = caSigAll{ccsi(cc)}';
    allsCaSig(cc,:) = sCaSig;
    minSCaSig(cc) = min(sCaSig); maxSCaSig(cc) = max(sCaSig);
    mintsp(cc) = min(tsp); maxtsp(cc) = max(tsp);
end
minCaSig = min(minSCaSig);
maxCaSig = max(maxSCaSig);

% startTimeFrame = b.frames_f(find(traceTime>300,1,'first'));
% endTimeFrame = b.frames_f(find(traceTime<550,1','last'));

xStart = traceTime(find(traceTime>4,1,'first'));
xEnd = traceTime(find(traceTime<230,1','last'));
xlims = [xStart xEnd];
timeLength = 0.0833*60;
air_puff_raw = b.air_puff_raw;
air_puff_raw(b.air_puff_raw < 0.987) = NaN;
air_puff_raw(b.air_puff_raw >= 0.987) = 1;
onsets = ei.plane{pp}.contexts(1).markers.air_onsets;
offsets = ei.plane{pp}.contexts(1).markers.air_offsets;
photo_sensor_raw = NaN(size(b.photo_sensor_raw));
photo_sensor_raw(b.photo_sensor_f-1) = 0;%max(b.fSpeed)-max(b.fSpeed)/8;
photo_sensor_raw(b.photo_sensor_f) = max(b.fSpeed)/3;

cc = 1;
while cc<=length(ccsi)+1
%     tsp = spSigAll{ccsi(cc)}';
%     caSig = signals(ccs(ccsi(cc)),:)';
%     sCaSig = caSigAll{ccsi(cc)}';
    axes(ff.h_axes(cc,1));
%     get(gca,'Position')
    changePosition(gca,[-0.035 0 0.06 0]);
    if cc<length(ccsi)+1
        tsp = alltsp(cc,:);
%         tsp(tsp==0) = NaN;
        caSig = allcaSig(cc,:);
        sCaSig = allsCaSig(cc,:);
%         plot(traceTime,caSig);hold on;
        plot(traceTime,sCaSig,'b');hold on;
        plot(traceTime,tsp,'r','linewidth',0.25);
        ylims = ylim;
%         plot(b.ts,air_puff_raw*(ylims(2)/3),'-','color','b','linewidth',0.25);hold on;
        ylim([minCaSig maxCaSig]);
        text(xStart+8.5,maxCaSig/2,sprintf('%.2d',ccsi(cc)),'FontSize',4.5);
    else
%         plot(b.ts,0.25*b.photo_sensor_raw+0.55,'color','r');
        plot(b.ts,b.fSpeed,'color','m');hold on;
%         plot(b.ts,air_puff_raw*(max(b.fSpeed)/2),'color','b','linewidth',0.25);
        plot(b.ts,photo_sensor_raw,'color',[0 0.7 0.3]);
%         plot(dpsr*max(b.fSpeed),'color','r');
        n = 0;
%         plot(b.ts,0.5*b.fSpeed/max(b.fSpeed),'color','m');
        ylim([-0.1 max(b.fSpeed)+1.1]);
        
        sdfxs = xEnd - 156; 
        sdfxe = sdfxs + timeLength; sdfys = max(b.fSpeed)/2; sdfye = sdfys + max(b.fSpeed)/4;
        plot([sdfxs sdfxs],[sdfys-10 sdfye-10],'k','linewidth',1);
        text(xStart+8.5,sdfys+5,sprintf('Speed'),'FontSize',4.5);
        text(sdfxs-(0.06*60),sdfye+3,sprintf('%d cm/sec',round(max(b.fSpeed)/4)),'fontsize',4.5);
%         text(sdfxs+(0.015*60),sdfye-15,sprintf('%d cm/sec',round(max(b.fSpeed)/4)),'fontsize',4.5);
    end
%     ylim([min(caSig) max(caSig)]);
    xlim(xlims);
    
    
    axis off;
%     ylabel(num2str(cc));
    text(-25,max(caSig)/4,num2str(cc),'rotation',90);
    if cc == 1
        NcaPerNtsp = max(tsp)/maxCaSig;
        tspbar = NcaPerNtsp * 100;
        sdfxs = xEnd-156; 
        sdfxe = sdfxs + timeLength; sdfys = maxCaSig/1.8; sdfye = sdfys + 100;
        plot([sdfxs sdfxe],[sdfys sdfys],'k','linewidth',0.5);
        plot([sdfxe sdfxe],[sdfys sdfye],'k','linewidth',0.5);
        text(sdfxs-(0.1*60),sdfys + 60,sprintf('%dsecs',ceil(timeLength)),'fontsize',4);
        text(sdfxe+(0.025*60),sdfys + 60,'100% DF/F (25 A.U. Firing Rate)','fontsize',4);
    end
    n = 0;
%     onsets = b.air_puff_r(b.sTrials);
%     offsets = b.air_puff_f(b.sTrials);
% %     tSignals = plotTrials(b,caSig,tsp,onsets,offsets);
% %     tSignals = plotInterTrials(ei,caSig);
%     A = getDistRaster_1(b,caSig,tsp,onsets,offsets,0);
%     lastOffset = offsets(end);
%     
%     
%     ccSignalA = A.distSigRaster./A.distDurRaster;
%     xValsA = A.dists;
%     
%     onsets = b.air_puff_f(b.sTrials(1:(end-1)));
%     offsets = b.air_puff_r(b.sTrials(2:end));
%     AI = getTimeRaster_1(b,caSig,tsp,onsets,offsets,0);
%     figure(101);clf;
%     imagesc(A.distSigRaster./A.distDurRaster);
%     figure(102);clf;
%     imagesc(AI.sigRaster);
    cc = cc + 1;
%     cc = keyboardInput(cc,[1 length(ccsi)],[1 5],'');
%     if cc < 0
%         break;
%     end
end
axes(ff.h_axes(1,1));ylims = ylim;
[TLx TLy] = ds2nfu(b.ts(onsets(1)),ylims(2)-140);
axes(ff.h_axes(length(ccsi)+1,1));ylims = ylim;
[BLx BLy] = ds2nfu(b.ts(onsets(1)),ylims(1));
aH = (TLy - BLy);
for ii = 1:length(onsets)
    [BRx BRy] = ds2nfu(b.ts(offsets(ii)),ylims(1));
    [BLx BLy] = ds2nfu(b.ts(onsets(ii)),ylims(1));
    aW = (BRx-BLx);
    annotation('rectangle',[BLx BLy aW aH],'facealpha',0.2,'linestyle','none','facecolor','b');
end
figure(100);
save2pdf('traces.pdf',gcf,600);
n = 0;
%%
ei.widthUM = ei.thorExp.widthUM;
ei.pixelX = ei.thorExp.pixelX;
ei.ops1{1} = ei.plane{pp}.tP.ops;
ei.areCells = find(ei.plane{pp}.tP.areCells);
ei.tP = ei.plane{pp}.tP;
h = figure(102);plot(0,0);ha = gca;
showCells(gca,ei,ccsi);
h = figure(102);
set(gcf,'units','inches');
set(gcf,'Position',[10 4 1.7 1.7]);
changePosition(gca,[-0.165 -0.167 0.267 0.27]);
save2pdf('cellsAvgImg.pdf',gcf,600);
