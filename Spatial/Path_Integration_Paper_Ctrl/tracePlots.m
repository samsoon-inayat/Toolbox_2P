function tracePlots (ei)

dCols = distinguishable_colors(15,'k');

if ~exist('ei','var')
    ei = evalin('base','ei{1}');
end
%%
b = ei.b;
signals = ei.tP.signals;
spSigAll = ei.deconv.spSigAll;
caSigAll = ei.deconv.caSigAll;
coiSI = find(ei.contexts(1).rasters.airD.SI>5);
rasters = ei.contexts(1).rasters.airD.rasters;
allSIs = ei.contexts(1).rasters.airD.SI(coiSI);
% labelCells(ei,find(coiSI),allSIs);
ccsi = coiSI(randi([1 length(coiSI)],1,6))
% ccsi = [488 136 319 39 295 77    78]; %137, 39, 77, 5, 136, 429, 319
% ccsi = [39,48,390,429,248,136];
ccsi = [21 24 7 18 40 105];
% ccsi = sort(ccsi);
ccs = ei.areCells;%(ccsi);
n = 0;
%%

numberOfRows = 3;
numberOfCols = 2;
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

ff = makeFigureRowsCols(105,[0.5 0.5 15 7],'RowsCols',[numberOfRows numberOfCols],...
    'spaceRowsCols',[0.05 0.15],'rightUpShifts',[0.07 0.12],'widthHeightAdjustment',...
    [-170 -100]);
gg = 1;
set(gcf,'color','w');
while gg<=numberOfGroups
    for rr = 1:numberOfRows
        for cc = 1:numberOfCols
            cn = ccsi(groupIndices(rr,cc));
            axes(ff.h_axes(rr,cc));
            thisRaster = rasters(:,:,cn);
            imagesc(thisRaster); hc = colorbar;
            colorbarPos = get(hc,'Position');
            colorbarPos = colorbarPos + [0.11 0 -0.01 0];
            set(hc,'Position',colorbarPos);
            if rr < 3
                set(gca,'XTick',[]);
            else
                set(gca,'XTick',[1 25 50],'XTickLabel',[0 71 142]);
                xlabel('Position (cm)');
            end
            if cc == 2
                set(gca,'YTick',[]);
            end
            if rr == 2
                if cc == 1
                    ylabel('Trials');
                end
                text(65,15,'Occupancy Normalized Spike Rate (Hz/sec)','Rotation',90,'FontWeight','Normal','FontSize',8);
            end
            text(20,-0.5,sprintf('Cell %d',cn),'FontSize',10);
            set(gca,'FontSize',12,'FontWeight','Bold','linewidth',1.5);
        end
    end
    gg = gg + 1;
end
figure(105);
save2pdf('rasters.pdf',gcf,600);
% return;

%%
ff = makeFigureRowsCols(100,[1 1 5.4 1.5],'RowsCols',[(length(ccsi))+1 1],'spaceRowsCols',[-0.009 0.0009],'rightUpShifts',[0.04 0.02],'widthHeightAdjustment',[-70 7]);
set(gcf,'color','w');
lengthSigs = size(signals,2);
traceTime = ei.b.ts(ei.b.frames_f(1:lengthSigs));
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
timeLength = 0.5*60;
air_puff_raw = b.air_puff_raw;
air_puff_raw(b.air_puff_raw < 0.987) = NaN;
air_puff_raw(b.air_puff_raw >= 0.987) = 1;
photo_sensor_raw = NaN(size(b.photo_sensor_raw));
photo_sensor_raw(b.photo_sensor_f-1) = max(b.fSpeed)-max(b.fSpeed)/8;
photo_sensor_raw(b.photo_sensor_f) = max(b.fSpeed)/2;

cc = 1;
while cc<=length(ccsi)+1
%     tsp = spSigAll{ccsi(cc)}';
%     caSig = signals(ccs(ccsi(cc)),:)';
%     sCaSig = caSigAll{ccsi(cc)}';
    axes(ff.h_axes(cc,1));
    
    if cc<length(ccsi)+1
        tsp = alltsp(cc,:);
        caSig = allcaSig(cc,:);
        sCaSig = allsCaSig(cc,:);
%         plot(traceTime,caSig);hold on;
        plot(traceTime,sCaSig);hold on;
        plot(traceTime,tsp);
        ylims = ylim;
        plot(b.ts,air_puff_raw*(ylims(2)/3),'-','color','b');hold on;
        ylim([minCaSig maxCaSig]);
        text(xStart,maxCaSig/1.5,sprintf('Cell %d',ccsi(cc)),'FontSize',7);
    else
%         plot(b.ts,0.25*b.photo_sensor_raw+0.55,'color','r');
        plot(b.ts,b.fSpeed,'color','m');hold on;
        plot(b.ts,air_puff_raw*(max(b.fSpeed)/2),'color','b');
        plot(b.ts,photo_sensor_raw,'color','r');
%         plot(dpsr*max(b.fSpeed),'color','r');
        n = 0;
%         plot(b.ts,0.5*b.fSpeed/max(b.fSpeed),'color','m');
        ylim([-0.1 max(b.fSpeed)+1.1]);
        
        sdfxs = xStart; 
        sdfxe = sdfxs + timeLength; sdfys = max(b.fSpeed)/2; sdfye = sdfys + max(b.fSpeed)/4;
        plot([sdfxs sdfxs],[sdfys sdfye],'k','linewidth',1);
        text(sdfxs+(0.05*60),sdfye,sprintf('%d cm/sec',round(max(b.fSpeed)/4)),'fontsize',7);
    end
%     ylim([min(caSig) max(caSig)]);
    xlim(xlims);
    
    
    axis off;
%     ylabel(num2str(cc));
    text(-25,max(caSig)/4,num2str(cc),'rotation',90);
    if cc == 1
        NcaPerNtsp = max(tsp)/maxCaSig;
        tspbar = NcaPerNtsp * 100;
        sdfxs = xEnd-200; 
        sdfxe = sdfxs + timeLength; sdfys = maxCaSig/1.8; sdfye = sdfys + 100;
        plot([sdfxs sdfxe],[sdfys sdfys],'k','linewidth',0.5);
        plot([sdfxe sdfxe],[sdfys sdfye],'k','linewidth',0.5);
        text(sdfxs-(0.1*60),sdfys + 60,sprintf('%dsecs',timeLength),'fontsize',7);
        text(sdfxe+(0.1*60),sdfys + 60,'100% DF/F (25 A.U. Firing Rate)','fontsize',7);
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
figure(100);
save2pdf('traces.pdf',gcf,600);

%%
showCells(102,ei,ccsi);
h = figure(102);
set(gcf,'units','inches');
set(gcf,'Position',[10 4 1.5 1.5]);
changePostion(gca,[-0.182 -0.175 0.28 0.28]);
save2pdf('cellsAvgImg.pdf',gcf,600);
