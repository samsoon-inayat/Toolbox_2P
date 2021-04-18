function tracePlots (ei)
%%
ei = evalin('base','ei10_C');
% data = evalin('base','data');
% dataB = evalin('base','datab');
% dataAOn = evalin('base','dataAOn010');
% dataAOff= evalin('base','dataAOff010');
mData = evalin('base','mData');
selAnimals = 1; pl = 1;
tei = ei{selAnimals};
planeNumbers = pl;
maxDistTime = [Inf Inf];
contextNumber = [1 1];
stimMarkers = {'air','air'};
rasterTypes = {'dist','time'};
n = 0;
%%
varName = '';
selCells = 'All';
for ss = 1:length(stimMarkers)
    distD = [];
    for jj = 1:length(selAnimals)
        [tempVals cns ACs] = getParamValues(varName,ei(selAnimals(jj)),planeNumbers,contextNumber(ss),stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
%         [clus] = getParamValues('cluster4',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),stimMarkers{ss},'dist',selCells,maxDistTime);
%         [placeCells5] = getParamValues('placeCells5',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),stimMarkers{ss},'dist',selCells,maxDistTime);
    end
    dataC{ss} = tempVals;
end
n = 0;
A = dataC{1}; B = dataC{2}; %C = dataC{3}; %D = dataC{4};
% rasters = A;
% PCs = ACs & placeCells5;
PCs = ACs;% & clus;
coiSI = find(PCs)
% % % % % sum(PCs)/length(PCs);
% % % % % [~,inds] = sort([A.centers(coiSI)']);
% % % % % coiSI = coiSI(inds); centers = A.centers(coiSI);
% % % % % [coiSI(inds) A.centers(coiSI(inds))' A.pws(coiSI(inds))'];

% coiSI =  coiSI(centers > 1 & centers < 40)';
% coiSI =  coiSI(centers > 40)';

% plot_trial_correlation(dataC,coiSI);return;
plotRasters(dataC,coiSI);return
ccsi = [200 46 152 37 22]; % animal 1 plane 1

%%
runthis = 1;
if runthis
    
cellList = ccsi;
numberOfRows = 2;
numberOfCols = 5;
graphsOnOneFigure = numberOfRows * numberOfCols;
numberOfData = length(ccsi);
numberOfGroups = ceil(numberOfData/graphsOnOneFigure);
numberOfFullGroups = floor(numberOfData/graphsOnOneFigure);
indices = NaN(1,(numberOfGroups*graphsOnOneFigure));
indices(1:numberOfData) = 1:numberOfData;
groupIndices = reshape(indices,numberOfRows,numberOfCols,numberOfGroups);
ff = makeFigureRowsCols(105,[0.5 0.5 4 1],'RowsCols',[numberOfRows numberOfCols],...
    'spaceRowsCols',[0.15 0.03],'rightUpShifts',[0.062 0.12],'widthHeightAdjustment',...
    [-40 -155]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[10 4 4.5 2.5]);
for rr = 1:2
    for cc = 1:5
        cn = cellList(cc);
        axes(ff.h_axes(rr,cc));
        A = dataC{rr};
        thisRaster = A.rasters(:,:,cn);
    end
end
for rr = 1:2
    for cc = 1:5
        cn = cellList(cc);
        axes(ff.h_axes(rr,cc));
        A = dataC{rr};
        thisRaster = A.rasters(:,:,cn);
        mSig = nanmean(thisRaster);
%         ft = fittype(A.formula);
        xs = 1:size(A.rasters,2);
%         coeff = A.coeff(:,cn);
        fitplot = gauss_fit(xs,A.gauss_fit_on_mean.coefficients_Rs_mean(cn,1:3),A.gauss_fit_on_mean.gauss1Formula);
        imagesc(thisRaster,[min(thisRaster(:)) 0.35*max(thisRaster(:))]);hold on;
        plot(size(thisRaster,1)*fitplot/max(fitplot),'linewidth',1.5,'color','r');
        box off;
        text(size(thisRaster,2)+5,0,sprintf('Max FR %d - MI = %.2f - Rs = %.2f',round(max(thisRaster(:))),A.SI(cn),...
            A.gauss_fit_on_mean.coefficients_Rs_mean(cn,4)),'FontSize',4,'color','k','rotation',90);
        if rr == 1
            text(size(thisRaster,1)/2,size(thisRaster,1)+1.25,sprintf('Cell %d',cn),'FontSize',6,'color','k');
        end
        set(gca,'FontSize',6,'FontWeight','Normal','linewidth',0.75,'Ydir','normal','TickDir','out');
        colormap jet;
        if rr == 1
            xticks = [1:20:size(thisRaster,2)];
            set(gca,'XTick',xticks,'XTickLabels',A.xs(xticks));
%                 if cc == 2
                hx = xlabel('Position (cm)');
%                     changePosition(hx,[-50 0 0]);                    
%                 end
        else
            xticks = [1:50:size(thisRaster,2)];
            set(gca,'XTick',xticks,'XTickLabels',A.xs(xticks));
            
                hx = xlabel('Time (sec)');
%                     changePosition(hx,[-20 0 0]);                    
        end

        if cc >1
            set(gca,'YTick',[]);
        end
        if cc == 1
            set(gca,'YTick',[1 5 10]);
            h = ylabel('Trials');
            changePosition(h,[3 0 0]);
            text
        end
        
        if cc == 4 && rr == 2
%             hca = gca;
%             hc = putColorBar(hca,[-5 10 0 0],[0 max(thisRaster(:))],4);
%             set(hc,'Orientation','horizontal');
            
%             pos = get(gca,'Position');
%             h = axes('Position',pos);
% %             changePosition(h,[0.01 0.22 0.02 -0.1]);
%             hc = colorbar('north'); axis off
%             set(hc,'YTick',[],'linewidth',0.01,'box','off');
%             changePosition(hc,[0 0 0 -0.01]);
%             ylims = ylim;
%             xlims = xlim;
%             text(xlims(1)+0.15,ylims(1),'0','FontSize',4.5);
%             text(xlims(2)-0.4,ylims(1),'Max FR','FontSize',4.5);
        end
        
    end
end
% while gg<=numberOfGroups
%     for rr = 1:numberOfRows
%         for cc = 1:numberOfCols
%             cn = ccsi(groupIndices(rr,cc));
%             axes(ff.h_axes(rr,cc));
%             thisRaster = rasters(:,:,cn);
%             
%             mSig = nanmean(thisRaster);
%             ft = fittype(A.formula);
%             xs = 1:size(A.rasters,2);
%             coeff = A.coeff(:,cn);
%             fitplot = ft(coeff(1),coeff(2),coeff(3),xs);
%             
%             imagesc(thisRaster);hold on;
%             plot(size(thisRaster,1)*fitplot/max(fitplot),'linewidth',0.5,'color','r');
% %             hc = colorbar;
%             box off;
% %             colorbarPos = get(hc,'Position');
% %             colorbarPos = colorbarPos + [0.11 0 -0.01 0];
% %             set(hc,'Position',colorbarPos);
%             if rr < 2
%                 set(gca,'XTick',[]);
%             else
%                 set(gca,'XTick',[1 24 48],'XTickLabel',[0 70 140]);
%                 if cc == 2
%                     hx = xlabel('Position (cm)');
% %                     changePosition(hx,[0 50 0]);                    
%                 end
%             end
%             if cc >1
%                 set(gca,'YTick',[]);
%             end
%             if cc == 1
%                 h = ylabel('Trials');
%                 changePosition(h,[3 0 0]);
%                 text
%             end
% %                 text(65,15,'Occupancy Normalized Spike Rate (Hz/sec)','Rotation',90,'FontWeight','Normal','FontSize',5);
%             text(1,11.2,sprintf('Cell %d - Max FR %d',cn,round(max(thisRaster(:)))),'FontSize',5,'color','k');
%             text(53,4,sprintf('MI = %.2f',A.SI(cn)),'FontSize',5,'color','k','rotation',90);
%             set(gca,'FontSize',7,'FontWeight','Normal','linewidth',0.75,'Ydir','normal','TickDir','out');
%             if rr == 1 & cc == 1
%                 text(-11,15,'Rasters - Occupancy normalized firing rate (FR)','FontSize',6.5);
%             end
%             if cc == 3 && rr == 1
%                 pos = get(gca,'Position');
%                 h = axes('Position',pos);
%                 changePosition(h,[0.01 0.27 0.02 -0.1]);
%                 hc = colorbar('north'); axis off
%                 set(hc,'YTick',[],'linewidth',0.01);
%                 changePosition(hc,[0 0 0 -0.01]);
%                 ylims = ylim;
%                 xlims = xlim;
%                 text(xlims(1)+0.15,ylims(1)+0.75,'0','FontSize',4.5);
%                 text(xlims(2)-0.4,ylims(1)+0.75,'Max FR','FontSize',4.5);
%             end
%         end
%     end
%     gg = gg + 1;
% end
figure(105);colormap parula;
fileName = fullfile(mData.pdf_folder,sprintf('rasters.pdf'));
save_pdf(gcf,mData.pdf_folder,sprintf('rasters.pdf'),600);
close(gcf);
winopen(fileName);
return;
end


%%
runthis = 0
if runthis
% ccsi = [19 85 127 129 126 128];
numberOfRows = 2;
numberOfCols = 3;
graphsOnOneFigure = numberOfRows * numberOfCols;
numberOfData = length(ccsi);
numberOfGroups = ceil(numberOfData/graphsOnOneFigure);
numberOfFullGroups = floor(numberOfData/graphsOnOneFigure);
indices = NaN(1,(numberOfGroups*graphsOnOneFigure));
indices(1:numberOfData) = 1:numberOfData;
groupIndices = reshape(indices,numberOfRows,numberOfCols,numberOfGroups);
ff = makeFigureRowsCols(105,[0.5 0.5 4 1],'RowsCols',[numberOfRows numberOfCols],...
    'spaceRowsCols',[0.07 0.06],'rightUpShifts',[0.09 0.17],'widthHeightAdjustment',...
    [-80 -200]);
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
%                     changePosition(hx,[0 50 0]);                    
                end
            end
            if cc >1
                set(gca,'YTick',[]);
            end
            if cc == 1
                h = ylabel('Trials');
                changePosition(h,[3 0 0]);
                text
            end
%                 text(65,15,'Occupancy Normalized Spike Rate (Hz/sec)','Rotation',90,'FontWeight','Normal','FontSize',5);
            text(1,11.2,sprintf('Cell %d - Max FR %d',cn,round(max(thisRaster(:)))),'FontSize',5,'color','k');
            text(53,4,sprintf('MI = %.2f',A.SI(cn)),'FontSize',5,'color','k','rotation',90);
            set(gca,'FontSize',7,'FontWeight','Normal','linewidth',0.75,'Ydir','normal','TickDir','out');
            if rr == 1 & cc == 1
                text(-11,15,'Rasters - Occupancy normalized firing rate (FR)','FontSize',6.5);
            end
            if cc == 3 && rr == 1
                pos = get(gca,'Position');
                h = axes('Position',pos);
                changePosition(h,[0.01 0.27 0.02 -0.1]);
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
save_pdf(gcf,mData.pdf_folder,sprintf('rasters.pdf'),600);
return;
end

%%
runthis = 1;
if runthis
ff = makeFigureRowsCols(100,[1 5 5.4 1.5],'RowsCols',[(length(ccsi))+1 1],'spaceRowsCols',[-0.009 0.0009],...
    'rightUpShifts',[0.04 0.02],'widthHeightAdjustment',[-70 7]);
set(gcf,'Position',[10 8 3.5 1]);
set(gcf,'color','w');
spSigAll = tei.plane{pl}.tP.deconv.spSigAll;
caSigAll = tei.plane{pl}.tP.deconv.caSigAll;
signals = tei.plane{pl}.tP.signals;
ccs = find(tei.plane{pl}.tP.iscell(:,1));
lengthSigs = size(signals,2);
b = tei.b;
traceTime = b.ts(tei.plane{pl}.b.frames_f(1:lengthSigs));
for cc = 1:length(ccsi)
    tsp = spSigAll{ccsi(cc)}';
    alltsp(cc,:) = tsp;
    caSig = caSigAll{ccsi(cc)}';
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

xStart = traceTime(find(traceTime>3,1,'first'));
xEnd = traceTime(find(traceTime<275,1','last'));
xlims = [xStart xEnd];
timeLength = 0.0833*24;
air_puff_raw = b.air_puff_raw;
air_puff_raw(b.air_puff_raw < 0.987) = NaN;
air_puff_raw(b.air_puff_raw >= 0.987) = 1;
onsets = tei.plane{pl}.contexts(1).markers.air_onsets;
offsets = tei.plane{pl}.contexts(1).markers.air_offsets;
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
    changePosition(gca,[-0.035 -0.01 0.06 0]);
    if cc<length(ccsi)+1
        tsp = alltsp(cc,:);
%         tsp(tsp==0) = NaN;
        caSig = allcaSig(cc,:);
        sCaSig = allsCaSig(cc,:);
%         plot(traceTime,caSig);hold on;
        plot(traceTime,sCaSig,'b');hold on;
        plot(traceTime,tsp,'r','linewidth',0.25);
        ylims = ylim;
%         ylim([min(tsp) max(tsp)]);
%         plot(b.ts,air_puff_raw*(ylims(2)/3),'-','color','b','linewidth',0.25);hold on;
        ylim([minCaSig maxCaSig]);
        text(xStart+16.5,maxCaSig/2,sprintf('%.2d',ccsi(cc)),'FontSize',4.5);
    else
%         plot(b.ts,0.25*b.photo_sensor_raw+0.55,'color','r');
        plot(b.ts,b.fSpeed,'color','m');hold on;
%         plot(b.ts,air_puff_raw*(max(b.fSpeed)/2),'color','b','linewidth',0.25);
        plot(b.ts,photo_sensor_raw,'color',[0 0.7 0.3]);
%         plot(dpsr*max(b.fSpeed),'color','r');
        n = 0;
%         plot(b.ts,0.5*b.fSpeed/max(b.fSpeed),'color','m');
        ylim([-0.1 max(b.fSpeed)+10]);
        
        sdfxs = xEnd - 156; 
        sdfxe = sdfxs + timeLength; sdfys = max(b.fSpeed)/2; sdfye = sdfys + max(b.fSpeed)/4;
        plot([sdfxs sdfxs],[sdfys-10 sdfye-10],'k','linewidth',0.5);
        text(xStart+16.5,sdfys+5,sprintf('Speed'),'FontSize',4.5);
        text(sdfxs-(0.15*100),sdfye+4,sprintf('%d cm/sec',round(max(b.fSpeed)/4)),'fontsize',4.5);
%         text(sdfxs+(0.015*60),sdfye-15,sprintf('%d cm/sec',round(max(b.fSpeed)/4)),'fontsize',4.5);
    end
%     ylim([min(caSig) max(caSig)]);
    xlim(xlims);
    
    
    axis off;
%     ylabel(num2str(cc));
    text(-25,max(caSig)/4,num2str(cc),'rotation',90);
    if cc == 3
        NcaPerNtsp = max(tsp)/maxCaSig;
        tspbar = NcaPerNtsp * 100;
        sdfxs = xEnd-175; 
        sdfxe = sdfxs + timeLength; sdfys = maxCaSig/2.5; sdfye = sdfys + 100;
        plot([sdfxs sdfxe],[sdfys sdfys],'k','linewidth',0.5);
        plot([sdfxe sdfxe],[sdfys sdfye],'k','linewidth',0.5);
        text(sdfxs-(0.2*60),sdfys + 60,sprintf('%dsecs',ceil(timeLength)),'fontsize',4);
        text(sdfxe+(0.05*20),sdfys + 60,'100% DF/F (25 A.U. Firing Rate)','fontsize',4);
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
    annotation('rectangle',[BLx BLy aW aH],'facealpha',0.2,'linestyle','none','facecolor','k');
end
figure(100);
save_pdf(gcf,mData.pdf_folder,sprintf('traces.pdf'),600);
% return;
end
%%
% showCells(102,ei,ccsi);
h = figure(102);clf;plot(0,0);
% ccsi = [19 85 127 129];
cellList = [];
showCells(gca,tei,pl,cellList)
set(h,'units','inches');
set(h,'Position',[10 4 1.5 1.5]);
changePosition(gca,[-0.165 -0.167 0.267 0.27]);
save_pdf(h,mData.pdf_folder,sprintf('cellsAvgImg.pdf'),600);

