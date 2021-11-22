function tracePlots_1_C (ei)
%%

ei = evalin('base','ei');
% data = evalin('base','data');
% dataB = evalin('base','datab');
% dataAOn = evalin('base','dataAOn010');
% dataAOff= evalin('base','dataAOff010');
mData = evalin('base','mData');
selAnimals = 4; pl = 1;
tei = ei{selAnimals};
planeNumbers = pl;
maxDistTime = [Inf Inf];
conditionNumber = 3;
contextNumber = [conditionNumber conditionNumber];
stimMarkers = {'air','air'};
rasterTypes = {'dist','time'};
n = 0;
%%
varName = '';
selCells = 'areCellsAll';
for ss = 1:length(stimMarkers)
    distD = [];
    for jj = 1:length(selAnimals)
        [tempVals cns ACs] = getParamValues(varName,ei(selAnimals(jj)),planeNumbers,contextNumber(ss),stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
%         [clus] = getParamValues('cluster4',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),stimMarkers{ss},'dist',selCells,maxDistTime);
        [placeCells] = getParamValues('placeCells5',ei(selAnimals(jj)),planeNumbers,contextNumber(ss),stimMarkers{ss},'dist',selCells,maxDistTime);
    end
    dataC{ss} = tempVals;
end
n = 0;
A = dataC{1}; B = dataC{2}; %C = dataC{3}; %D = dataC{4};
% rasters = A;
placeCellNums = find(ACs);
PCs = placeCells;
% PCs = ACs;% & clus;
coiSI = find(PCs);
sum(PCs)/length(PCs);
[~,inds] = sort([A.place_field_properties.centers(coiSI)']);
coiSI = coiSI(inds); 
centers = A.place_field_properties.centers(coiSI);
[coiSI(inds) A.place_field_properties.centers(coiSI(inds))' A.place_field_properties.pws(coiSI(inds))'];

% coiSI =  coiSI(centers > 25 & centers < 80)';
% coiSI =  coiSI(centers > 40)';

% plot_trial_correlation(dataC,coiSI);return;
% plotRasters(dataC,coiSI);return
% ccsi = [132 374 59 18 181];% 51]; % animal 1 plane 1
ccsi = [168 356 662 406 157 75 66];% 51]; % animal 1 plane 1
ccsi = [231 586 66 282 81 146];% 51]; % animal 3 plane 1
ccsi = coiSI;
n=0;

%%
if 0
% showCells(102,ei,ccsi);
% close(figure(102))
h = figure(102);clf;plot(0,0);
% ccsi = [19 85 127 129];
% ccsi = [48 93 77 18 4]; % animal 1 plane 1
cellList = ccsi;
showCells(gca,tei,pl,cellList,[10 1])
set(h,'units','inches');
set(h,'Position',[10 4 1.5 1.5]);
changePosition(gca,[-0.165 -0.167 0.267 0.27]);
save_pdf(h,mData.pdf_folder,sprintf('cellsAvgImg.pdf'),600);
return;
end
%%
runthis = 0;
if runthis
cellList = ccsi;
numberOfRows = 1;
numberOfCols = 5;
graphsOnOneFigure = numberOfRows * numberOfCols;
numberOfData = length(ccsi);
numberOfGroups = ceil(numberOfData/graphsOnOneFigure);
numberOfFullGroups = floor(numberOfData/graphsOnOneFigure);
indices = NaN(1,(numberOfGroups*graphsOnOneFigure));
indices(1:numberOfData) = 1:numberOfData;
groupIndices = reshape(indices,numberOfRows,numberOfCols,numberOfGroups);
ff = makeFigureRowsCols(105,[0.5 0.5 4 1],'RowsCols',[numberOfRows numberOfCols],...
    'spaceRowsCols',[0.15 0.03],'rightUpShifts',[0.052 0.2],'widthHeightAdjustment',...
    [-40 -375]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[10 4 4.85 1.25]);
ddtt = 1;
for cc = 1:5
    cn = cellList(cc);
    A = dataC{ddtt};
    thisRaster = A.rasters(:,:,cn);
    minRasters(cc) = min(thisRaster(:));
    maxRasters(cc) = max(thisRaster(:));
end
minmin = min(minRasters(:));
maxmax = max(maxRasters(:));
for rr = 1
    for cc = 1:5
        cn = cellList(cc);
        axes(ff.h_axes(rr,cc));
        A = dataC{ddtt};
%         thisRaster = A.rasters(:,:,cn);
        thisRaster = A.sp_rasters1(:,:,cn);
        mSig = nanmean(thisRaster);
%         ft = fittype(A.formula);
        xs = 1:size(A.rasters,2);
%         coeff = A.coeff(:,cn);
        fitplot = gauss_fit(xs,A.gauss_fit_on_mean.coefficients_Rs_mean(cn,1:3),A.gauss_fit_on_mean.gauss1Formula);
        imagesc(thisRaster,[min(thisRaster(:)) max(thisRaster(:))]);hold on;
        if rr ==1 
            plot(size(thisRaster,1)*fitplot/max(fitplot),'linewidth',1,'color','m');
        end
        box off;
        if rr == 1
        text(size(thisRaster,2)+size(thisRaster,2)/20,0,sprintf('Max FR %d - zMI = %.2f - Rs = %.2f',round(max(thisRaster(:))),A.SI(cn),...
            A.gauss_fit_on_mean.coefficients_Rs_mean(cn,4)),'FontSize',4,'color','k','rotation',90);
        else
       text(size(thisRaster,2)+size(thisRaster,2)/15,2,sprintf('Max FR %d - zMI = %.2f',round(max(thisRaster(:))),A.SI(cn)),'FontSize',4,'color','k','rotation',90);
        end
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
            xticks = [1:10:size(thisRaster,2)];
            xs = 0:0.2:100;
            set(gca,'XTick',xticks,'XTickLabels',xs(xticks));
            for bb = 1:length(A.lastBin)
                xvalbin = A.lastBin(bb)-1;
                plot([xvalbin xvalbin]+0.15,[bb-0.75 bb+0.75],'r','linewidth',1.5);
            end
            
                hx = xlabel('Time (sec)');
%                     changePosition(hx,[-20 0 0]);                    
        end

        if cc >1
            set(gca,'YTick',[]);
        end
        if cc == 1
            set(gca,'YTick',[1 5 10]);
            h = ylabel('Trials');
%             changePosition(h,[3 0 0]);
%             text
        end
        
        if cc == 5 && rr == 1
            hca = gca;
            hc = putColorBar(hca,[0 0 -0.05 0],{'0','Max FR'},6,'northoutside',[0.15 0.22 0.05 0.22]);
        end
    end
end
figure(105);colormap parula;
save_pdf(gcf,mData.pdf_folder,sprintf('rasters_dist.pdf'),600);
return;
end

%%
runthis = 0;
if runthis
cellList = ccsi;
numberOfRows = 1;
numberOfCols = 5;
graphsOnOneFigure = numberOfRows * numberOfCols;
numberOfData = length(ccsi);
numberOfGroups = ceil(numberOfData/graphsOnOneFigure);
numberOfFullGroups = floor(numberOfData/graphsOnOneFigure);
indices = NaN(1,(numberOfGroups*graphsOnOneFigure));
indices(1:numberOfData) = 1:numberOfData;
groupIndices = reshape(indices,numberOfRows,numberOfCols,numberOfGroups);
ff = makeFigureRowsCols(105,[0.5 0.5 4 1],'RowsCols',[numberOfRows numberOfCols],...
    'spaceRowsCols',[0.15 0.03],'rightUpShifts',[0.052 0.2],'widthHeightAdjustment',...
    [-40 -375]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[10 4 4.85 1.25]);
ddtt = 2;
for cc = 1:5
    cn = cellList(cc);
    A = dataC{ddtt};
    thisRaster = A.rasters(:,:,cn);
    minRasters(cc) = min(thisRaster(:));
    maxRasters(cc) = max(thisRaster(:));
end
minmin = min(minRasters(:));
maxmax = max(maxRasters(:));
for rr = 2
    for cc = 1:5
        cn = cellList(cc);
        axes(ff.h_axes(1,cc));
        A = dataC{ddtt};
%         thisRaster = A.rasters(:,:,cn);
        thisRaster = A.sp_rasters1(:,:,cn);
        mSig = nanmean(thisRaster);
%         ft = fittype(A.formula);
        xs = 1:size(A.rasters,2);
%         coeff = A.coeff(:,cn);
        fitplot = gauss_fit(xs,A.gauss_fit_on_mean.coefficients_Rs_mean(cn,1:3),A.gauss_fit_on_mean.gauss1Formula);
        imagesc(thisRaster,[min(thisRaster(:)) max(thisRaster(:))]);hold on;
        if rr ==1 
            plot(size(thisRaster,1)*fitplot/max(fitplot),'linewidth',1,'color','m');
        end
        box off;
        if rr == 1
            text(size(thisRaster,2)+size(thisRaster,2)/20,0,sprintf('Max FR %d - zMI = %.2f - Rs = %.2f',round(max(thisRaster(:))),A.SI(cn),...
            A.gauss_fit_on_mean.coefficients_Rs_mean(cn,4)),'FontSize',4,'color','k','rotation',90);
        else
        text(size(thisRaster,2)+size(thisRaster,2)/15,2,sprintf('Max FR %d - zMI = %.2f',round(max(thisRaster(:))),A.SI(cn)),'FontSize',4,'color','k','rotation',90);
        end
%         if rr == 1
            text(size(thisRaster,1)/2,size(thisRaster,1)+1.25,sprintf('Cell %d',cn),'FontSize',6,'color','k');
%         end
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
            xticks = [1:10:size(thisRaster,2)];
            xs = 0:0.2:100;
            set(gca,'XTick',xticks,'XTickLabels',xs(xticks));
            for bb = 1:length(A.lastBin)
                xvalbin = A.lastBin(bb)-1;
                plot([xvalbin xvalbin]+0.15,[bb-0.75 bb+0.75],'r','linewidth',1);
            end
            
                hx = xlabel('Time (sec)');
%                     changePosition(hx,[-20 0 0]);                    
        end

        if cc >1
            set(gca,'YTick',[]);
        end
        if cc == 1
            set(gca,'YTick',[1 5 10]);
            h = ylabel('Trials');
%             changePosition(h,[3 0 0]);
%             text
        end
        
        if cc == 5 %&& rr == 1
            hca = gca;
            hc = putColorBar(hca,[0 0 -0.05 0],{'0','Max FR'},6,'northoutside',[0.15 0.22 0.05 0.22]);
        end
    end
end
figure(105);colormap parula;
save_pdf(gcf,mData.pdf_folder,sprintf('rasters_time.pdf'),600);
return;
end

%%
runthis = 0;
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
    'spaceRowsCols',[0.17 0.04],'rightUpShifts',[0.1 0.15],'widthHeightAdjustment',...
    [-60 -220]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[5 6 3.45 1.95]);
for rr = 1:2
    for cc = 1:5
        cn = cellList(cc);
        axes(ff.h_axes(rr,cc));
        A = dataC{rr};
        thisRaster = A.rasters(:,:,cn);
        minRasters(rr,cc) = min(thisRaster(:));
        maxRasters(rr,cc) = max(thisRaster(:));
    end
end
minmin = min(minRasters(:));
maxmax = max(maxRasters(:));
mean_color = [0.85 0.85 0.85];
for cc = 1:5
    cn = cellList(cc);
    axes(ff.h_axes(1,cc));
    A = dataC{2};
    thisRaster = A.sp_rasters1(:,:,cn);
    mSig = nanmean(thisRaster);
    xs = 1:size(A.rasters,2);
    fitplot = gauss_fit(xs,A.gauss_fit_on_mean.coefficients_Rs_mean(cn,1:3),A.gauss_fit_on_mean.gauss1Formula);
    imagesc(thisRaster,[min(thisRaster(:)) max(thisRaster(:))]);hold on;
    box off;
    plot(size(thisRaster,1)*mSig/max(mSig),'linewidth',0.25,'color',mean_color);
    text(size(thisRaster,2)+size(thisRaster,2)/15,0.5,sprintf('MFR = %d , zMI = %.2f',round(max(thisRaster(:))),A.SI(cn)),'FontSize',4,'color','k','rotation',90);
    text(size(thisRaster,1)/2,size(thisRaster,1)+1.5,sprintf('Cell %d',cn),'FontSize',6,'color','k');
    set(gca,'FontSize',5,'FontWeight','Normal','linewidth',0.75,'Ydir','normal','TickDir','out');
    xticks = [1:10:size(thisRaster,2)];
    xs = 0:0.2:100;
    set(gca,'XTick',xticks,'XTickLabels',xs(xticks));
    for bb = 1:length(A.lastBin)
        xvalbin = A.lastBin(bb)-1;
        plot([xvalbin xvalbin]+0.15,[bb-0.75 bb+0.75],'r','linewidth',0.5);
    end
    hx = xlabel('Time (sec)');
    if cc >1
        set(gca,'YTick',[]);
    end
    if cc == 1
        set(gca,'YTick',[1 5 10]);
        h = ylabel({'Air-Time (AT)','Trials'});
        text(-5,13.5,sprintf('Condition %d',contextNumber(1)),'FontSize',6,'color','k','FontWeight','Bold');
    end
    if cc == 5
        hca = gca;
        hc = putColorBar(hca,[0.01 0.01 -0.05 0],{'0','MFR'},5,'northoutside',[0.25 0.3 0.1 0.3]);
    end
    

    
    axes(ff.h_axes(2,cc));
    A = dataC{1};
    thisRaster = A.sp_rasters1(:,:,cn);
    mSig = nanmean(thisRaster);
    xs = 1:size(A.rasters,2);
    fitplot = gauss_fit(xs,A.gauss_fit_on_mean.coefficients_Rs_mean(cn,1:3),A.gauss_fit_on_mean.gauss1Formula);
    imagesc(thisRaster,[min(thisRaster(:)) max(thisRaster(:))]);hold on;
    box off;
    if  cc == 4
        plot(size(thisRaster,1)*fitplot/max(fitplot),'linewidth',1,'color','m');
        text(size(thisRaster,1)/2,size(thisRaster,1)+1.25,sprintf('Rs = %.2f',A.gauss_fit_on_mean.coefficients_Rs_mean(cn,4)),'FontSize',4,'color','k');
    end
    plot(size(thisRaster,1)*mSig/max(mSig),'linewidth',0.25,'color',mean_color);
%     text(size(thisRaster,2)+size(thisRaster,2)/20,0,sprintf('MFR %d - zMI = %.2f - Rs = %.2f',round(max(thisRaster(:))),A.SI(cn),...
%     A.gauss_fit_on_mean.coefficients_Rs_mean(cn,4)),'FontSize',4,'color','k','rotation',90);
    text(size(thisRaster,2)+size(thisRaster,2)/15,0.5,sprintf('MFR = %d , zMI = %.2f',round(max(thisRaster(:))),A.SI(cn)),'FontSize',4,'color','k','rotation',90);
    set(gca,'FontSize',5,'FontWeight','Normal','linewidth',0.75,'Ydir','normal','TickDir','out');
    xticks = [1  round(size(thisRaster,2)/2) size(thisRaster,2)];
    xsd = 0:3:1000; xsd = xsd(1:size(thisRaster,2));
    set(gca,'XTick',xticks,'XTickLabels',xsd(xticks));
    hx = xlabel('Position (cm)');
    if cc >1
        set(gca,'YTick',[]);
    end
    if cc == 1
        set(gca,'YTick',[1 5 10]);
        h = ylabel({'Air-Distance (AD)','Trials'});
    end
    
end

figure(105);colormap parula;
fileName = fullfile(mData.pdf_folder,sprintf('rasters.pdf'));
save_pdf(gcf,mData.pdf_folder,sprintf('rasters_%d.pdf',contextNumber(1)),600);
return;
end

%%
runthis = 1;
if runthis
ff = makeFigureRowsCols(100,[1 5 5.4 1.5],'RowsCols',[(length(ccsi))+1 1],'spaceRowsCols',[-0.009 0.0009],...
    'rightUpShifts',[0.04 0.02],'widthHeightAdjustment',[-70 7]);
set(gcf,'Position',[10 8 4.8 1]);
set(gcf,'color','w');
spSigAll = tei.plane{pl}.tP.deconv.spSigAll;
% caSigAll = tei.plane{pl}.tP.deconv.caSigAll;
signals = get_calcium_data(tei,pl);
caSigAll = signals;
ccs = find(tei.plane{pl}.tP.iscell(:,1));
lengthSigs = size(signals,2);
b = tei.b;
traceTime = b.ts(tei.plane{pl}.b.frames_f(1:lengthSigs));
for cc = 1:length(ccsi)
    tsp = spSigAll(ccsi(cc),:);
    alltsp(cc,:) = tsp;
    caSig1 = caSigAll(ccsi(cc),:);
    caSig = signals(placeCellNums(ccsi(cc)),:);
    allcaSig(cc,:) = caSig;
    sCaSig = caSigAll(ccsi(cc),:);
    allsCaSig(cc,:) = sCaSig;
    minSCaSig(cc) = min(sCaSig); maxSCaSig(cc) = max(sCaSig);
    mintsp(cc) = min(tsp); maxtsp(cc) = max(tsp);
end
minCaSig = min(minSCaSig);
maxCaSig = max(maxSCaSig);

% startTimeFrame = b.frames_f(find(traceTime>300,1,'first'));
% endTimeFrame = b.frames_f(find(traceTime<550,1','last'));
condN = conditionNumber;
onsets = tei.plane{pl}.contexts(condN).markers.air_onsets;
% onsets1 = tei.plane{pl}.contexts(condN+1).markers.air_onsets;
offsets = tei.plane{pl}.contexts(condN).markers.air_offsets;
% offsets1 = tei.plane{pl}.contexts(condN+1).markers.air_offsets;
xStart = traceTime(find(traceTime>(b.ts(onsets(1))-5),1,'first'));
xEnd = traceTime(find(traceTime<(b.ts(offsets(end))+5),1','last'));
xlims = [xStart xEnd];
timeLength = 0.0833*24;
air_puff_raw = b.air_puff_raw;
air_puff_raw(b.air_puff_raw < 0.987) = NaN;
air_puff_raw(b.air_puff_raw >= 0.987) = 1;

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
        plot(traceTime,caSig,'b');hold on;
        plot(traceTime,tsp,'r','linewidth',0.25);
        ylims = ylim;
%         ylim([min(tsp) max(tsp)]);
%         plot(b.ts,air_puff_raw*(ylims(2)/3),'-','color','b','linewidth',0.25);hold on;
        ylim([minCaSig maxCaSig]);
%         ylim([min(caSig) max(caSig)]);
        text(xStart+21.5,maxCaSig/2,sprintf('%.2d',ccsi(cc)),'FontSize',4.5);
    else
%         plot(b.ts,0.25*b.photo_sensor_raw+0.55,'color','r');
        plot(b.ts,b.fSpeed,'color','m');hold on;
%         plot(b.ts,air_puff_raw*(max(b.fSpeed)/2),'color','b','linewidth',0.25);
        plot(b.ts,photo_sensor_raw,'color',[0 0.7 0.3]);
%         plot(dpsr*max(b.fSpeed),'color','r');
        n = 0;
%         plot(b.ts,0.5*b.fSpeed/max(b.fSpeed),'color','m');
        ylim([-0.1 max(b.fSpeed)+10]);
        
        sdfxs = xEnd - 136; 
        sdfxe = sdfxs + timeLength; sdfys = max(b.fSpeed)/2; sdfye = sdfys + max(b.fSpeed)/4;
        plot([sdfxs sdfxs],[sdfys-10 sdfye-10],'k','linewidth',0.5);
        text(xStart+21.5,sdfys+10,sprintf('Speed'),'FontSize',4.5);
        text(sdfxs-(0.01*100),sdfye+4,sprintf('%d cm/sec',round(max(b.fSpeed)/4)),'fontsize',4.5);
%         text(sdfxs+(0.015*60),sdfye-15,sprintf('%d cm/sec',round(max(b.fSpeed)/4)),'fontsize',4.5);
    end
%     ylim([min(caSig) max(caSig)]);
    xlim(xlims);
    
    
    axis off;
%     ylabel(num2str(cc));
    text(-25,max(caSig)/4,num2str(cc),'rotation',90);
    if cc == 2
        NcaPerNtsp = max(tsp)/maxCaSig;
        tspbar = NcaPerNtsp * 100;
        sdfxs = xEnd-125; 
        sdfxe = sdfxs + timeLength; sdfys = maxCaSig/2.5; sdfye = sdfys + 100;
        plot([sdfxs sdfxe],[sdfys sdfys],'k','linewidth',0.5);
        plot([sdfxe sdfxe],[sdfys sdfye],'k','linewidth',0.5);
        text(sdfxs-(0.15*60),sdfys + 60,sprintf('%dsecs',ceil(timeLength)),'fontsize',4);
        text(sdfxe+(0.05*20),sdfys + 60,'100% DF/F (25 A.U. Firing Rate)','fontsize',4);
    end
    n = 0;
    cc = cc + 1;
end
axes(ff.h_axes(1,1));ylims = ylim;
[TLx TLy] = ds2nfu(b.ts(onsets(1)),ylims(2)-0);
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
save_pdf(gcf,mData.pdf_folder,sprintf('traces_cond_%d.pdf',condN),600);
% return;
end

%%
runthis = 0;
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
    caSig1 = caSigAll{ccsi(cc)}';
    caSig = signals(ccsi(cc),:);
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
xEnd = traceTime(find(traceTime<210,1','last'));
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
        text(xStart+78.5,maxCaSig/2,sprintf('%.2d',ccsi(cc)),'FontSize',4.5);
    else
%         plot(b.ts,0.25*b.photo_sensor_raw+0.55,'color','r');
        plot(b.ts,b.fSpeed,'color','m');hold on;
%         plot(b.ts,air_puff_raw*(max(b.fSpeed)/2),'color','b','linewidth',0.25);
        plot(b.ts,photo_sensor_raw,'color',[0 0.7 0.3]);
%         plot(dpsr*max(b.fSpeed),'color','r');
        n = 0;
%         plot(b.ts,0.5*b.fSpeed/max(b.fSpeed),'color','m');
        ylim([-0.1 max(b.fSpeed)+10]);
        
        sdfxs = xEnd - 146; 
        sdfxe = sdfxs + timeLength; sdfys = max(b.fSpeed)/2; sdfye = sdfys + max(b.fSpeed)/4;
        plot([sdfxs sdfxs],[sdfys-10 sdfye-10],'k','linewidth',0.5);
        text(xStart+28.5,sdfys+10,sprintf('Speed'),'FontSize',4.5);
        text(sdfxs-(0.15*100),sdfye+4,sprintf('%d cm/sec',round(max(b.fSpeed)/4)),'fontsize',4.5);
%         text(sdfxs+(0.015*60),sdfye-15,sprintf('%d cm/sec',round(max(b.fSpeed)/4)),'fontsize',4.5);
    end
%     ylim([min(caSig) max(caSig)]);
    xlim(xlims);
    
    
    axis off;
%     ylabel(num2str(cc));
    text(-25,max(caSig)/4,num2str(cc),'rotation',90);
    if cc == 2
        NcaPerNtsp = max(tsp)/maxCaSig;
        tspbar = NcaPerNtsp * 100;
        sdfxs = xEnd-155; 
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
[TLx TLy] = ds2nfu(b.ts(onsets(1)),ylims(2)-0);
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


