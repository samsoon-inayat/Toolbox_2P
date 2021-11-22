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
selContexts = [3 4 5]; rasterNames = {'airD','airD','airD'};
oD = get_data(ei,selContexts,rasterNames);

R = oD.Rs{selAnimals,conditionNumber};
[rs,MFR,centersA,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width); 
n = 0;

%%
props1 = get_props_Rs(oD.Rs,70);
cLi = (props1.good_FR_and_Gauss{selAnimals,conditionNumber});
cL = find(cLi);
centers = centersA(cLi);
[~,inds] = sort([centersA(cLi)']);
cL = cL(inds); 
% cL'
cL(centers > 0 & centers < 30)'
cL(centers > 30 & centers < 60)'
cL(centers > 60 & centers < 90)'
cL(centers > 90 & centers < 120)'
cL(centers > 120 & centers < 150)'
ccsi = [97 203 189 180 270];
% plotRasters_simplest(R,cL)
% return
n=0;


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
    caSig = signals((ccsi(cc)),:);
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
%         ylim([minCaSig maxCaSig]);
        ylim([min(caSig) max(caSig)]);
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

