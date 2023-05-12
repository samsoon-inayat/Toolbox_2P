function tracePlots_1_C (ei)
%%
ei = evalin('base','ei10_C');

% selContexts = [1 2 3 4];
% rasterNames = {'airD','airD','airD','airD'};
% RsC = get_rasters_data(ei,selContexts,rasterNames);
% mRsC = calc_mean_rasters(RsC,1:10);
% RsC = find_responsive_rasters(RsC,1:10);
% % view_population_vector(Rs,mRs,300);
% [resp_fractionC,resp_valsC,OIC,mean_OIC] = get_responsive_fraction(RsC)

si = [C1_t_D C2_t_D C3_t_D C4_t_D];
RsC = oC.Rs(:,si);
ntrials = 50;
props_C = get_props_Rs(RsC,ntrials);
pop_var_name = {'vals','good_zMI'};
pop_var_name = {'vals','good_zMI','good_Gauss'};
sel_pop_C = cell_list_op(props_C,pop_var_name); 


mData = evalin('base','mData');
selAnimals = 3; pl = 1;
tei = ei{selAnimals};
planeNumbers = pl;
maxDistTime = [Inf Inf];
conditionNumber = 1;
contextNumber = [conditionNumber conditionNumber];
stimMarkers = {'air','air'};
rasterTypes = {'dist','time'};
n = 0;
%%

ccsi = [235 158 251 147 97 18];% 51]; % animal 3 plane 1

ccsi = [235 158 251 260 97 18];% 51]; % animal 3 plane 1

n = 0;
%%
spSigAll = tei.plane{pl}.tP.deconv.spSigAll;
T_C2 = evalin('base','T_C');
temp = get_cluster_timecourses(cell2mat(T_C2{selAnimals,7}))
caSigAll = temp.ratio';
n = 0;
%%
alltsp = []; allcaSig = []; allsCaSig = [];

signals = caSigAll;%tei.plane{pl}.tP.signals;
ccs = find(tei.plane{pl}.tP.iscell(:,1));
lengthSigs = size(signals,2);
b = tei.b;
traceTime = b.ts(tei.plane{pl}.b.frames_f(1:lengthSigs));
for cc = 1:length(ccsi)
    tsp = spSigAll(ccsi(cc),:)';
    alltsp(cc,:) = tsp;
    caSig1 = caSigAll(ccsi(cc),:)';
    caSig = signals(ccsi(cc),:);
    allcaSig(cc,:) = caSig;
    sCaSig = caSigAll(ccsi(cc),:)';
    allsCaSig(cc,:) = sCaSig;
    minSCaSig(cc) = min(sCaSig); maxSCaSig(cc) = max(sCaSig);
    mintsp(cc) = min(tsp); maxtsp(cc) = max(tsp);
end
minCaSig = min(minSCaSig);
maxCaSig = max(maxSCaSig);
minCaC = minCaSig;
maxCaC = maxCaSig;

ff = makeFigureRowsCols(100,[1 5 5.4 1.5],'RowsCols',[(length(ccsi))+1 1],'spaceRowsCols',[-0.009 0.0009],...
    'rightUpShifts',[0.04 0.02],'widthHeightAdjustment',[-70 7]);
set(gcf,'Position',[10 5 4.8 1]);
set(gcf,'color','w');

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
        plot(traceTime,caSig,'b','linewidth',0.25);hold on;
        plot(traceTime,2*tsp,'r','linewidth',0.25);
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
        text(sdfxe+(0.05*20),sdfys + 60,'100% DF/F (50 A.U. Firing Rate)','fontsize',4);
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


%% raw calcium traces and deconvolved spikes
% hf = figure(100);clf;set(gcf,'Units','Inches');set(gcf,'Position',[1 5 6.97 5.5],'color','w'); hold on;
ff = makeFigureRowsCols(100,[10 4 6.97 5],'RowsCols',[7 1],'spaceRowsCols',[0.04 0.02],'rightUpShifts',[0.07 0.08],'widthHeightAdjustment',[-140 -50]);
selcells = [11 66  3  205   255   196];
% selcells = [11 66  3  205   233   75];
for ii = 1:6
    axes(ff.h_axes(ii,1));
    yyaxis left
    plot(traceTime,signalsR(selcells(ii),:)','b');hold on;
    yyaxis right
    plot(traceTime,signalsS(selcells(ii),:)','r');
    ylims = ylim; ylim([0 ylims(2)]);
    yyaxis left
    ylims = ylim;
%     for ii = 1:length(onsets)
%        tt = b.ts(onsets(ii));
%        plot([tt tt],[ylims(1) ylims(2)],'color',colors{2},'linewidth',lwdth,'marker','none');
%        tt = b.ts(offsets(ii));
%        plot([tt tt],[ylims(1) ylims(2)],'color',colors{3},'linewidth',lwdth,'marker','none');
%    end
    % changePosition(gca,[-0.075 -0.005 0.16 -0.07]);
    xlim([1.35 4.05]);
    xlim([0 traceTime(end)]);
    box off
    if ii >2
        ha = gca;
%         ha.XAxis.Visible = 'off';
    end
    format_axes(gca);
    ylabel('DF/F (%)');
    ylabel('Raw Ca Sig');
    yyaxis right;
    ylabel('FR (A.U.)');
end
[xoff,yoff] = get_ca_motion_correction_data(pd_rec,pl);
spSigAll = sqrt(xoff.^2 + yoff.^2) * pixelSize(an);

axes(ff.h_axes(ii+1,1));
plot(traceTime,spSigAll','m');hold on;
xlim([1.35 4.05]);
xlim([0 traceTime(end)]);
format_axes(gca);
ylabel('MC Dist (\mum)');
xlabel('Time (min)');

axes(ff.h_axes(1,1));ylims = ylim;
[TLx TLy] = ds2nfu(b.ts(onsets(1)),ylims(2)-0);
axes(ff.h_axes(7,1));ylims = ylim;
[BLx BLy] = ds2nfu(b.ts(onsets(1)),ylims(1));
aH = (TLy - BLy);
% for ii = 1:length(onsets)
for ii = [1:10 41:50]%length(onsets)
    [BRx BRy] = ds2nfu(b.ts(offsets(ii)),ylims(1));
    [BLx BLy] = ds2nfu(b.ts(onsets(ii)),ylims(1));
    aW = (BRx-BLx);
    annotation('rectangle',[BLx BLy aW aH],'facealpha',0.2,'linestyle','none','facecolor','k');
end

save_pdf(ff.hf,mData.pdf_folder,sprintf('overall_pop.pdf'),600);
