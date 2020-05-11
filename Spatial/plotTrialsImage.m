function plotTrialsImage(b,markers1,markers2,fn)
%%
n = 0;
%%
ei = evalin('base','ei10');
mData = evalin('base','mData');
ii = 8; cc = 1;
b = ei{ii}.b;
markers1i = ei{ii}.plane{1}.contexts(cc).markers.air_onsets;
markers2i = ei{ii}.plane{1}.contexts(cc).markers.air_offsets;
fn = 101;

speed = b.fSpeed;
% figure(1000);clf;plot(b.speed);
ts = b.ts;

timeBefore = 0;
timeAfter = 15;

markers1 = markers1i - round(1e6 * timeBefore/b.si);
markers2 = markers2i + round(1e6 * timeAfter/b.si);

for ii = 1:length(markers1)
    st = markers1(ii);
    se = markers2(ii);
    sp{ii} = speed(st:se);
    lsp(ii) = length(sp{ii});
    t{ii} = ts(st:se)-ts(st);
    ind(ii) = find((st:se)-markers2i(ii)>0,1,'first');
end

mlsp = min(lsp);
for ii = 1:length(markers1)
    isp(ii,:) = sp{ii}(1:mlsp);
end
it = t{1}(1:mlsp);

% ff = makeFigureWindow__one_axes_only(5,[6 4 5 5],[0.1 0.012 0.85 0.99]);
figure(101);clf;
set(gcf,'units','inches');
set(gcf,'Position',[15 4 1 1]);
set(gcf,'color','w');

imagesc(isp);

% hc = colorbar;changePosition(hc,[0.19 -0.1 -0.02 0.02]);
% set(hc,'linewidth',0.25,'TickDirection','out');
changePosition(gca,[-0.01 -0.2 -0.2 0.07]);
box off
xtt = 0:5:20;
for ii = 1:length(xtt)
    xt(ii) = find(it>xtt(ii),1,'first');
end
txt = cellstr(num2str(ceil(xtt)'));

% axis equal
set(gca,'Ydir','normal','XTick',xt,'XTickLabel',txt,'FontSize',6,'FontWeight','Bold','TickDir','out','linewidth',1);
h = ylabel('Trials'); changePosition(h,[20000 0 0]);
h = xlabel('Time (sec)');changePosition(h,[0 0.7 0]);
text(length(it)+50000,1,'Speed (cm/sec)','rotation',90,'FontSize',5);
% n = 0;
hold on;
for ii = 1:length(markers1)
    plot([ind(ii) ind(ii)],[ii-0.5 ii+0.5],'r','linewidth',0.75)
end
% [x1 y1] = ds2nfu(-8000,13.5);
% [x2 y2] = ds2nfu(0,11.75);
% annotation('textarrow',[x1 x2],[y1 y2],'String','Air onset','HeadLength',3,'HeadWidth',3,'FontSize',5)
% [x1 y1] = ds2nfu(ind(ii)+6000,13);
% [x2 y2] = ds2nfu(ind(ii),11.6);
% annotation('textarrow',[x1 x2],[y1 y2],'String','Air offset','HeadLength',3,'HeadWidth',3,'FontSize',5,'color','r')
hc = putColorBar(gca,[-0.1 0.1 0 -0.2],[min(isp(:)) max(isp(:))],5);
colormap parula;

save_pdf(gcf,mData.pdf_folder,sprintf('speedTrialsImage_%d.pdf',cc),600);
return;

%%
raster = getDistRaster_forSpeed(b,markers1i,markers2i+round(1e6 * 0/b.si));
figure(101);clf;
set(gcf,'units','inches');
set(gcf,'Position',[25 4 4.5 1.5]);
plot(b.ts,b.air_puff_raw,'color','b','linewidth',1.5);hold on;
% plot(b.ts,0.25*b.photo_sensor_raw+0.55,'color','r','linewidth',1.5);
% plot(b.ts,0.5*b.fSpeed/max(b.fSpeed),'color','m','linewidth',1.5);
xlabel('Time (sec)');
ylim([0-0.1 1+0.1]);
xlim([0 100]);
save2pdf(sprintf('trialStructure_%d.pdf',cc),gcf,600);
% return;
% imagesc(raster.distSpeedRaster);
% hc = colorbar;changePosition(hc,[0.25 -0.02 -0.03 -0.03]);
% set(hc,'linewidth',0.25,'TickDirection','out');
% changePosition(gca,[0.035 -0.02 -0.15 -0.03]);
% box off
% % axis equal
% set(gca,'Ydir','normal','FontSize',8,'FontWeight','Bold','TickDir','out');
% h = ylabel('Trials'); changePosition(h,[8000 0 0]);
% h = xlabel('Time (sec)');changePosition(h,[0 0.1 0]);
% text(length(it)+30000,3,'Speed (cm/sec)','rotation',90,'FontSize',7);
% save2pdf(sprintf('speedRaster_%d.pdf',cc),gcf,600);
% n = 0;


function out =  getDistRaster_forSpeed(b,onsets,offsets)
trialNumbers = 1:length(onsets);
for iii = 1:length(trialNumbers)
    ii = trialNumbers(iii);
    st = onsets(ii);
    se = offsets(ii);
%     frames = (find(b.frames_f >= st & b.frames_f <= se)); % find frame numbers
%     oSignals{iii} = caSig(frames);
%     spoSignals{iii} = spSignal(frames);
%     maxSignals(iii) = max(oSignals{iii});
%     minSignals(iii) = min(oSignals{iii});
%     ts{iii} = b.ts(b.frames_f(frames))-b.ts(st);
    temp = b.fSpeed(st:se);
    temp(temp<0) = 0;
    speed{iii} = temp;
    tss{iii} = b.ts(st:se)-b.ts(st);
    dist{iii} = b.dist(st:se)-b.dist(st);
    distVal(iii) = dist{iii}(end);
%     if isfield(b,'tone_light_stim')
%         tempTL = (find(b.tone_light_stim_r >= st & b.tone_light_stim_r <= se)); % find frame numbers
%         try
%             tone_light{iii} = tempTL(3:4);
%         catch
%             tone_light{iii} = tempTL(1:2);
%         end
%         tstl{iii} = b.ts(b.tone_light_stim_r(tone_light{iii}))-b.ts(st);
%         dstl(iii,:) = b.dist(b.tone_light_stim_r(tone_light{iii}))-b.dist(st);
%     end
end
minDist = min(distVal);
nbins = 50;
bins = linspace(0,minDist,(nbins+1));
sbins = bins(1:nbins);
ebins = bins(2:end);

%%% Method %%%
% 1) Divide the minDist into nbins

for iii = 1:length(trialNumbers)
    thisDist = dist{iii};
    thisT = tss{iii};
%     thisTSig = ts{iii};
%     thisspSign = spoSignals{iii};
    thisSpeed = speed{iii};
    for jj = 1:nbins
        st = sbins(jj);  se = ebins(jj);
        dVals = find(thisDist >= st & thisDist <se);
        tVals = thisT(dVals);
        stt = tVals(1);
        sett = tVals(end);
%         tValsInd = find(thisTSig >= stt & thisTSig < sett);
%         if isnan(mean(thisspSign(tValsInd))/(sett-stt))
%             error('NaN in raster calculation');
%         end
%         distSigRaster(iii,jj) = mean(thisspSign(tValsInd));
%         distSpikesRaster(iii,jj) = sum(thisspSign(tValsInd))*(sett-stt);
%         distSigRaster(iii,jj) = mean(thisspSign(tValsInd))/(sett-stt);
        distDurRaster(iii,jj) = sett-stt;
        distSpeedRaster(iii,jj) = mean(thisSpeed(dVals));
    end
%     dSigF(iii,:) = applyGaussFilt(distSigRaster(iii,:),3);
end


out.distSpeedRaster = distSpeedRaster;
out.distDurRaster = distDurRaster;
out.dist = sbins;
