function trace_plot_all(pd_rec,flag)
if ~exist('flag','var')
    flag = [1 1 0 1];
end
if ~exist('pd_rec','var')
    pd_rec = evalin('base','ei{4}');
end
an = 4;
pl = 1;
try
    signalsS = pd_rec.deconv';
    frames_f = pd_rec.b.frames_f;
catch
    signalsS = pd_rec.plane{pl}.tP.deconv.spSigAll;
    frames_f = pd_rec.plane{pl}.b.frames_f;
end

b = pd_rec.b;
b.ts = b.ts/60;
traceTime = b.ts(frames_f);

% onsets = b.photo_sensor_f;
% offsets = b.photo_sensor_r;
% st = find(traceTime>(b.ts(onsets(1))-5),1,'first'); en = find(traceTime<(b.ts(offsets(end))+5),1','last');
% xStart = traceTime(st);
% xEnd = traceTime(en);
% xlims = [xStart xEnd];
onsets = b.air_puff_r;
offsets = b.air_puff_f;
light_onsets = b.stim_r;
mData = evalin('base','mData');
signals = get_calcium_data(pd_rec,pl);
signals = signalsS;
% [signals,signalsR] = get_calcium_data_raw(pd_rec,pl);
% [xoff,yoff] = get_ca_motion_correction_data(pd_rec,pl);
% r = sqrt(xoff.^2 + yoff.^2);
% figure(33);
% plot(traceTime,r);

% for ii = 1:5
%     pixelSize(ii) = ei{ii}.thorExp.widthUM/ei{ii}.thorExp.pixelX;
% end
n = 0;
%%
while 0
    hf = figure(100);clf;set(gcf,'Units','Inches');set(gcf,'Position',[1 5 6.95 3],'color','w'); hold on;
    spSigAllN = normalizeSignal(signals,2);
    [maxVal,maxLoc] = max(spSigAllN,[],2);
    [sorted,locs] = sort(maxLoc);
    maxSig = 0.5;
    hi = imagesc([traceTime(1), traceTime(end)],[1 size(spSigAllN,1)],spSigAllN(locs,:),[0 maxSig]);
    % colorbar;
    hold on;
    % plot(b.ts(frames_f),spSigAll(:,:));
    ylims = ylim;
    hold on;
    lwdth = 0.1;
    colors = {[1 0 1]/1.5;[0 0.3 1];[1 1 0]/1.75;[0 1 1]/1.25};
    if flag(1)
    if isfield(b,'air_puff_raw')
    %    plot(b.ts,(b.air_puff_raw-0.5)*ylims(2)*4,'color',mData.colors{3},'linewidth',0.5);
       for ii = 1:length(onsets)
           tt = b.ts(onsets(ii));
           plot([tt tt],[ylims(1)-5 ylims(2)+5],'color',colors{2},'linewidth',lwdth);
           tt = b.ts(offsets(ii));
           plot([tt tt],[ylims(1)-5 ylims(2)+5],'color',colors{3},'linewidth',lwdth);
       end
    end

    end
    if flag(2)
    if isfield(b,'stim_raw')
    %    plot(b.ts,(b.stim_raw-0.5)*ylims(2)*4,'m','linewidth',0.5);
       for ii = 1:length(light_onsets)
           tt = b.ts(light_onsets(ii));
           plot([tt tt],[ylims(1)-5 ylims(2)+5],'color',colors{1},'linewidth',lwdth);
       end
    end
    end
    if flag(3)
    plot(b.ts,b.photo_sensor_raw*ylims(2),mData.colors{3},'linewidth',0.5)
    end
    if flag(4)
    plot(b.ts,b.fSpeed*0.25,'color',colors{4},'linewidth',lwdth);
    end
    numcells = 100;
    ylim([0 numcells]); xlim([0 traceTime(end)]);
    ylims = ylim;
    upfac = 3;
    text(b.ts(light_onsets(1)),ylims(2)+upfac,'L','FontSize',6);
    text(b.ts(light_onsets(21)),ylims(2)+upfac,'L*','FontSize',6);
    text(b.ts(onsets(1)),ylims(2)+upfac,'A','FontSize',6);
    text(b.ts(onsets(11)),ylims(2)+upfac,'R','FontSize',6);
    text(b.ts(onsets(21)),ylims(2)+upfac,'V','FontSize',6);
    text(b.ts(onsets(31)),ylims(2)+upfac,'R*','FontSize',6);
    text(b.ts(onsets(41)),ylims(2)+upfac,'A*','FontSize',6);
    text(b.ts(light_onsets(1))-1.45,ylims(2)+upfac+2,'Conditions','FontSize',6);
    set(gca,'Ydir','Normal');
    xlabel('Time (min)');
    ylabel('Cell Number (arranged by peaks)');
    set(gca,'FontSize',6,'FontWeight','Normal');
    changePosition(gca,[-0.075 -0.025 0.16 -0.01]);
    ht = title(sprintf('Normalized Calcium signal (%cF/Fo) of %d of %d cells from a representative animal',916,numcells,size(signals,1))); 
    changePosition(ht,[0 9 0]);
    set(ht,'FontSize',6,'FontWeight','Normal');
    legs = {'Light onset','Air onset','Air offset','Speed',[0.5 0.1 109 0.1]};
    putLegendH(gca,legs,'colors',colors);
    hc = putColorBar(gca,[0.75 0.62 -0.8 -0.6],{'0',sprintf('>%.1f',maxSig)},6,'northoutside',[0.15 0.27 0.045 0.27]);
%     colormap_ig
    colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('overall_pop.pdf'),600);
    break;
end

%%
%%
    CRc = corr(spSigAllN);
    figure(1000);clf;
    subplot(2,1,1);imagesc(spSigAllN);colorbar;set(gca,'Ydir','Normal');
    subplot(2,1,2);imagesc(CRc);colorbar;set(gca,'Ydir','Normal');
    
    figure(1003);clf;
    imagesc(CRc);colorbar;set(gca,'Ydir','Normal');
    %%
    CRc1 = corr(spSigAllN');
    oM = ones(size(CRc1));
    mask1 = (triu(oM,0) & tril(oM,0)); CRc1(mask1==1) = NaN;
    
    figure(1001);clf;
    subplot(2,1,1);imagesc(spSigAllN');colorbar;set(gca,'Ydir','Normal');
    subplot(2,1,2);imagesc(CRc1);colorbar;set(gca,'Ydir','Normal');
   
%%
while 1
    hf = figure(100);clf;set(gcf,'Units','Inches');set(gcf,'Position',[1 5 6.97 2.4],'color','w'); hold on;
    spSigAllN = normalizeSignal(signals,2);
    [maxVal,maxLoc] = max(spSigAllN,[],2);
    [sorted,locs] = sort(maxLoc);
    maxSig = 0.5;
    hi = imagesc([traceTime(1), traceTime(end)],[1 size(spSigAllN,1)],spSigAllN(locs,:),[0 maxSig]);
    % colorbar;
    hold on;
    % plot(b.ts(frames_f),spSigAll(:,:));
    ylims = ylim;
    hold on;
    lwdth = 0.1;
    colors = {[1 0 1]/1.5;[0 0.3 1];[1 1 0]/1.75;[0 1 1]/1.25};
    if flag(1)
    if isfield(b,'air_puff_raw')
    %    plot(b.ts,(b.air_puff_raw-0.5)*ylims(2)*4,'color',mData.colors{3},'linewidth',0.5);
       for ii = 1:length(onsets)
           tt = b.ts(onsets(ii));
           plot([tt tt],[ylims(1)-5 ylims(2)+5],'color',colors{2},'linewidth',lwdth);
           tt = b.ts(offsets(ii));
           plot([tt tt],[ylims(1)-5 ylims(2)+5],'color',colors{3},'linewidth',lwdth);
       end
    end

    end
    if flag(2)
    if isfield(b,'stim_raw')
    %    plot(b.ts,(b.stim_raw-0.5)*ylims(2)*4,'m','linewidth',0.5);
       for ii = 1:length(light_onsets)
           tt = b.ts(light_onsets(ii));
           plot([tt tt],[ylims(1)-5 ylims(2)+5],'color',colors{1},'linewidth',lwdth);
       end
    end
    end
    if flag(3)
    plot(b.ts,b.photo_sensor_raw*ylims(2),mData.colors{3},'linewidth',0.5)
    end
    if flag(4)
    plot(b.ts,b.fSpeed*0.25,'color',colors{4},'linewidth',lwdth);
    plot([traceTime(end) traceTime(end)],[0 0.25*30],'color',[0 0 1],'linewidth',lwdth+1.5);
    end
    numcells = 100;
    ylim([0 numcells+1]); xlim([0 traceTime(end)]);
    ylims = ylim; rtfac = 0.5;
    plot([b.ts(light_onsets(1))-0.1 b.ts(onsets(10))+0.2],[ylims(2) ylims(2)],'r','linewidth',1.5)
    plot([b.ts(light_onsets(21))-0.1 b.ts(onsets(50))+0.2],[ylims(2) ylims(2)],'r','linewidth',1.5)
    upfac = 5;
    text(b.ts(light_onsets(1)),ylims(2)+upfac,'1-Light','FontSize',6);
    text(b.ts(light_onsets(21)),ylims(2)+upfac,'6-Light','FontSize',6);
    text(b.ts(onsets(1)),ylims(2)+upfac,'2-Air','FontSize',6);
    text(b.ts(onsets(11)),ylims(2)+upfac,'3-Air','FontSize',6);
    text(b.ts(onsets(21)),ylims(2)+upfac,'4-Air-Light','FontSize',6);
    text(b.ts(onsets(31)),ylims(2)+upfac,'5-Air','FontSize',6);
    text(b.ts(onsets(41)),ylims(2)+upfac,'7-Air','FontSize',6);
    text(b.ts(light_onsets(1))-1.45,ylims(2)+upfac+5,'Configurations','FontSize',6);
    set(gca,'Ydir','Normal');
    xlabel('Time (min)');
    ylabel('Cell Number (arranged by peaks)');
    set(gca,'FontSize',6,'FontWeight','Normal');
    changePosition(gca,[-0.075 -0.005 0.16 -0.07]);
    ht = title(sprintf('Normalized Calcium signal (%cF/Fo) of %d of %d cells from a representative animal',916,numcells,size(signals,1))); 
    changePosition(ht,[0 12 0]);
    set(ht,'FontSize',6,'FontWeight','Normal');
    legs = {'Light onset','Air onset','Air offset','Speed',[0.85 0.15 111 0.1]};
    putLegendH(gca,legs,'colors',colors);
    hc = putColorBar(gca,[0.75 0.62 -0.8 -0.6],{'0',sprintf('>%.1f',maxSig)},6,'northoutside',[0.15 0.5 0.045 0.5]);
    colormap_ig
    
%     colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('overall_pop.pdf'),600);
    break;
end
return;
%%
while 0
    hf = figure(100);clf;set(gcf,'Units','Inches');set(gcf,'Position',[1 5 5.69 2.45],'color','w'); hold on;
    spSigAllN = normalizeSignal(signals,2);
    [maxVal,maxLoc] = max(spSigAllN,[],2);
    [sorted,locs] = sort(maxLoc);
    maxSig = 0.5;
    hi = imagesc([traceTime(1), traceTime(end)],[1 size(spSigAllN,1)],spSigAllN(locs,:),[0 maxSig]);
    % colorbar;
    hold on;
    % plot(b.ts(frames_f),spSigAll(:,:));
    ylims = ylim;
    hold on;
    lwdth = 0.1;
    colors = {[1 0 1]/1.5;[0 0.3 1];[1 1 0]/1.75;[0 1 1]/1.25};
    if flag(1)
    if isfield(b,'air_puff_raw')
    %    plot(b.ts,(b.air_puff_raw-0.5)*ylims(2)*4,'color',mData.colors{3},'linewidth',0.5);
       for ii = 1:length(onsets)
           tt = b.ts(onsets(ii));
           plot([tt tt],[ylims(1)-5 ylims(2)+5],'color',colors{2},'linewidth',lwdth);
           tt = b.ts(offsets(ii));
           plot([tt tt],[ylims(1)-5 ylims(2)+5],'color',colors{3},'linewidth',lwdth);
       end
    end

    end
    if flag(2)
    if isfield(b,'stim_raw')
    %    plot(b.ts,(b.stim_raw-0.5)*ylims(2)*4,'m','linewidth',0.5);
       for ii = 1:length(light_onsets)
           tt = b.ts(light_onsets(ii));
           plot([tt tt],[ylims(1)-5 ylims(2)+5],'color',colors{1},'linewidth',lwdth);
       end
    end
    end
    if flag(3)
%     plot(b.ts,b.photo_sensor_raw*ylims(2),mData.colors{3},'linewidth',0.5)
    end
    if flag(4)
    plot(b.ts,b.fSpeed*0.25,'color',colors{4},'linewidth',lwdth);
    end
    numcells = 100;
%     numcells = size(signals,1) % added for Bruce
    ylim([0 numcells]); xlim([0 traceTime(end)]);
%     ylims = ylim;
    upfac = 3;
    upfac = -1; % added for Bruce
    yte = numcells+5; % used to be ylims(2)+upfac
    text(b.ts(light_onsets(1)),yte,'1-Light','FontSize',6);
    text(b.ts(light_onsets(21)),yte,'6-Light','FontSize',6);
    text(b.ts(onsets(1)),yte,'2-Air','FontSize',6);
    text(b.ts(onsets(11)),yte,'3-Air','FontSize',6);
    text(b.ts(onsets(21)),yte,'4-Air-Light','FontSize',6);
    text(b.ts(onsets(31)),yte,'5-Air','FontSize',6);
    text(b.ts(onsets(41)),yte,'7-Air','FontSize',6);
    text(b.ts(light_onsets(1))-1.45,yte+5,'Configurations','FontSize',6);
    set(gca,'Ydir','Normal');
    xlabel('Time (min)'); 
    ylabel('Cell Number (arranged by peak firing time)');
    set(gca,'FontSize',6,'FontWeight','Normal');
    changePosition(gca,[-0.075 -0.005 0.16 -0.07]);
    ht = title(sprintf('Normalized Calcium signal (%cF/Fo) of %d of %d cells from a representative animal',916,numcells,size(signals,1))); 
    changePosition(ht,[0 12 0]);
    set(ht,'FontSize',6,'FontWeight','Normal');
%     legs = {'Light onset','Air onset','Air offset','Speed',[0.5 0.1 111 0.1]};
    legs = {'Light onset','Air onset','Air offset','Speed',[2 0.1 numcells+11 0.1]};
    putLegendH(gca,legs,'colors',colors);
    hc = putColorBar(gca,[0.75 0.62 -0.8 -0.6],{'0',sprintf('>%.1f',maxSig)},6,'northoutside',[0.15 0.4 0.045 0.4]);
    colormap_ig
%     colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('overall_pop.pdf'),600);
    break;
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
