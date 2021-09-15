function trace_plot_all(pd_rec,flag)
if ~exist('flag','var')
    flag = [1 1 0 1];
end
if ~exist('pd_rec','var')
    pd_rec = evalin('base','ei{4}');
end
try
signals = pd_rec.deconv';
frames_f = pd_rec.b.frames_f;
catch
    signals = pd_rec.plane{1}.tP.deconv.spSigAll;
    frames_f = pd_rec.plane{1}.b.frames_f;
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
signals = get_calcium_data(pd_rec,1);
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
    text(b.ts(light_onsets(1)),ylims(2)+upfac,'C1','FontSize',6);
    text(b.ts(light_onsets(21)),ylims(2)+upfac,'C1''','FontSize',6);
    text(b.ts(onsets(1)),ylims(2)+upfac,'C2','FontSize',6);
    text(b.ts(onsets(11)),ylims(2)+upfac,'C3','FontSize',6);
    text(b.ts(onsets(21)),ylims(2)+upfac,'C4','FontSize',6);
    text(b.ts(onsets(31)),ylims(2)+upfac,'C3''','FontSize',6);
    text(b.ts(onsets(41)),ylims(2)+upfac,'C2''','FontSize',6);
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
    colormap_ig
    save_pdf(hf,mData.pdf_folder,sprintf('overall_pop.pdf'),600);
    break;
end

%%
%%
while 1
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
    text(b.ts(light_onsets(1)),ylims(2)+upfac,'C1','FontSize',6);
    text(b.ts(light_onsets(21)),ylims(2)+upfac,'C1''','FontSize',6);
    text(b.ts(onsets(1)),ylims(2)+upfac,'C2','FontSize',6);
    text(b.ts(onsets(11)),ylims(2)+upfac,'C3','FontSize',6);
    text(b.ts(onsets(21)),ylims(2)+upfac,'C4','FontSize',6);
    text(b.ts(onsets(31)),ylims(2)+upfac,'C3''','FontSize',6);
    text(b.ts(onsets(41)),ylims(2)+upfac,'C2''','FontSize',6);
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
    colormap_ig
    save_pdf(hf,mData.pdf_folder,sprintf('overall_pop.pdf'),600);
    break;
end