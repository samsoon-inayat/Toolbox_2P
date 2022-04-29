function trace_plot_all_dist(pd_rec,pl)

if ~exist('flag','var')
    flag = [1 1 0 1];
end

if ~exist('pd_rec','var')
    pd_rec = evalin('base','ei{4}');
    pl = 1;
end
try
    signals = pd_rec.deconv';
    frames_f = pd_rec.b.frames_f;
catch
    signals = pd_rec.plane{pl}.tP.deconv.spSigAll;
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
tei = pd_rec; pp = pl; markersOn = frames_f(1); markersOff = frames_f(end); thisRasterType = 'dist'; binwidths = [0.3 3]; owr = 0;
rasters = make_rasters_whole(tei,pp,markersOn,markersOff,thisRasterType,binwidths,signals,owr);
% [nbins,binWidth] = get_nbins_F(b,markersOn,markersOff,binwidths(2),thisRasterType);
onsets = update_markers_dist(b,onsets(11:40),rasters);
offsets = update_markers_dist(b,offsets(11:40),rasters);
light_onsets = update_markers_dist(b,light_onsets(11:20),rasters);
signals = squeeze(rasters.sp_rasters1);
signals = signals';
b.ts_t = b.ts;
b.ts = rasters.distL1/100;
traceTime = b.ts;
b.fSpeed_old = b.fSpeed;
b.fSpeed = rasters.speed(1:size(signals,2));
n = 0;
%%
while 1
    hf = figure(100);clf;set(gcf,'Units','Inches');set(gcf,'Position',[1 5 6.9 5],'color','w'); hold on;
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
    numcells = size(signals,1); % added for Bruce
    ylim([0 numcells]); xlim([0 traceTime(end)]);
%     ylims = ylim;
    upfac = 3;
    upfac = -1; % added for Bruce
    yte = numcells+5; % used to be ylims(2)+upfac
%     text(b.ts(light_onsets(1)),yte,'1-Light','FontSize',6);
%     text(b.ts(light_onsets(21)),yte,'6-Light','FontSize',6);
%     text(b.ts(onsets(1)),yte,'2-Air','FontSize',6);
    text(b.ts(onsets(11-10)),yte,'3-Air','FontSize',6);
    text(b.ts(onsets(21-10)),yte,'4-Air-Light','FontSize',6);
    text(b.ts(onsets(31-10)),yte,'5-Air','FontSize',6);
%     text(b.ts(onsets(41)),yte,'7-Air','FontSize',6);
    text(b.ts(onsets(11-10))-1.45,yte+5,'Conditions','FontSize',6);
    set(gca,'Ydir','Normal');
    xlabel('Distance (m)');
    ylabel('Cell Number (arranged by peaks)');
    set(gca,'FontSize',6,'FontWeight','Normal');
    changePosition(gca,[-0.075 -0.005 0.16 -0.07]);
    ht = title(sprintf('Normalized Calcium signal (%cF/Fo) of %d of %d cells from a representative animal',916,numcells,size(signals,1))); 
    changePosition(ht,[0 12 0]);
    set(ht,'FontSize',6,'FontWeight','Normal');
%     legs = {'Light onset','Air onset','Air offset','Speed',[0.5 0.1 111 0.1]};
    legs = {'Light onset','Air onset','Air offset','Speed',[0.5 0.1 numcells+20 0.1]};
    putLegendH(gca,legs,'colors',colors);
    hc = putColorBar(gca,[0.75 0.62 -0.8 -0.6],{'0',sprintf('>%.1f',maxSig)},6,'northoutside',[0.15 0.4 0.045 0.4]);
    colormap_ig
%     colormap jet
    save_pdf(hf,mData.pdf_folder,sprintf('overall_pop_dist_%s_%d.pdf',pd_rec.recordingFolder(61:66),pl),600);
    disp(pd_rec.recordingFolder(61:66));
    break;
end

function onsetso = update_markers_dist(b,onsets,rasters)
markers = rasters.bin_markers;
st = markers(:,1); se = markers(:,2);
for ii = 1:length(onsets)
    to = onsets(ii);
    indst = find((st-to)>0,1,'first');
    if length(indst) > 1 || isempty(indst)
        onsetso(ii) = onsetso(ii-1);
        continue;
    end
    onsetso(ii) = indst;
end


