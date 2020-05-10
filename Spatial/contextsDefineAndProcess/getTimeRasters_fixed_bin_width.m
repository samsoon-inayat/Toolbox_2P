function rasters =  getTimeRasters_fixed_bin_width(ei,pp,onsets,offsets)

binWidth = 0.2; % unit of bin width is sec
bins = 0:binWidth:50;
maxbins = 75;
ccs = 1:length(ei.tP.deconv.spSigAll);
b = ei.b;
b.frames_f = ei.plane{pp}.b.frames_f;
b.frameRate = ei.thorExp.frameRate;
spSigAll = zeros(length(ccs),length(ei.deconv.spSigAll{1}));
% caSigAll = zeros(length(ccs),length(ei.deconv.caSigAll{1}));
for ii = 1:length(ccs)
    spSigAll(ii,:) = ei.deconv.spSigAll{ii}';
%     caSigAll(ii,:) = ei.deconv.caSigAll{ii}';
end
rasters = getTimeRaster_1(b,spSigAll,onsets,offsets,binWidth,maxbins);
rasters.xs = bins;
rasters.onsets = onsets; 
rasters.offsets = offsets;
rasters.cell_history = getCellHistory(ei,b,onsets,rasters,binWidth);
out =  getTimeRaster_2(b,spSigAll,onsets,offsets,binWidth);
rasters.wholeData = out;

function out =  getTimeRaster_1(b,spSignal,onsets,offsets,bin_width,maxbins)
n_samples = floor(bin_width/(b.si*1e-6));
nbins = maxbins;
total_samples = nbins * n_samples;
raster = NaN(length(onsets),maxbins,size(spSignal,1));
speed = NaN(length(onsets),maxbins);
dist = NaN(length(onsets),maxbins);
space = NaN(length(onsets),maxbins);
mean_belt_length = mean(diff(b.dist(b.photo_sensor_f)));
duration = NaN(length(onsets),maxbins);
timeStart = duration;
timeEnd = duration;
for trial = 1:length(onsets)
    st = onsets(trial);
    se = st + total_samples;
    st1 = st:n_samples:(se-n_samples);
    se1 = ((st+n_samples):n_samples:se)-1;
    for bin = 1:nbins
        if se1(bin) > offsets(trial)
            break;
        end
        frames = b.frames_f > st1(bin) & b.frames_f < se1(bin);
        temp = spSignal(:,frames); temp = mean(temp,2);
        raster(trial,bin,:) = temp;

        speed(trial,bin) = mean(b.fSpeed(st1(bin):se1(bin)));
        dist(trial,bin) = mean(b.dist(st1(bin):se1(bin)))-b.dist(st);
        bsl = st1(bin) + floor((se1(bin)-st1(bin))/2);
        inds_p = find(b.photo_sensor_f < bsl,1,'last');
        photo_sensor_place = 'behind';
        if isempty(inds_p)
            inds_p = find(b.photo_sensor_f > bsl,1,'first');
            photo_sensor_place = 'ahead';
        end
        photo_sensor = b.photo_sensor_f(inds_p);
        if strcmp(photo_sensor_place,'ahead')
            space(trial,bin) = mean_belt_length - (b.dist(photo_sensor) - b.dist(bsl));
        else
            space(trial,bin) = -(b.dist(photo_sensor) - b.dist(bsl));
        end
        duration(trial,bin) = b.ts(se1(bin))-b.ts(st1(bin));
        timeStart(trial,bin) = b.ts(st1(bin));
        timeEnd(trial,bin) = b.ts(se1(bin));
    end
end
out.sp_rasters = raster; out.speed = speed; out.dist = dist; out.space = space; out.duration = duration;
out.timeStart = timeStart; out.timeEnd = timeEnd;

function out =  getTimeRaster_2(b,spSignal,onsets,offsets,bin_width)
n_samples = floor(bin_width/(b.si*1e-6));
S1 = onsets(1) - (10*n_samples);
S2 = offsets(end) + (10*n_samples);
nbins = floor((S2-S1)/n_samples);
total_samples = nbins * n_samples;
S2 = S1+total_samples;
st = S1;
se = S2;
st1 = st:n_samples:(se-n_samples);
se1 = ((st+n_samples):n_samples:se)-1;
raster = NaN(nbins,size(spSignal,1));
speed = NaN(nbins,1);
dist = NaN(nbins,1);
space = NaN(nbins,1);
air = space;
mean_belt_length = mean(diff(b.dist(b.photo_sensor_f)));
duration = NaN(nbins,1);
timeStart = duration;
timeEnd = duration;
for bin = 1:nbins
    frames = (b.frames_f > st1(bin) & b.frames_f < se1(bin));
    temp = spSignal(:,frames); raster(bin,:) = mean(temp,2)';
    speed(bin) = mean(b.fSpeed(st1(bin):se1(bin)));
    dist(bin) = mean(b.dist(st1(bin):se1(bin)))-b.dist(st);
    air(bin) = mean(b.air_puff_raw(st1(bin):se1(bin))) > 0.5;
    bsl = st1(bin) + floor((se1(bin)-st1(bin))/2);
    inds_p = find(b.photo_sensor_f < bsl,1,'last');
    photo_sensor_place = 'behind';
    if isempty(inds_p)
        inds_p = find(b.photo_sensor_f > bsl,1,'first');
        photo_sensor_place = 'ahead';
    end
    photo_sensor = b.photo_sensor_f(inds_p);
    if strcmp(photo_sensor_place,'ahead')
        space(bin) = mean_belt_length - (b.dist(photo_sensor) - b.dist(bsl));
    else
        space(bin) = -(b.dist(photo_sensor) - b.dist(bsl));
    end
    duration(bin) = b.ts(se1(bin))-b.ts(st1(bin));
    timeStart(bin) = b.ts(st1(bin))-b.ts(st1(1));
    timeEnd(bin) = b.ts(se1(bin))-b.ts(st1(1));
    timeMid(bin) = timeStart(bin) + ((timeEnd(bin)-timeStart(bin))/2);
end
out.sp_rasters = raster; out.speed = speed; out.dist = dist; out.space = space; out.air = air; out.duration = duration;
out.timeStart = timeStart; out.timeEnd = timeEnd; out.xs = timeMid;

n_samples = floor(bin_width/(b.si*1e-6));
nbins = 6;
total_samples = 6 * n_samples;
se = st1(1);
st = se - total_samples;
st1 = st:n_samples:(se-n_samples);
se1 = ((st+n_samples):n_samples:se)-1;
cell_history = NaN(size(spSignal,1),nbins);
for bin = 1:nbins
    frames = find(b.frames_f > st1(bin) & b.frames_f < se1(bin));
    temp = spSignal(:,frames); 
    cell_history(:,bin) = mean(temp,2);
end
out.cell_history = fillmissing(cell_history,'linear',2,'endValues','nearest');


function cell_history = getCellHistory(ei,b,onsets,rasters,bin_width)
n_samples = floor(bin_width/195e-6);
nbins = 6;
total_samples = 6 * n_samples;
cell_history = NaN(length(onsets),nbins,size(rasters.sp_rasters,3));
for trial = 1:length(onsets)
    se = onsets(trial);
    st = se - total_samples;
    st1 = st:n_samples:(se-n_samples);
    se1 = ((st+n_samples):n_samples:se)-1;
    for bin = 1:nbins
        frames = find(b.frames_f > st1(bin) & b.frames_f < se1(bin));
        for cn = 1:size(rasters.sp_rasters,3)
            spSig = ei.deconv.spSigAll{cn}';
            if ~isempty(frames)
                cell_history(trial,bin,cn) = mean(spSig(frames));
            else
                cell_history(trial,bin,cn) = NaN;
            end
        end
    end
end
cell_history = fillmissing(cell_history,'linear',2,'endValues','nearest');

% function out =  getSpeedRaster(b,onsets,offsets,bins,maxbins)
% sbins = bins(1:(end-1));
% ebins = bins(2:end);
% trialNumbers = 1:length(onsets);
% tss = cell(length(trialNumbers),1);
% dist = cell(length(trialNumbers),1);
% speed = cell(length(trialNumbers),1);
% for iii = 1:length(trialNumbers)
%     ii = trialNumbers(iii);
%     st = onsets(ii);
%     se = offsets(ii);
%     temp = b.fSpeed(st:se);
%     temp(temp<0) = 0;
%     speed{iii} = temp;
%     tss{iii} = b.ts(st:se)-b.ts(st);
%     dist{iii} = b.dist(st:se)-b.dist(st);
% end
% speed_raster = NaN(length(trialNumbers),maxbins);
% dist_raster = NaN(length(trialNumbers),maxbins);
% space_raster = NaN(length(trialNumbers),maxbins);
% duration_raster = NaN(length(trialNumbers),maxbins);
% mean_belt_length = mean(diff(b.dist(b.photo_sensor_f)));
% for iii = 1:length(trialNumbers)
%     thisT = tss{iii};
%     thisSpeed = speed{iii};
%     thisDist = dist{iii};
%     st_trial = onsets(trialNumbers(iii));
%     for jj = 1:maxbins
%         st = sbins(jj);  se = ebins(jj);
%         dVals = (thisT >= st & thisT<se);
%         if sum(dVals) == 0
%             break;
%         end
%         tVals = thisT(dVals);
%         stt = tVals(1);
%         sett = tVals(end);
%         speed_raster(iii,jj) = mean(thisSpeed(dVals));
%         dist_raster(iii,jj) = max(thisDist(dVals));
%         duration_raster(iii,jj) = sett-stt;
%         bin_sample_location = st_trial + find(dVals);
%         bsl1 = bin_sample_location(1); bsl2 = bin_sample_location(end);
%         bsl = bsl1 + floor((bsl2-bsl1)/2);
%         inds_p = find(b.photo_sensor_f < bsl,1,'last');
%         photo_sensor_place = 'behind';
%         if isempty(inds_p)
%             inds_p = find(b.photo_sensor_f > bsl,1,'first');
%             photo_sensor_place = 'ahead';
%         end
%         photo_sensor = b.photo_sensor_f(inds_p);
%         if strcmp(photo_sensor_place,'ahead')
%             space_raster(iii,jj) = mean_belt_length - (b.dist(photo_sensor) - b.dist(bsl));
%         else
%             space_raster(iii,jj) = -(b.dist(photo_sensor) - b.dist(bsl));
%         end
% %         plotPoints = 1:10:maxbins;
% %         if sum(plotPoints == jj) > 0
% %             figure(1000);clf;plot(b.photo_sensor_raw);hold on;
% %             plot(b.air_puff_raw);
% %             plot(bsl,0.5,'r*');
% %             plot(b.fSpeed/max(b.fSpeed),'m');
% %             xlim([0 b.air_puff_f(3)]);
% %             title(space_raster(iii,jj));
% %             pause(0.001);
% %         end
%     end
% end
% out.speed_raster = speed_raster;
% out.dist_raster = dist_raster;
% out.duration_raster = duration_raster;
% out.space_raster = space_raster;
% out.xs = bins(1:maxbins);


