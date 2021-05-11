function rasters =  getRasters_fixed_bin_width_ctrl(ei,pp,onsets,offsets,rasterType)

if strcmp(rasterType,'time')
    binWidth = 0.15; % unit of bin width is sec
    bins = 0:binWidth:50;
    maxbins = 500;
end
if strcmp(rasterType,'dist')
    binWidth = 1.75; % unit of bin width is cm
    bins = 0:binWidth:1000;
    maxbins = 95;
end

ccs = 1:length(ei.tP.deconv.spSigAll);
% ccs = find(ei.tP.iscell(:,1));%    1:length(ei.tP.deconv.spSigAll);
b = ei.b;
b.frames_f = ei.plane{pp}.b.frames_f;
b.frameRate = ei.thorExp.frameRate;
% spSigAll = zeros(length(ccs),length(ei.deconv.spSigAll{1}));
% for ii = 1:length(ccs)
%     cellNum = ccs(ii);
%     spSigAll(ii,:) = ei.deconv.spSigAll{cellNum}';
% end
spSigAll = ei.deconv.spSigAll;
rasters = getRasters(b,spSigAll,onsets,offsets,binWidth,maxbins,rasterType);
rasters.onsets = onsets; 
rasters.offsets = offsets;
rasters = fixRastersForNaN(rasters);
rasters.xs = bins(1:(min(rasters.lastBin)-1));
if strcmp(rasterType,'time')
    rasters.cell_history = getCellHistory(ei,b,onsets,rasters,binWidth);
    if isfield(b,'air_puff_raw')
        out =  getTimeRaster_2(b,spSigAll,onsets,offsets,binWidth);
        rasters.wholeData = out;
    end
    out =  getRasters_from_frames(b,spSigAll,onsets,offsets);
    rasters.fromFrames = out;
end


function out =  getRasters(b,spSignal,onsets,offsets,bin_width,nbins,rasterType)
raster = NaN(length(onsets),nbins,size(spSignal,1));
speed = NaN(length(onsets),nbins);
dist = NaN(length(onsets),nbins);
space = NaN(length(onsets),nbins);
duration = NaN(length(onsets),nbins);
timeStart = duration;
timeEnd = duration;
mean_belt_length = mean(diff(b.dist(b.photo_sensor_f)));
lastBin = NaN(length(onsets),1);
for trial = 1:length(onsets)
    st = onsets(trial);
    se = offsets(trial);
    if strcmp(rasterType,'time')
        [st1,se1] = getTimeBins(b,bin_width,nbins,st);
    end
    if strcmp(rasterType,'dist')
        [st1,se1] = getDistBins(b,bin_width,nbins,st,se);
    end
    for bin = 1:nbins
        if se1(bin) > se || isnan(st1(bin))
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
    lastBin(trial) = bin;
end
out.sp_rasters = raster; out.speed = speed; out.dist = dist; out.space = space; out.duration = duration;
out.timeStart = timeStart; out.timeEnd = timeEnd; out.lastBin = lastBin;

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
%             spSig = ei.deconv.spSigAll{cn}';
            spSig = ei.deconv.spSigAll(cn,:);
            if ~isempty(frames)
                cell_history(trial,bin,cn) = mean(spSig(frames));
            else
                cell_history(trial,bin,cn) = NaN;
            end
        end
    end
end
cell_history = fillmissing(cell_history,'linear',2,'endValues','nearest');

function [st1,se1] = getDistBins(b,bin_width,nbins,st,se)
st1 = NaN(nbins,1);
se1 = st1;
for jj = 1:nbins
    if jj == 1
        st1(jj) = st;
    else
        st1(jj) = se1(jj-1);
    end
    temp = b.dist(st1(jj)) + bin_width;
    temp = (b.dist(st1(jj):se) - temp);
    if temp(end) < 0
        st1(jj) = NaN;
        break;
    end
    temp = temp.^2;
    [~,ind] = min(temp);
    se1(jj) = ind + st1(jj);
end

function [st1,se1] = getTimeBins(b,bin_width,nbins,st)
n_samples = floor(bin_width/(b.si*1e-6));
total_samples = nbins * n_samples;
se = st + total_samples;
st1 = st:n_samples:(se-n_samples);
se1 = ((st+n_samples):n_samples:se)-1;

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

function out =  getRasters_from_frames(b,spSignal,onsets,offsets)
for trial = 1:length(onsets)
    st = onsets(trial);
    se = offsets(trial);
    frames = b.frames_f > st & b.frames_f < se;
    nFrames(trial) = sum(frames);
end
nbins = max(nFrames);
raster = NaN(length(onsets),nbins,size(spSignal,1));
cell_history = NaN(length(onsets),6,size(spSignal,1));
speed = NaN(length(onsets),nbins);
dist = NaN(length(onsets),nbins);
space = NaN(length(onsets),nbins);
duration = NaN(length(onsets),nbins);
timeStart = duration;
timeEnd = duration;
mean_belt_length = mean(diff(b.dist(b.photo_sensor_f)));
for trial = 1:length(onsets)
    st = onsets(trial);
    se = offsets(trial);
    frames = b.frames_f > st & b.frames_f < se;
    nF = nFrames(trial); 
    temp = spSignal(:,frames);
    raster(trial,1:nF,:) = temp';
    speed(trial,1:nF) = b.fSpeed(b.frames_f(frames));
    dist(trial,1:nF) = b.dist(b.frames_f(frames))-b.dist(st);
    frames_list = find(frames);
    for bin = 1:nF
        bsl = b.frames_f(frames_list(bin));
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
    end
    duration(trial,1:nF) = b.ts(b.frames_f(frames))-b.ts(st);
    timeStart(trial,1:nF) = b.ts(b.frames_f(frames));
    timeEnd(trial,1:nF) = b.ts(b.frames_f(frames));
    hist_frames = b.frames_f < st;
    temp = spSignal(:,hist_frames);
    cell_history(trial,:,:) = temp(:,(size(temp,2)-5):end)';
end
out.sp_rasters = raster; out.speed = speed; out.dist = dist; out.space = space; out.duration = duration;
out.timeStart = timeStart; out.timeEnd = timeEnd; out.cell_history = cell_history;