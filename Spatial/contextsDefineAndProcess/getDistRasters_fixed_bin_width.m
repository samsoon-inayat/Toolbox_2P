function rasters =  getDistRasters_fixed_bin_width(ei,pp,onsets,offsets)

binWidth = 3; % unit of bin width is cm
bins = 0:binWidth:1000;
maxbins = 55;
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
rasters = getDistRaster_2(b,spSigAll,onsets,offsets,binWidth,maxbins);
rasters.xs = bins;
rasters.onsets = onsets; 
rasters.offsets = offsets;

function out =  getDistRaster_2(b,spSignal,onsets,offsets,bin_width,maxbins)
raster = NaN(length(onsets),maxbins,size(spSignal,1));
speed = NaN(length(onsets),maxbins);
dist = NaN(length(onsets),maxbins);
space = NaN(length(onsets),maxbins);
mean_belt_length = mean(diff(b.dist(b.photo_sensor_f)));
duration = NaN(length(onsets),maxbins);
timeStart = duration;
timeEnd = duration;
% Dist = removeDistOutliers(Dist);
for trial = 1:length(onsets)
    st = onsets(trial);
    se = offsets(trial);
    st1 = NaN(maxbins,1);
    se1 = st1;
    for jj = 1:maxbins
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
    for bin = 1:maxbins
        if isnan(st1(bin))
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

