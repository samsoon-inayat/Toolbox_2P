function allA_P = getRasters_motion_offsets(ei,trials,ow)

if ~exist('ow','var')
    ow = 0;
end

if iscell(ei)
    for ii = 1:length(ei)
        temp = getRasters_motion_offsets(ei{ii},trials,ow);
        if ii == 1
            cellRaster = temp.cells;
            rasters = temp.rasters;
        else
            cellRaster = [cellRaster temp.cells];
            rasters = cat(3,rasters,temp.rasters);
        end
    end
    allA_P = temp;
    allA_P.cells = cellRaster;
    allA_P.rasters = rasters;
    return;
end

if ~exist('ei','var')
    ei = evalin('base','ei{1}');    
end

onsets = ei.b.air_puff_r(trials);
offsets = ei.b.air_puff_f(trials);
% temp = ei.b.fSpeed < 0.5;
% temp1 = ei.b.fSpeed > 0.5;
motionOnsets = find_rising_edge(ei.b.fSpeed > 0.5,0.1,5000);
motionOffsets = find_falling_edge(ei.b.fSpeed > 0.5,-0.1,5000);
temp = find(diff(motionOnsets)<=10000);
motionOnsets(temp+1) = [];
temp = find(diff(motionOffsets)<=10000);
motionOffsets(temp+1) = [];

motionOnsets = motionOnsets(motionOnsets>onsets(1));
motionOffsets = motionOffsets(motionOffsets>onsets(1));

for ii = 1:length(onsets)
    inds = find(motionOnsets>onsets(ii) & motionOnsets<offsets(ii));
    motionOnsets(inds) = [];
end
temp = motionOffsets;
motionOffsets = [];
for ii = 1:length(motionOnsets)
    inds = find(temp > motionOnsets(ii),1,'first');
    motionOffsets(ii) = temp(inds);
end

diffs = motionOffsets - motionOnsets;
inds = find(diffs < 10000);
motionOnsets(inds) = [];
motionOffsets(inds) = [];

bothM = sort([motionOnsets motionOffsets]);
temp = find(diff(bothM)<=10000);
bothM(temp+1) = [];
motionOnsets = intersect(bothM,motionOnsets);
motionOffsets = intersect(bothM,motionOffsets);

temp = motionOffsets;
motionOffsets = [];
for ii = 1:length(motionOnsets)
    inds = find(temp > motionOnsets(ii),1,'first');
    motionOffsets(ii) = temp(inds);
end

%check
% figure(101);clf
% plot(ei.b.ts,ei.b.fSpeed);hold on;
% tMOn = zeros(size(ei.b.fSpeed));
% tMOn(motionOnsets) = 15;
% plot(ei.b.ts,tMOn,'m');
% tMOn = zeros(size(ei.b.fSpeed));
% tMOn(motionOffsets) = 15;
% plot(ei.b.ts,tMOn,'c');


timeBefore = 2; % sec
timeAfter = 2; % sec

a_p_onsets = motionOffsets - round(1e6 * timeBefore/ei.b.si);
a_p_offsets = motionOffsets + round(1e6 * timeAfter/ei.b.si);
temp = find(a_p_onsets < 0);
a_p_onsets(temp) = [];
a_p_offsets(temp) = [];
% plotMarkers(ei,a_p_onsets,a_p_offsets,101);
allA_P = getTimeRaster(ei,a_p_onsets,a_p_offsets,ow,0);




function edge = find_rising_edge(signal,threshold,minimum_diff)
edge = find(diff(signal) >= threshold);
temp = find(diff(edge)<=minimum_diff);
edge(temp+1) = [];

function edge = find_falling_edge(signal,threshold,minimum_diff)
edge = find(diff(signal) <= threshold);
temp = find(diff(edge)<=minimum_diff);
edge(temp+1) = [];