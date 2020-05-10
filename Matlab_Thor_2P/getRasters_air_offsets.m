function allA_P = getRasters_air_offsets(ei,trials,ow)

if ~exist('ow','var')
    ow = 0;
end

if iscell(ei)
    for ii = 1:length(ei)
        temp = getRasters_air_onsets(ei{ii},trials,ow);
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

% temp = ei.b.fSpeed < 0.5;
% temp1 = ei.b.fSpeed > 0.5;
motionOnsets = find_rising_edge(ei.b.fSpeed > 0.5,0.1,5000);
motionOffsets = find_falling_edge(ei.b.fSpeed > 0.5,-0.1,5000);
bothM = sort([motionOnsets motionOffsets]);
temp = find(diff(bothM)<=20000);
bothM(temp+1) = [];
motionOnsets = intersect(bothM,motionOnsets);
motionOffsets = intersect(bothM,motionOffsets);
% %check
% figure(101);clf
% plot(ei.b.ts,ei.b.fSpeed);hold on;
% tMOn = zeros(size(ei.b.fSpeed));
% tMOn(motionOnsets) = 15;
% plot(ei.b.ts,tMOn,'m');
% tMOn = zeros(size(ei.b.fSpeed));
% tMOn(motionOffsets) = 15;
% plot(ei.b.ts,tMOn,'c');
onsets = ei.b.air_puff_r(trials);
offsets = ei.b.air_puff_f(trials);
% motionOnsets = motionOnsets(motionOnsets>onsets(1) & motionOnsets<offsets(end));
% motionOffsets = motionOffsets(motionOffsets>onsets(1) & motionOffsets<offsets(end));

timeBefore = 2; % sec
timeAfter = 2; % sec

a_p_onsets = offsets - round(1e6 * timeBefore/ei.b.si);
a_p_offsets = offsets + round(1e6 * timeAfter/ei.b.si);
temp = find(a_p_onsets < 0);
a_p_onsets(temp) = [];
a_p_offsets(temp) = [];
% plotMarkers(ei,a_p_onsets,a_p_offsets,101);
allA_P = getTimeRaster(ei,a_p_onsets,a_p_offsets,ow,0);

% a_p_onsets = motionOffsets - round(1e6 * timeBefore/ei.b.si);
% a_p_offsets = motionOffsets + round(1e6 * timeAfter/ei.b.si);
% temp = find(a_p_onsets < 0);
% a_p_onsets(temp) = [];
% a_p_offsets(temp) = [];
% % plotMarkers(ei,a_p_onsets,a_p_offsets,101);
% allA_P.motionOffset = getTimeRaster(ei,a_p_onsets,a_p_offsets,ow,0);



function edge = find_rising_edge(signal,threshold,minimum_diff)
edge = find(diff(signal) >= threshold);
temp = find(diff(edge)<=minimum_diff);
edge(temp+1) = [];

function edge = find_falling_edge(signal,threshold,minimum_diff)
edge = find(diff(signal) <= threshold);
temp = find(diff(edge)<=minimum_diff);
edge(temp+1) = [];