function allA_P = getRasters_air_onsets_offsets(ei,trials,ow)

if ~exist('ow','var')
    ow = 0;
end

if iscell(ei)
    for ii = 1:length(ei)
        temp = getRasters_air_onsets_offsets(ei{ii},trials,ow);
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

timeBefore = 2; % sec
timeAfter = 2; % sec

a_p_onsets = onsets - round(1e6 * timeBefore/ei.b.si);
a_p_offsets = onsets + round(1e6 * timeAfter/ei.b.si);
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