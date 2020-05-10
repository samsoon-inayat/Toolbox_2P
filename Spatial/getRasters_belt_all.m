function allA_P = getRasters_belt_all(ei,trials,ow)

if ~exist('ow','var')
    ow = 0;
end

if iscell(ei)
    for ii = 1:length(ei)
        temp = getRasters_belt_all(ei{ii},trials,ow);
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

onsets = ei.b.air_puff_r;
offsets = ei.b.air_puff_f;
photo_sensor = ei.b.photo_sensor_f(ei.b.photo_sensor_f>onsets(1) & ei.b.photo_sensor_f<offsets(end));
dists = ei.b.dist(photo_sensor);
diff_dists = diff(dists);
inds = find(diff_dists < 100);
temp_photo_sensor = photo_sensor;
temp_photo_sensor(inds+1) = [];


a_p_onsets = temp_photo_sensor(1:(length(temp_photo_sensor)-1));
a_p_offsets = temp_photo_sensor(2:end);
% plotMarkers(ei,a_p_onsets,a_p_offsets,101);
allA_P = getDistRaster(ei,a_p_onsets,a_p_offsets,ow,0);


