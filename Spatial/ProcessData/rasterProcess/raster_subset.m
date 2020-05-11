function rasters = raster_subset(rasters,lastcol)

rasters.sp_rasters_subset = rasters.sp_rasters;%(:,1:lastcol,:);
rasters.duration_subset = rasters.duration;%(:,1:lastcol);
return;

Rs = rasters.sp_rasters;
Dur = rasters.duration;

for ii = 1:size(Rs,3)
    thisR = Rs(:,:,ii);
    [rnan,cnan] = find(isnan(thisR));
    for jj = 1:length(rnan)
        if cnan(jj) == 1
            continue;
        end
        if cnan(jj) == size(thisR,2)
            continue;
        end
        valB = thisR(rnan(jj),cnan(jj)-1);
        valA = thisR(rnan(jj),cnan(jj)+1);
        if ~isnan(valB) && ~isnan(valA)
            thisR(rnan(jj),cnan(jj)) = mean([valB valA]);
        end
    end
    for jj = 1:length(rnan)
    [rnan,cnan] = find(isnan(thisR));
    dcnan = diff(cnan);
    end
end




