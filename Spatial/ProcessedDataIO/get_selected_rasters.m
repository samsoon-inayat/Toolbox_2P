function rasters = get_selected_rasters(rasters,all_ccs)

for rr = 1:size(rasters,1)
    for cc = 1:size(rasters,2)
        thisMR = rasters{rr,cc};
        rasters{rr,cc} = thisMR(all_ccs{rr,cc},:);
    end
 end
