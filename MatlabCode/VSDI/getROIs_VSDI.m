function [rois,roi_names] = getROIs (ags,gg,aa)

coord = get_ROI_coordinates(ags,gg,aa);
roi_names = fieldnames(coord); roi_names(1) = [];
pixels = get_ROI_pixels(coord,[2 2]);

for rr = 1:length(pixels);
    rois(rr).name = roi_names{rr};
    rois(rr).pixels = pixels{rr};
end