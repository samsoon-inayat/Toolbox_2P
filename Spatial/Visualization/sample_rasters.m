function ff = sample_rasters(R,ccs,ff)

if iscell(R)
    for cns = 1:length(ccs)
        tff = ff;
        for rr = 1:length(R)
            tff.h_axes = ff.h_axes(cns,rr);
            tff.axesPos = ff.axesPos(cns,rr);
            sample_rasters(R{rr},ccs(cns),tff);
        end
    end
    return;
end

n = 0;
if isfield(R,'fromFrames')
    ff = plot_time_rasters(R,ccs,ff);
else
    ff = plot_dist_rasters(R,ccs,ff);
end
cm = colormap(gray);
cm = flipud(cm(1:size(cm,1),:));
colormap(cm);


