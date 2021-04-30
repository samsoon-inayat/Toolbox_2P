function [all_ccs,fracResp] = get_responsive_cells(rasters)

for rr = 1:size(rasters,1)
    for cc = 1:size(rasters,2)
        tRs = rasters{rr,cc};
        ccs = tRs.resp.p < 0.05;
        all_ccs{rr,cc} = ccs;
        fracResp(rr,cc) = sum(ccs)/length(ccs);
    end
end
