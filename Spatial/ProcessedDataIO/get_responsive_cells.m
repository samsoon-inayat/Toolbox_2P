function [all_ccs,fracResp] = get_responsive_cells(rasters,thr)

for rr = 1:size(rasters,1)
    for cc = 1:size(rasters,2)
        tRs = rasters{rr,cc};
        try
            ccs = tRs.resp.p < thr(1);
        catch
            ccs = tRs.resp.zMIs > thr(2);
        end
        all_ccs{rr,cc} = ccs;
        fracResp(rr,cc) = sum(ccs)/length(ccs);
    end
end
