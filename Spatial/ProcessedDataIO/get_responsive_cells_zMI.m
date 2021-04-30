function [all_ccs,fracResp] = get_responsive_cells_zMI(rasters)

for rr = 1:size(rasters,1)
    for cc = 1:size(rasters,2)
        tRs = rasters{rr,cc};
        try
            ccs = tRs.zMIs > 1.96;
        catch
            ccs = tRs.resp.zMIs > 1.96;
        end
        all_ccs{rr,cc} = ccs;
        fracResp(rr,cc) = sum(ccs)/length(ccs);
    end
end
