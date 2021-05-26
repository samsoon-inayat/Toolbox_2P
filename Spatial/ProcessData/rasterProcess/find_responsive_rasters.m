function resp = find_responsive_rasters(Rs,th)

for rr = 1:size(Rs,1)
    for cc = 1:size(Rs,2)
        zMIs = Rs{rr,cc}.info_metrics.ShannonMI_Zsh;
        resp{rr,cc} = zMIs > th;
    end
end

    

