function o = get_props_Rs(Rs)

for rr = 1:size(Rs,1)
    for cc = 1:size(Rs,2)
        R = Rs{rr,cc}; 
        o.zMI{rr,cc} = R.info_metrics.ShannonMI_Zsh';
        [rs,MFR,centers,PWs] = get_gauss_fit_parameters(R.gauss_fit_on_mean,R.bin_width);
        o.rs{rr,cc} = rs'; o.MFR{rr,cc} = MFR'; o.centers{rr,cc} = centers'; o.PWs{rr,cc} = PWs';
        o.peak_locations{rr,cc} = R.peak_location';
    end
end


