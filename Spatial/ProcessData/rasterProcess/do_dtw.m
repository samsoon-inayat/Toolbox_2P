function Rs = do_dtw(Rs)

for rr = 1:size(Rs,1)
    for cc = 1:size(Rs,2)
        R = Rs{rr,cc};
        rasters = permute(R.sp_rasters1,[1 2 3]);
        
    end
end

n = 0;