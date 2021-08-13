function Rs = filterRasters(Rs)

for rr = 1:size(Rs,1)
    for cc = 1:size(Rs,2)
        R = Rs{rr,cc};
        rasters = permute(R.sp_rasters1,[2 1 3]);
        [B,T,N] = size(rasters);
        rasters = reshape(rasters,B*T,N);
        fRasters = fast_smooth(rasters,2);
%         rasters = reshape(rasters,B,T,N);
        fRasters = reshape(fRasters,B,T,N);
%         rasters = permute(rasters,[2 1 3]);
        fRasters = permute(fRasters,[2 1 3]);
        R.sp_rasters1_unfiltered = R.sp_rasters1;
        R.sp_rasters1 = fRasters;
        Rs{rr,cc} = R;
%         figure(1000);clf;imagesc(squeeze(rasters(:,:,1)));
%         figure(1001);clf;imagesc(squeeze(fRasters(:,:,1)));
%         meanFR = nanmean(rasters,[2 3]);
    end
end