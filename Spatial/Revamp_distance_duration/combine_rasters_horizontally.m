function Rs = combine_rasters_horizontally(Rsi)

for ii = 1:size(Rsi,1)
    rasters = Rsi{ii,1}.sp_rasters1;
    for jj = 2:size(Rsi,2)
        thisR = Rsi{ii,jj};
        rasters = cat(2,rasters,thisR.sp_rasters1);
    end
    Rs{ii,1} = Rsi{ii,1};
    Rs{ii,1}.sp_rasters1 = rasters;
end