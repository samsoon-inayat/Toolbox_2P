function Rs = combine_rasters_conditions(Rsi)

for ii = 1:size(Rsi,1)
    rasters = Rsi{ii,1}.sp_rasters;
    dur = Rsi{ii,1}.duration;
    for jj = 2:size(Rsi,2)
        thisR = Rsi{ii,jj};
        rasters = cat(1,rasters,thisR.sp_rasters);
        dur = cat(1,dur,thisR.duration);
    end
    Rs{ii,1} = Rsi{ii,1};
    Rs{ii,1}.sp_rasters1 = rasters;
    Rs{ii,1}.duration1 = dur;
end