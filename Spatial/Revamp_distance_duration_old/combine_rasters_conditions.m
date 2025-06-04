function Rs = combine_rasters_conditions(Rsi)

for ii = 1:size(Rsi,1)
    for jj = 1:size(Rsi,2)
        thisR = Rsi{ii,jj};
        cols(ii,jj) = size(thisR.sp_rasters,2);
    end
end
mcols = max(cols,[],2);

for ii = 1:size(Rsi,1)
    rasters = Rsi{ii,1}.sp_rasters;
    dur = Rsi{ii,1}.duration;
    lb = Rsi{ii,1}.lastBin;
    tmcol = mcols(ii);
    if size(rasters,2) < tmcol
        diffcols = tmcol - size(rasters,2);
        rasters = padarray(rasters,[0 diffcols],0,'post');
        dur = padarray(dur,[0 diffcols],0,'post');
    end
    
    for jj = 2:size(Rsi,2)
        thisR = Rsi{ii,jj};
        if size(thisR.sp_rasters,2) < tmcol
            diffcols = tmcol - size(thisR.sp_rasters,2);
            thisR.sp_rasters = padarray(thisR.sp_rasters,[0 diffcols],0,'post');
            thisR.duration = padarray(thisR.duration,[0 diffcols],0,'post');
        end
        try
        rasters = cat(1,rasters,thisR.sp_rasters);
        dur = cat(1,dur,thisR.duration);
        lb = cat(1,lb,thisR.lastBin);
        catch
            n = 0;
        end
    end
    Rs{ii,1} = Rsi{ii,1};
    Rs{ii,1}.sp_rasters1 = rasters;
    Rs{ii,1}.duration1 = dur;
    Rs{ii,1}.lastBin = lb;
end