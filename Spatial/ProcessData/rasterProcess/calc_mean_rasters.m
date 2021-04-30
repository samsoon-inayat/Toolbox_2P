function mRs = calc_mean_rasters(rasters,trials)
if ~iscell(trials)
    for ii = 1:size(rasters,1)
        for cc = 1:size(rasters,2)
            mRsi{ii,cc} = (squeeze(nanmean(rasters{ii,cc}.rasters(trials,:,:),1)))';
        end
    end
    mRs = mRsi;
else
    for tt = 1:length(trials)
        mRs{tt} = calc_mean_rasters(rasters,trials{tt});
    end
end
