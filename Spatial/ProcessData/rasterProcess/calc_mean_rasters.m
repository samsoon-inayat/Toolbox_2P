function mRs = calc_mean_rasters(rasters,trials)
if ~iscell(trials)
    for ii = 1:length(rasters)
        mRs{ii,1} = (squeeze(nanmean(rasters{ii}.rasters(trials,:,:),1)))';
    end
else
    for ii = 1:length(rasters)
        for tt = 1:length(trials)
            mRs{ii,tt} = (squeeze(nanmean(rasters{ii}.rasters(trials{tt},:,:),1)))';
        end
    end
end
