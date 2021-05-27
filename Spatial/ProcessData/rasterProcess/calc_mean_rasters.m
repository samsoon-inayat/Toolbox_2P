function mRs = calc_mean_rasters(Rs,trials)
if ~iscell(trials)
    for ii = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            R = Rs{ii,cc};
            if strcmp(R.marker_name,'light22T') || strcmp(R.marker_name,'air55T')
                temp1 = (squeeze(nanmean(Rs{ii,cc}.fromFrames.sp_rasters(trials,:,:),1)))';
                mRsi{ii,cc} = normalizeSignal(temp1,2);
            end
            if strcmp(R.marker_name,'airD')
                temp1 = (squeeze(nanmean(Rs{ii,cc}.sp_rasters1(trials,:,:),1)))';
                mRsi{ii,cc} = normalizeSignal(temp1,2);
            end
        end
    end
    mRs = mRsi;
    
else
    for tt = 1:length(trials)
        mRs{tt} = calc_mean_rasters(Rs,trials{tt});
    end
end
