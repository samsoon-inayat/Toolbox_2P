function [nbins,binWidth] = get_nbins_F(b,onsets,offsets,binWidth,rasterType)
for trial = 1:length(onsets)
    st = onsets(trial);
    se = offsets(trial);
    if strcmp(rasterType,'time')
        trial_vals(trial) = b.ts(se)-b.ts(st);
    end
    if strcmp(rasterType,'dist')
        trial_vals(trial) = b.dist(se)-b.dist(st);
    end
end

max_trial_vals = max(trial_vals);
nbins = ceil(max_trial_vals/binWidth);