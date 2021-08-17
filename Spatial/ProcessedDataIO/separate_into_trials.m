function [RsC,mRsCT] = separate_into_trials(RsC,trials)

RsC = repmat(RsC,1,length(trials));
mRsCT = cell(size(RsC,1),length(trials));
for ii = 1:length(trials)
    [mRsCT(:,ii),~] = calc_mean_rasters(RsC(:,1),trials{ii});
end
RsC = find_responsive_rasters(RsC,1:10);