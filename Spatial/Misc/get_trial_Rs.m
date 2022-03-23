function [allRsC,allmRsT] = get_trial_Rs(o,si,trialNs)

Rs = o.Rs(:,si);mR = o.mR(:,si); 
for cn = 1:length(si)
    trials = mat2cell([trialNs]',ones(size([trialNs]')));
    RsC = repmat(Rs(:,cn),1,length(trialNs));
    mRsCT = cell(size(RsC,1),length(trials));
    for ii = 1:length(trials)
        ii;
        [mRsCT(:,ii),~] = calc_mean_rasters(RsC(:,1),trials{ii});
    end
    allmRsT{cn} = mRsCT;
    allRsC{cn} = RsC;
end
