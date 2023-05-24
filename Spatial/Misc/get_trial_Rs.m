function [allRsC,allmRsT,allmRsT_Raw] = get_trial_Rs(o,si,trialNs)

Rs = o.Rs(:,si);mR = o.mR(:,si); 
for cn = 1:length(si)
    trials = mat2cell([trialNs]',ones(size([trialNs]')));
    RsC = repmat(Rs(:,cn),1,length(trialNs));
    mRsCT = cell(size(RsC,1),length(trials)); mRsCT_Raw = mRsCT;
    for ii = 1:length(trials)
        ii;
        [mRsCT(:,ii),~,mRsCT_Raw(:,ii)] = calc_mean_rasters(RsC(:,1),trials{ii});
    end
    allmRsT{cn} = mRsCT; allmRsT_Raw{cn} = mRsCT_Raw;
     allRsC{cn} = RsC;
end
