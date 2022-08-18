function [mRs,Rs,mRs1] = calc_mean_rasters(Rs,trialsi)
if ~iscell(trialsi)
    for ii = 1:size(Rs,1)
        for cc = 1:size(Rs,2)
            R = Rs{ii,cc};
            if isempty(trialsi)
                trials = 1:size(R.sp_rasters,1);
            else
                trials = trialsi;
            end
            if strcmp(R.marker_name,'light22T') || strcmp(R.marker_name,'air55T') || strcmp(R.marker_name,'airIT') || ...
                    strcmp(R.marker_name,'air77T')|| strcmp(R.marker_name,'air33T') || strcmp(R.marker_name,'air44T') || strcmp(R.marker_name,'tone22T') ...
                    || strcmp(R.marker_name,'airT') || strcmp(R.marker_name,'airRT') || strcmp(R.marker_name,'airIRT') ...
                    || strcmp(R.marker_name,'motionOnsets') || strcmp(R.marker_name,'motionOffsets') || strcmp(R.marker_name,'airOnsets22T') || strcmp(R.marker_name,'airOffsets22T')...
                    || strcmp(R.marker_name,'airOnsets55T') || strcmp(R.marker_name,'airOffsets55T') || strcmp(R.marker_name,'airOffsets22_C') || strcmp(R.marker_name,'airOnsets22_C')...
                    || strcmp(R.marker_name,'airOnsets22P') || strcmp(R.marker_name,'airOnsets11T') || strcmp(R.marker_name,'airOffsets11T') || strcmp(R.marker_name,'light11T')
%                 temp1 = (squeeze(nanmean(Rs{ii,cc}.fromFrames.sp_rasters(trials,:,:),1)))';
                temp1 = (squeeze(nanmean(Rs{ii,cc}.sp_rasters1(trials,:,:),1)))';
                mRsi{ii,cc} = normalizeSignal(temp1,2);
                mRsi1{ii,cc} = temp1;
                % I want to check if normalization makes a difference in
                % later calculations of population vector correlation and
                % remapping calculcations
%                 temp2 = mRsi{ii,cc};
%                 [ctemp1,~] = corr(temp1); % find correlation
%                 [ctemp2,~] = corr(temp2); % find correlation
%                 figure(1000);clf;subplot 211;imagesc(temp1);colorbar;subplot 212;imagesc(temp2);colorbar;
%                 figure(1000);clf;subplot 211;imagesc(ctemp1);colorbar;subplot 212;imagesc(ctemp2);colorbar;
            end
            if strcmp(R.marker_name,'airD') || strcmp(R.marker_name,'beltD') || strcmp(R.marker_name,'airID')
                temp1 = (squeeze(nanmean(Rs{ii,cc}.sp_rasters1(trials,:,:),1)))';
                mRsi{ii,cc} = normalizeSignal(temp1,2);
                mRsi1{ii,cc} = temp1;
            end
            if length(trials) < 5
                Rs{ii,cc} = reduce_R(Rs{ii,cc},trials);
            end
        end
    end
    mRs = mRsi;
    mRs1 = mRsi1;
else
    for tt = 1:length(trialsi)
        [mRs{tt},~,mRs1{tt}] = calc_mean_rasters(Rs,trialsi{tt});
    end
end

function R = reduce_R(R,trials)
fields = fieldnames(R);
num_trials = size(R.sp_rasters,1);
for ii = 1:length(fields)
    thisField = fields{ii};
    cmdTxt = sprintf('TT = R.%s;',thisField); eval(cmdTxt);
    if size(TT,1) == num_trials
        if size(TT,3) > 1
            TT = TT(trials,:,:);
        else
            if size(TT,2) > 1
                TT = TT(trials,:);
            else
                TT = TT(trials,1);
            end
        end
        cmdTxt = sprintf('R.%s = TT;',thisField); eval(cmdTxt);
    end
end

