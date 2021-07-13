function [mRs,Rs] = calc_mean_rasters(Rs,trialsi)
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
                    strcmp(R.marker_name,'air77T')|| strcmp(R.marker_name,'air33T') || strcmp(R.marker_name,'airT')...
                    || strcmp(R.marker_name,'motionOnsets') || strcmp(R.marker_name,'motionOffsets')
%                 temp1 = (squeeze(nanmean(Rs{ii,cc}.fromFrames.sp_rasters(trials,:,:),1)))';
                temp1 = (squeeze(nanmean(Rs{ii,cc}.sp_rasters1(trials,:,:),1)))';
                mRsi{ii,cc} = normalizeSignal(temp1,2);
            end
            if strcmp(R.marker_name,'airD') || strcmp(R.marker_name,'beltD')
                temp1 = (squeeze(nanmean(Rs{ii,cc}.sp_rasters1(trials,:,:),1)))';
                mRsi{ii,cc} = normalizeSignal(temp1,2);
            end
            if length(trials) < 5
                Rs{ii,cc} = reduce_R(Rs{ii,cc},trials);
            end
        end
    end
    mRs = mRsi;
else
    for tt = 1:length(trialsi)
        mRs{tt} = calc_mean_rasters(Rs,trialsi{tt});
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

