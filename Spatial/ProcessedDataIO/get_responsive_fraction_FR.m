function [resp,respei] = get_responsive_fraction_FR(ei,Rs,sp_threshold)

% for rr = 1:size(Rs,1)
%     for cc = 1:size(Rs,2)
%         R = Rs{rr,cc};
%         rasters = permute(R.sp_rasters1,[2 1 3]);
%         sR = sum(squeeze(nansum(rasters,1))>0);
%         R.resp.FR_based = (sR)>3;
%         Rs{rr,cc} = R;
%     end
% end

% resp_exc_inh = NaN;

for rr = 1:size(Rs,1)
%     respT = zeros(size(Rs{rr,1}.resp.FR_based'));
    for cc = 1:size(Rs,2)
        resp{rr,cc} = Rs{rr,cc}.resp.FR_based';
    end
end

for ii = 1:length(ei)
    spSignals = ei{ii}.plane{1}.tP.deconv.spSigAll;
    if length(ei{ii}.plane) > 1
        spSignals = [spSignals;ei{ii}.plane{2}.tP.deconv.spSigAll];
    end
    minSp = min(spSignals,[],2);
    maxSp = max(spSignals,[],2);
    respei{ii,1} = maxSp < sp_threshold;
    all_maxSp(ii) = max(maxSp);
end
respei = repmat(respei,1,size(Rs,2));


