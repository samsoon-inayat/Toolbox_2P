function rasters = fixRastersForNaN(rasters)

Rs = rasters.sp_rasters;
dur = rasters.duration;

lastBin = min(rasters.lastBin)-1;
dur = dur(:,1:lastBin);
Rs = Rs(:,1:lastBin,:);

if size(Rs,2) == 1
    rasters.sp_rasters_nan_corrected = Rs;
    rasters.duration_nan_corrected = dur;
    return;
end
n = 0;
parfor ii = 1:size(Rs,3)
    thisR = Rs(:,:,ii);
%     [rnan,~] = find(isnan(thisR));
%     for jj = 1:length(rnan)
% %         try
%             thisR(rnan(jj),:) = fillmissing(thisR(rnan(jj),:),'linear','EndValues','nearest');%inpaint_nans(thisR(rnan(jj),:),4);
% %         catch
% %             n = 0;
% %         end
%     end
%     [rnan,~] = find(isnan(thisR));
% %     if ~isempty(rnan)
% %         n = 0;
% %     end
    thisR = fillmissing(thisR,'linear',2,'EndValues','nearest');%inpaint_nans(thisR(rnan(jj),:),4);
    inds = thisR<0;
    thisR(inds) = 0;
    Rs(:,:,ii) = thisR;
end


rasters.sp_rasters_nan_corrected = Rs;
rasters.duration_nan_corrected = dur;