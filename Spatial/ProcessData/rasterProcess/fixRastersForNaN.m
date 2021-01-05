function rasters = fixRastersForNaN(rasters)

Rs = rasters.sp_rasters;
dur = rasters.duration;

lastBin = min(rasters.lastBin)-1;
dur = dur(:,1:lastBin);
Rs = Rs(:,1:lastBin,:);

lastBin1 = max(rasters.lastBin)-1;
dur1 = rasters.duration(:,1:lastBin1);
Rs1 = rasters.sp_rasters(:,1:lastBin1,:);
rasters.sp_rasters1 = Rs1;
rasters.duration1 = dur1;

if size(Rs,2) == 1
    rasters.sp_rasters_nan_corrected = Rs;
    rasters.duration_nan_corrected = dur;
    return;
end
n = 0;
parfor ii = 1:size(Rs,3)
    thisR = Rs(:,:,ii);
    thisR1 = Rs1(:,:,ii);
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
    
    thisR1 = fillmissing(thisR1,'linear',2,'EndValues','nearest');%inpaint_nans(thisR(rnan(jj),:),4);
    inds = thisR1<0;
    thisR1(inds) = 0;
    Rs1(:,:,ii) = thisR1;
end

rasters.sp_rasters_nan_corrected = Rs;
rasters.duration_nan_corrected = dur;
rasters.sp_rasters_nan_corrected1 = Rs1;
rasters.duration_nan_corrected1 = dur1;


