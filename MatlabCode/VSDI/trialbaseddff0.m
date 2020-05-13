
function [img_hi_out,img_lo_out]=trialbaseddff0(img_hi,img_lo,img_no,nFrames,baseline_from,baseline_to)

img_x_size = size (img_hi,1);
img_y_size = size (img_hi,2);
img_hi = reshape(single(img_hi),img_x_size * img_y_size,nFrames);
img_lo = reshape(single(img_lo),img_x_size * img_y_size,nFrames);
img_no = reshape(single(img_no),img_x_size * img_y_size,nFrames);

img_hi_out = img_hi-img_no;
img_lo_out = img_lo-img_no;


% % % img_hi = img_hi ./ img_no;
% % % img_lo = img_lo ./ img_no;
% % % 
% % % img_hi_base = img_hi';
% % % img_hi_base = (mean(img_hi_base(baseline_from:baseline_to,:)))';
% % % img_hi_base = repmat(img_hi_base,1,nFrames);
% % % img_hi_out = ((img_hi - img_hi_base) ./ img_hi_base) * 100;
% % % 
% % % img_lo_base = img_lo';
% % % img_lo_base = (mean(img_lo_base(baseline_from:baseline_to,:)))';
% % % img_lo_base = repmat(img_lo_base,1,nFrames);
% % % img_lo_out = ((img_lo - img_hi_base) ./ img_lo_base) * 100;
   

% 
% img_ave_hi_out = (img_ave_hi_out / stimto)';
% img_ave_lo_out = (img_ave_lo_out / stimto)';
% 
% mn_lo = min(min(min(img_ave_lo_out)));
% x_lo = 10 - mn_lo;
% img_ave_lo_out = img_ave_lo_out + x_lo;
% 
% mn_hi = min(min(min(img_ave_hi_out)));
% x_hi = 10 - mn_hi;
% img_ave_hi_out = img_ave_hi_out + x_hi;
% 
% baseline_hi = (mean(img_ave_hi_out(:,baseline_from:baseline_to)'))';
% baseline_hi = repmat(baseline_hi,1,nframes);
% img_ave_hi_out = (( img_ave_hi_out - baseline_hi ) ./ baseline_hi) * 100;
% img_ave_hi_out = reshape(img_ave_hi_out, [img_x_size img_y_size nframes]);
% 
% baseline_lo = (mean(img_ave_lo_out(:,baseline_from:baseline_to)'))';
% baseline_lo = repmat(baseline_lo,1,nframes);
% img_ave_lo_out = (( img_ave_lo_out - baseline_lo ) ./ baseline_lo) * 100;
% img_ave_lo_out = reshape(img_ave_lo_out, [img_x_size img_y_size nframes]);
% 
% 
% % sigma = 2.5;
% % H = fspecial('gaussian',filtsize,sigma);
% % img_gauss01 = zeros (size(img_norm_01_ave));
% % img_gauss02 = zeros (size(img_norm_02_ave));
% % 
% % for i = 1:nframes
% %     img_gauss01(:,:,i) = imfilter(img_norm_01_ave(:,:,i),H);
% %     img_gauss02(:,:,i) = imfilter(img_norm_02_ave(:,:,i),H);
% % end
% % 
% % clear img_01;
% % clear img_no;
% % clear img_norm_01;
% % clear img_02;
% % clear img_norm_02;
% 
% for i = 1:nframes
%     img_ave_hi_out(:,:,i) = rot90(img_ave_hi_out(:,:,i));
%     img_ave_hi_out(:,:,i) = flipud(img_ave_hi_out(:,:,i));
%     img_ave_lo_out(:,:,i) = rot90(img_ave_lo_out(:,:,i));
%     img_ave_lo_out(:,:,i) = flipud(img_ave_lo_out(:,:,i));
% end
% 
% 
% filename_save = [pathname,'hi_DF_F0.raw'];
% fid = fopen(filename_save,'w', 'b');
% fwrite(fid, img_ave_hi_out, 'float32');
% fclose(fid);
% clear img_ave_hi_out;
% 
% filename_save = [pathname,'lo_DF_F0.raw'];
% fid = fopen(filename_save,'w', 'b');
% fwrite(fid, img_ave_lo_out, 'float32');
% fclose(fid);
% clear img_ave_lo_out;
% 
% 
% 
% 
