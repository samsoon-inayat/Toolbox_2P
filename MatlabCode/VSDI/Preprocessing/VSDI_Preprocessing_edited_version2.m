clear all;  
clc
tic;
for fn=1:4
nFrames =90000;
%img_stk = imreadallraw('F:\Data\January 20_2015\02_Spon_200Hz_bilateral_Urethane_20 min after injection.raw',128,128,96000,'*uint16');
camera_srate = 100;

% file = 'Mask_Mean_G1.tif';
%change path IMPORTANT
pathname = 'E:\Users\surjeet.singh\Documents\20160407\';
file={'02_spon_90000fr_100Hz_raw.raw','09_spon_90000fr_100Hz_raw.raw','11_spon_90000fr_100Hz_raw.raw','14_spon_90000fr_100Hz_raw.raw'};
filename = [pathname, char(file(fn))];


img_stk = imreadallraw(filename,128,128,nFrames,'*uint16');
img_stk = single(img_stk);

%change path IMPORTANT
Mask = imread('E:\Users\surjeet.singh\Documents\20160407\Mask.tif','tiff');

% Mask(Mask == 255) = 0;
% Mask = Mask / 254;
Mask = single(Mask);
NonZeroPixelsIndex = find(reshape(Mask,128*128,1) == 1);
%%
% for ii =1:nFrames
%     img_stk(:,:,ii) = img_stk(:,:,ii) .* Mask;
% end    
rows = size(img_stk,1);
columns = size(img_stk,2);
trim = 0;
% img_stk = imreadallraw('F:\Data\January 16_2015_VSDI_LFP\VSDI\Preprocessed\03_spon_96000fr_200Hz_PCA_detrended_deharmonized_dff0.raw',57,57,96000,'*float32');
%X = reshape(img_stk,128*128, nFrames);
X = reshape(img_stk,rows*columns, nFrames);
X = X';
% X = double(X);
X = single(X);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SVD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
[U, S, V] = svd(X(:,NonZeroPixelsIndex),0);
% toc
% [U_filt S_filt V_filt] = svd(X_filtered_trimed,0);
% [dU1, dS, dV] = svd(dXrecon,0);
% [Urecon, Srecon, Vrecon] = svd(Xrecon,0);
% [U_T, S_T, V_T] = svd(T,0);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% plot singular values on semilogarithmic scale %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; semilogy(diag(S)); grid on;
% figure; semilogy(diag(dS)); grid on;
% figure; semilogy(diag(S_filt)); grid on;
title('Singular values spectrum');
xlabel('components');
ylabel('singular value');
% removing the harmonics
% U_noharm = rmlinesc(U(:,1:20), params);
% dU_noharm = rmlinesmovingwinc(dU1(:,1:30),movingwin,10,params);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Spectrogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params.Fs = 2000;
% params.tapers=[3 5];
% params.fpass=[0.1 5];
% params.pad=6;
% movingwin=[1 0.5];
% 
% [Specgram,t,freq]=mtspecgramc(reshape(img_stk_recon(50,40,:),1,96000),movingwin,params);
% [Specgram,t,freq]=mtspecgramc(dXrecon(:,1000),movingwin,params);
% figure; plot_matrix(Specgram,t,freq);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Frequency content of principal components %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Spec, f] = mtspectrumc(U(:,1:32), params);
% [Spec_noharm, f_noharm] = mtspectrumc(U_noharm(:,1:29), params);
% % plot the spectrum
% figure; pcolor([1:29], f, 10*log10(Spec)); shading flat;
% title('Frequency content of components');
% xlabel('components');
% ylabel('frequency (Hz)');
% 
% figure; pcolor([1:29], f_noharm, 10*log10(Spec_noharm)); shading flat;
% title('Frequency content of components');
% xlabel('components');
% ylabel('frequency (Hz)');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Deharmonization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U_noharm = rmlinesc(U(:,1:M),params);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% reconstruction of data matrix using first M eigenvectors %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M =50;
% X_filtered_01_4_recon = U_filtered_trimed(:,1:M) * S(1:M,1:M) * V(:,1:M)';
% Xrecon = U_IMF1_test * S(1:M,1:M) * V(:,1:M)';
% Xrecon = U_IMF1(:,1:M) * S(1:M,1:M) * V(:,1:M)';
Xrecon = U(:,1:M) * S(1:M,1:M) * V(:,1:M)';
clear U S V X img_stk
%Xrecon_noharm = dU_noharm(:,1:M) * dS(1:M,1:M) * dV(:,1:M)';
%Xrecon_slow = dU_slow(:,1:M) * dS(1:M,1:M) * dV(:,1:M)';
%X_slow = U(:,1:size(S,1)) * S(1:size(S,1),1:size(S,1)) * V(:,1:size(S,1))'; 
% X_filtered_01_4_recon = U_filt(:,1:M) * S_filt(1:M,1:M) * V_filt(:,1:M)';
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Detrending %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Xrecon orginal (reconstructed signal after PCA/SVD)
%dXrecon is orginal minus the trend
dXrecon = locdetrend(Xrecon,camera_srate,[2 1]);
% figure;
% plot(Xrecon(:,NonZeroPixelsIndex(1000)) - dXrecon(:,NonZeroPixelsIndex(1000)),'r');
% hold on; plot(Xrecon(:,NonZeroPixelsIndex(1000)));
%%
% dU1 = locdetrend(U_noharm(:,1),200,[5 2.5]);
% figure; plot(0:0.005:95999*0.005, U_noharm(:,1) - dU1, 'r'); grid on; hold on;
% plot(0:0.005:95999*0.005, U_noharm(:,1), 'g');
%dX_noharm = locdetrend(Xrecon_noharm ,150,[17 1]);
dXrecon_dff0 = single(zeros(nFrames,rows*columns));
dXrecon_dff0(:,NonZeroPixelsIndex) = 100 * dXrecon./(Xrecon-dXrecon);
% dXrecon_dff0 = single(dXrecon_dff0);
% figure; plot(mean(dXrecon_dff0(:,NonZeroPixelsIndex),2)); grid on
clear dXrecon Xrecon
% T = Xrecon - dXrecon;
%T_noharm = Xrecon_noharm - dX_noharm;
%dU = locdetrend(U(:,1:20) ,150,[.1 .05]);
% make pseudo-color plot of matrix
%Y = double(dX);
%figure; pcolor(dX); shading flat;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   nFrames = number of frames
%   FR = Sampling rate
%   Lowpasscutoff = Lower band freq. in Hz (eg. 0.1 Hz); 
%   Highpasscutt = Higher end freq in Hz (eg. 6 Hz)

% FR = 200; 
%ImageSize = 128;
% Lowpasscutoff = 0.1; 
% Highpasscutoff = 4;

% pathname = 'E:\Data\';
% file = 'M101811_03_sponVSD_150Hz_5000fr.raw';
% filename = [pathname file];
% img = imreadallraw(filename,ImageSize,ImageSize,nFrames,'*uint16');

% Nyquist = FR/2;    
% Bandpass = [Lowpasscutoff Highpasscutoff]/Nyquist;
% FILTERORDER = 4;
% [z,p,k] = cheby1(FILTERORDER,0.1,Bandpass); %bandpass gives you 0 mean
% [sos,g] = zp2sos(z,p,k);

%U_slow = zeros(size(dU,1),M);
% U_filtered = [];
% % wait_bar = waitbar(0,'two-way chebyshev Band-pass filter is applied...');
% for ii = 1:M
% %       U_filtered(:,ii) = filtfilt(sos,g,double(U(:,ii)));
%       U_filtered(:,ii) = filtfilt(filter_coefficients,1,U_test(:,ii));
%       waitbar(ii/M, wait_bar);
% end
% close(wait_bar);
% U_filtered_trimed = U_filtered(trim+1:nFrames - trim,:);

filter_coefficients = load('Filter_FIR_05_6_100Hz_order400.txt');
% filter_coefficients = load('Filter_FIR_InverseSincLowPass_Fc6Hz_Fs20Hz_order5.txt');
dXrecon_dff0_filtered = zeros(nFrames - 2*trim, rows*columns);
for ii = 1:size(dXrecon_dff0,2)
     if ~all(dXrecon_dff0(:,ii) == 0)
%       X_filtered(:,ii) = filtfilt(sos,g,X(:,ii));
%       X_filtered(:,ii) = filtfilt(filter_coefficients,1,X(:,ii));  
        dXrecon_dff0_filtered(:,ii) = filtfilt(filter_coefficients,1,double(dXrecon_dff0(:,ii)));
     end  
end
dXrecon_dff0_filtered = single(dXrecon_dff0_filtered);
clear dXrecon_dff0
% X_filtered_trimed = X_filtered(trim+1:nFrames - trim,:);


% VSDI_filtered = filtfilt(sos,g,double(reshape(mean(mean(img_stk(4:5,20:21,:))),1,96000)));
% VSDI_filtered = filtfilt(sos,g,double(reshape(mean(mean(img_stk(17:18,19:20,:))),1,96000)));
% CTX_LFP_filtered_7_16 = filtfilt(sos,g,eeg(:,4));
% HPC_LFP_filtered = filtfilt(sos,g,eeg(:,3));
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% saving the preprocessed stack %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 06_spon_90000fr_100Hz_raw
% imwriteallraw('F:\Data\January 16_2015_VSDI_LFP\VSDI\Preprocessed\04_spon_96000fr_200Hz_filtered_05_20_PCA_dff0.raw',reshape(X_filtered_01_4_recon_dff0',57,57,95000),'*float32');
% imwriteallraw('F:\Data\20150401\Preprocessed\07_spon_96000fr_200hz_80by74_PCA_peakremoved_filtered_01_4_dff0_501_95500fr.raw',reshape(X_filtered_01_4_recon_dff0',rows,columns,size(U_filtered_trimed,1)),'*float32');
% imwriteallraw('C:\Users\surjeet.singh\Documents\MATLAB\Data for analysis\20161104\06_spon_90000fr_100Hz_detrend_dff0.raw',reshape(dXrecon_dff0',rows,columns,size(dXrecon_dff0,1)),'*float32');

% imwriteallraw('E:\Data\20160819\Preprocessed\00_spon_100Hz_10000fr_PCA_detrend_dff0_LowPassFiltered_Below5Hz.raw',reshape( dXrecon_dff0_filtered',rows,columns,size(dXrecon_dff0_filtered,1)),'*float32');

% tic
newfile={'02_spon_90000fr_100Hz_detrend_dff0_filtered_05_6Hz.raw','09_spon_90000fr_100Hz_detrend_dff0_filtered_05_6Hz.raw','11_spon_90000fr_100Hz_detrend_dff0_filtered_05_6Hz.raw', '14_spon_90000fr_100Hz_detrend_dff0_filtered_05_6Hz.raw'};
filtered_file = [pathname, char(newfile(fn))];
imwriteallraw(filtered_file,reshape( dXrecon_dff0_filtered',rows,columns,size(dXrecon_dff0_filtered,1)),'*float32');

end
TimeSpent = toc;
% imwriteallraw('F:\Data\20150410\Preprocessed\04_spon_96000fr_200hz_70by76_PCA_peakremoved_detrend_96000fr.raw',reshape(dXrecon_dff0_filtered',rows,columns,size(dXrecon_dff0,1)),'*float32');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% for ii =1:4
%     subplot(2,2,ii);
%     %plot(X(:,600*ii + 30),'g');
%     plot(Xrecon(:,600*ii + 30));
%     hold on
%     plot(T(:,600*ii + 30),'r');
%     grid on
% end
% 
% trend = zeros(nFrames,rows*columns);
% for ii = 1:size(Xrecon,2)
%     if ~all(Xrecon(:,ii) == 0)
%         a = 1:nFrames;
%         [maxs, maxs_locs] = findpeaks(Xrecon(:,ii));
%         [mins, mins_locs] = findpeaks(-Xrecon(:,ii));
%         mins = -mins;
%         
%         % figure; plot(maxs_locs,maxs,a(~ismember(1:96000,maxs_locs)),maxs_envelope)
%         b = zeros(nFrames,1);
%         b(maxs_locs) = maxs;
%         b(a(~ismember(1:nFrames,maxs_locs))) = interp1(maxs_locs,maxs,a(~ismember(1:nFrames,maxs_locs)),'pchip');
%         
%         c = zeros(nFrames,1);
%         c(mins_locs) = mins;
%         c(a(~ismember(1:nFrames,mins_locs))) = interp1(mins_locs,mins,a(~ismember(1:nFrames,mins_locs)),'pchip');
%         
%         trend(:,ii) = (b+c)/2;
%     end
%     % figure; plot(1:96000,X(:,3000),1:96000,trend,'r');grid on;
% end
% dXrecon_dff0 = 100 * (Xrecon - trend)./trend;
% 
% 
% U_IMF1 = zeros(nFrames,M);
% for ii = 1:M
% %     if ~all(Xrecon(:,ii) == 0)
%         a = 1:nFrames;
%         [maxs, maxs_locs] = findpeaks(U(:,ii));
%         [mins, mins_locs] = findpeaks(-U(:,ii));
%         mins = -mins;
%         
%         % figure; plot(maxs_locs,maxs,a(~ismember(1:96000,maxs_locs)),maxs_envelope)
%         b = zeros(nFrames,1);
%         b(maxs_locs) = maxs;
%         b(a(~ismember(1:nFrames,maxs_locs))) = interp1(maxs_locs,maxs,a(~ismember(1:nFrames,maxs_locs)),'pchip');
%         
%         c = zeros(nFrames,1);
%         c(mins_locs) = mins;
%         c(a(~ismember(1:nFrames,mins_locs))) = interp1(mins_locs,mins,a(~ismember(1:nFrames,mins_locs)),'pchip');
%         
%         U_IMF1(:,ii) = (b+c)/2;
% %     end
%     % figure; plot(1:96000,U(:,1),1:96000,U_IMF1(:,1),'r');grid on;
% end
% % dU_dff0 = 100 * (U(:,1:M) - trend)./trend;
% U_IMF2 = zeros(nFrames,M);
% for ii = 1:M
% %     if ~all(Xrecon(:,ii) == 0)
%         a = 1:nFrames;
%         [maxs, maxs_locs] = findpeaks(U_IMF1(:,ii));
%         [mins, mins_locs] = findpeaks(-U_IMF1(:,ii));
%         mins = -mins;
%         
%         % figure; plot(maxs_locs,maxs,a(~ismember(1:96000,maxs_locs)),maxs_envelope)
%         b = zeros(nFrames,1);
%         b(maxs_locs) = maxs;
%         b(a(~ismember(1:nFrames,maxs_locs))) = interp1(maxs_locs,maxs,a(~ismember(1:nFrames,maxs_locs)),'pchip');
%         
%         c = zeros(nFrames,1);
%         c(mins_locs) = mins;
%         c(a(~ismember(1:nFrames,mins_locs))) = interp1(mins_locs,mins,a(~ismember(1:nFrames,mins_locs)),'pchip');
%         
%         U_IMF2(:,ii) = (b+c)/2;
% %     end
%     % figure; plot(1:96000,X(:,3000),1:96000,trend,'r');grid on;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%% df/f_0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %dXrecon_base = repmat(mean(dXrecon(2000:end,:), 1),94000, 1);
% % dXrecon_base = repmat(mean(dXrecon),nFrames, 1) + repmat(Xrecon(1,:),nFrames, 1);
% X_filtered_01_4_recon_base = repmat(mean(X_filtered_01_4_recon),size(X_filtered_trimed,1), 1) + repmat(X(1,:),size(X_filtered_trimed,1), 1);
% X_filtered_01_4_recon_base = repmat(mean(X_filtered_01_4_recon),size(U_filtered_trimed,1), 1) + repmat(X(1,:),size(U_filtered_trimed,1), 1);
% 
% %dXrecon_dff0 = (dXrecon(2001:end,:) - repmat(mean(dXrecon(2001:end,:)),94000, 1)) ./ (repmat(mean(dXrecon(2001:end,:),94000, 1) + repmat(dXrecon(2001,:),94000, 1)));
% % dXrecon_dff0 = 100 * (dXrecon - repmat(mean(dXrecon),96000, 1))  ./ dXrecon_base;
% X_filtered_01_4_recon_dff0 = 100 * (X_filtered_01_4_recon - repmat(mean(X_filtered_01_4_recon),size(X_filtered_trimed,1), 1))  ./ X_filtered_01_4_recon_base;
% X_filtered_01_4_recon_dff0 = 100 * (X_filtered_01_4_recon - repmat(mean(X_filtered_01_4_recon),size(U_filtered_trimed,1), 1))  ./ X_filtered_01_4_recon_base;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%% Ploting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% for ii = 3:6
%     plot(0+19.27:0.005:95999*0.005+19.27,reshape(grid_stk_5(9-ii,ii,:),1,nFrames) + (2*ii-1)* 20 .* ones(1,nFrames));
%     hold on
%     plot(0+19.27:0.005:95999*0.005+19.27,reshape(grid_stk_5(ii,ii,:),1,nFrames) + 2*ii* 20 .* ones(1,nFrames));
%     hold on
%     grid on
% end
% 
% img_stk_recon = reshape(dXrecon',128,128,96000);
% figure;
% for ii = 8:10:128
%     %plot(0:0.005:95999*0.005,reshape(img_stk_recon(136-ii,ii,:),1,nFrames) + (2*ii-1) .* ones(1,nFrames));
%     hold on
%     plot(0:0.005:95999*0.005,reshape(img_stk_recon(ii,ii,:),1,nFrames) + 2*ii .* ones(1,nFrames));
%     hold on
%     grid on
% end
% 
% 
% figure;
% for ii = 1:25
%     %plot(0:0.005:95999*0.005,reshape(img_stk_recon(136-ii,ii,:),1,nFrames) + (2*ii-1) .* ones(1,nFrames));
%     hold on
%     plot(0:0.005:95999*0.005,reshape(img_stk_recon(5*ii,40,:),1,nFrames) + 2*ii*20 .* ones(1,nFrames));
%     hold on
%     grid on
% end
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure;
% plot(0+19.2756:0.005:95999*0.005+19.2756,reshape(grid_stk_5(4,4,:),1,nFrames) +  20 .* ones(1,nFrames));
% grid on
% hold on
% plot(0:0.0005:1020671*0.0005,eeg(:,1))
% 
% figure;
% ax1 = subplot(3,1,1);
% plot(0+19.2756:0.005:95999*0.005+19.2756,reshape(mean(mean(img_stk_filtered_01_8(5:6,6:7,:))),1,96000));
% grid on
% ax2 = subplot(3,1,2);
% plot(0:0.0005:1020671*0.0005,eeg(:,1))
% grid on
% ax3 = subplot(3,1,3);
% plot(0:0.0005:1020671*0.0005,eeg(:,3))
% grid on
% linkaxes([ax1 ax2 ax3], 'x')
% 
% figure;
% ax1 = subplot(2,1,1);
% plot(0+19.2756:0.005:95999*0.005+19.2756, reshape(mean(mean(img_stk_reduced(5:6,6:7,:))),1,96000));
% grid on
% ax2 = subplot(2,1,2);
% plot_matrix(Spec,tvec+19.2756*ones(1,size(tvec,2)),f);
% grid on
% linkaxes([ax1 ax2], 'x')
% 
% 
% 
% 
% grid_stk_2 = [];
% n = 2;
% for jj = 1:nFrames
%     for kk = 1:n:57-n
%         for ll =1:n:57-n
%             grid_stk_2((kk+n-1)/n,(ll+n-1)/n,jj) = mean(mean(img_stk(kk:kk+n-1,ll:ll+n-1,jj)));
%         end
%     end
% end    
% 
% 
% 
% 
% for ii = 3:4
%     
%     X = reshape(imreadallraw(['F:\Data\20150508\' '0' num2str(ii) '_spon_96000fr_200Hz_94by72.raw'],72,94,96000,'*uint16'), 94*72,96000);
%     X = double(X');
%     [U, S, V] = svd(X,0);
%     save(['F:\Data\20150508\Preprocessed\SVD_' num2str(ii) '.mat'], 'U', 'S', 'V', '-v7.3')
%     
%     clear U S V X
% end
