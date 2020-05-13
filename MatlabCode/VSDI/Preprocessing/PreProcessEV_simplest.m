function dfbyfo = PreProcessEV_simplest(file,filepath,savePath,Fs,stimulusFrame)
trial = file(3:end);
% Load ImgSeq_raw
filetype = '.tif';
filenameImgSeq = [filepath, file, filetype];
ImgSeq_raw = readTiff(filenameImgSeq); ImgSeq_raw = double(ImgSeq_raw);
% ImgSeq_raw(:,:,32) = (ImgSeq_raw(:,:,31)+ImgSeq_raw(:,:,33))/2;
[width,height,~] = size(ImgSeq_raw);
% Load no trials and use it to calc F0
filetype = '.tif';
filenameNo = [filepath, 'no',trial, filetype];
F0 = readTiff(filenameNo); F0 = double(F0);
% parfor r = 1:width
%     for c = 1:height
%         sig = F0(r,c,:);sig =sig(:);
%         yhat = sig - locdetrend(sig,Fs,[0.3,0.1]);
%         F0(r,c,:) = yhat;
%     end
% end
dfbyfo = ImgSeq_raw./F0;
baselineFrames = (stimulusFrame-25):(stimulusFrame-3);
baseline = dfbyfo(:,:,baselineFrames);
baseline = mean(baseline,3);
baseline = repmat(baseline,1,1,size(ImgSeq_raw,3));
dfbyfo = (dfbyfo-baseline)./baseline;
% Calc DF/F0
% dfbyfo = 100*(ImgSeq_raw-F0)./F0; 
dfbyfo(isnan(dfbyfo))=0; dfbyfo(isinf(dfbyfo))=0;
% save the results
% imwriteallraw([savePath, file, '_DF_F0_LPF_25Hz_G1_G1.raw'],dfbyfo,'*float32');
save([savePath, 'dfbyfo.mat'],'dfbyfo','-v7.3');


function ImgSeq = readTiff(fileName)
info = imfinfo(fileName);
num_images = numel(info);
ImgSeq = zeros(info(1).Width,info(1).Height,num_images,'uint16');
for k = 1:num_images
    ImgSeq(:,:,k) = imread(fileName,'tif',k);
end