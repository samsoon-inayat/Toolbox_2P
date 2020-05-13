function PreProcessEV(file,filepath,savePath,Fs,stimulusFrame)
% Load df
filetype = '.tif';
filenamedfSeq = [filepath, file, filetype];
df = readTiff(filenamedfSeq); df = double(df);
% df(:,:,32) = (df(:,:,31)+df(:,:,33))/2;
[height,width,depth] = size(df);
% Load no trials and use it to calc F0

baselineFrames = 1:(stimulusFrame-1);
% baseline = df(:,baselineFrames);
% baseline = mean(baseline,2);
% baseline = repmat(baseline,1,depth);
% df = (df-baseline)./baseline;
% df(isnan(df))=0; df(isinf(df))=0;

df = reshape(df,height*width,depth);
parfor r = 1:(width*height)
    sig = df(r,:)';
    dsig = sig - locdetrend(sig,Fs,[0.3,0.1]);
    dsig = sig./dsig;
    Fo = mean(dsig(baselineFrames));
    df(r,:) = (dsig - Fo)/Fo;
%     figure(1);clf;%plot(sig);
%     hold on; plot(df(r,:),'r');plot(dsig,'m');
%     pause(0.1);
%     n = 0;
end

% baselineFrames = (stimulusFrame-25):(stimulusFrame-3);
% baseline = df(:,baselineFrames);
% baseline = mean(baseline,2);
% baseline = repmat(baseline,1,depth);
% df = (df-baseline)./baseline;
% df(isnan(df))=0; df(isinf(df))=0;
df = reshape(df,height,width,depth);
save([savePath '\' file '.mat'],'df','-v7.3');


