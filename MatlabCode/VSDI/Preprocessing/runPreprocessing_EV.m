function runPreprocessing_EV (ags,pags,stim_folders_info,indices,frames_info,maskFactor)

for ii = 1:length(indices)
    display(sprintf('Processing recording %d of %d',ii,length(indices)));
    indx = indices(ii);
    preProcessFolder(ags,pags,stim_folders_info,indx,frames_info,maskFactor);
end


function preProcessFolder(ags,pags,stim_folders_info,indx,frames_info,maskFactor)
gg = stim_folders_info.list(indx,2);
aa = stim_folders_info.list(indx,3);
rr = stim_folders_info.list(indx,4);
thisRecording = ags(gg).animals{aa}.eRecordings{rr};
pThisRecording = pags(gg).animals{aa}.eRecordings{rr};
disp(thisRecording.root_folder);
disp(pThisRecording.root_folder);
mask = getMask(pags(gg).animals{aa}.root_folder,maskFactor);

% rois = getROIs_VSDI(ags,gg,aa);

hitrials = thisRecording.hiTrials;
lotrials = thisRecording.loTrials;
notrials = thisRecording.noTrials;
disp('Loading data');
trials = notrials;
root_folder = thisRecording.root_folder;
noallImg = loadTrials(notrials,root_folder);
hiallImg = loadTrials(hitrials,root_folder);
loallImg = loadTrials(lotrials,root_folder);

disp('Processing data');
noallImg = processTrials(notrials,noallImg,frames_info,mask);
loallImg = processTrials(lotrials,loallImg,frames_info,mask);
hiallImg = processTrials(hitrials,hiallImg,frames_info,mask);

disp('Writing results to disk');
% writeToDisk(notrials,noallImg,pThisRecording.root_folder,rois,mask);
% writeToDisk(lotrials,loallImg,pThisRecording.root_folder,rois,mask);
% writeToDisk(hitrials,hiallImg,pThisRecording.root_folder,rois,mask);

writeToDiskAlldf(notrials,noallImg,pThisRecording.root_folder);
writeToDiskAlldf(hitrials,hiallImg,pThisRecording.root_folder);
writeToDiskAlldf(lotrials,loallImg,pThisRecording.root_folder);

disp('Finding averages ...');
df = findAverageOfTrials(hitrials,hiallImg);
df_hi = df;
% fileName = makeName('mean_hi.mat',pThisRecording.root_folder);
% save(fileName,'df','-v7.3');


df = findAverageOfTrials(lotrials,loallImg);
df_lo = df;
% fileName = makeName('mean_lo.mat',pThisRecording.root_folder);
% save(fileName,'df','-v7.3');


df = findAverageOfTrials(notrials,noallImg);
fileName = makeName('mean_no.mat',pThisRecording.root_folder);
save(fileName,'df','-v7.3');

df = (df_hi + df_lo)/2;
fileName = 'mean_hi_lo.mat';
fileName = makeName(fileName,pThisRecording.root_folder);
save(fileName,'df','-v7.3');

disp('Done with this folder');

function allImg = loadTrials(trials,root_folder)
parfor ii = 1:length(trials)
    filenamedfSeq = makeName([trials{ii} '.tif'],root_folder);
    img = readTiff(filenamedfSeq); img = double(img);
    allImg{ii}.img = img;
end

function allImg = processTrials(trials,allImg,frames_info,mask)
baseline_from = frames_info(2);
baseline_to = frames_info(3);
baselineFrames = baseline_from:baseline_to;
% startFrame = frames_info(4);
% endFrame = frames_info(5);
maskI = find(mask);
parfor kk = 1:length(trials)
    img = allImg{kk}.img;
    img = reshape(img,size(img,1)*size(img,2),size(img,3));
    mVals = img(maskI,:);
%     [U, S, V] = svd(mVals,0);
%     M = 150;
%     Xrecon = U(:,1:M) * S(1:M,1:M) * V(:,1:M)';
    dsig = findBleachingTrend(mVals');
%     dXrecon = locdetrend(Xrecon,150,[0.3 0.1]);
%     sig = 100 * dXrecon./(Xrecon-dXrecon);
    sig = mVals'./dsig;
    baselines = mean(sig(baselineFrames,:));
    baselines = repmat(baselines,size(sig,1),1);
    sig = 100*(sig-baselines)./baselines;
    sig(isnan(sig))=0; sig(isinf(sig))=0;
%     sig = sig(startFrame:endFrame,:);
    allImg{kk}.img = sig';
end

function writeToDiskAlldf (trials,allImg,root_folder)
for kk = 1:length(trials)
    df = allImg{kk}.img;
    fileName = makeName([trials{kk} '_df.mat'],root_folder);
    save(fileName,'df','-v7.3');
end

function writeToDisk (trials,allImg,root_folder,rois,mask)
pixels = rois.pixels;
for kk = 1:length(trials)
    dfbf1 = allImg{kk}.img;
    df1 = convertToImgSeq(dfbf1,mask);
    df1 = applySpatialFilter(df1);
    for jj = 1:length(pixels)
        thisPixels = pixels{jj};
        thisMask = makeMask(size(df1(:,:,1)),thisPixels);
        thisdfbyf = mean(getMaskedValues(df1,thisMask));
        df(jj,:) = thisdfbyf';
    end
    fileName = makeName([trials{kk} '.mat'],root_folder);
    save(fileName,'df','-v7.3');
end


function df = findAverageOfTrials(trials,allImg)
df = zeros(size(allImg{1}.img));
for kk = 1:length(trials)
    df = df + allImg{kk}.img;
end
df = df/length(trials);


function dfSeq = readTiff(fileName)
info = imfinfo(fileName);
num_images = numel(info);
dfSeq = zeros(info(1).Width,info(1).Height,num_images,'uint16');
parfor k = 1:num_images
    dfSeq(:,:,k) = imread(fileName,'tif',k);
end
