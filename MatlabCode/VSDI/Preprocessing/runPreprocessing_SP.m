function runPreprocessing_SP (ags,pags,spon_files_info,indices,sp_frames,maskFactor)

for ii = 1:length(indices)
    indx = indices(ii);
    preProcessFolder(ags,pags,spon_files_info,indx,sp_frames,maskFactor);
end


function preProcessFolder(ags,pags,spon_files_info,indx,sp_frames,maskFactor)
gg = spon_files_info.list(indx,2);
aa = spon_files_info.list(indx,3);
rr = spon_files_info.list(indx,4);
thisRecording = ags(gg).animals{aa}.sRecordings{rr};
pThisRecording = pags(gg).animals{aa}.sRecordings{rr};
disp(thisRecording.root_folder);
disp(pThisRecording.root_folder);
mask = getMask(ags(gg).animals{aa}.root_folder,maskFactor);

% coord = get_ROI_coordinates(ags,gg,aa);
% roi_names = fieldnames(coord); roi_names(1) = [];
% pixels = get_ROI_pixels(coord,[2 2]);
% rois.names = roi_names;
% rois.pixels = pixels;
disp('Loading data');
trials = {thisRecording.name};
notrials = trials;
if ischar(sp_frames)
    sp_frames = [1 thisRecording.props.number_of_frames];
end
noallImg = cell(1,length(trials));
root_folder = thisRecording.root_folder;
for ii = 1:length(trials)
    filenamedfSeq = makeName([trials{ii} '.tif'],root_folder);
    img = readTiff(filenamedfSeq,sp_frames);
    img = double(img);
    noallImg{ii}.img = img;
end

disp('Processing data');
noallImg = processTrials(notrials,noallImg,mask);

disp('Writing results to disk');
writeToDisk(notrials,noallImg,pThisRecording.root_folder);
% writeToDisk(lotrials,loallImg,pThisRecording.root_folder);
% writeToDisk(hitrials,hiallImg,pThisRecording.root_folder);

disp('Done with this folder');

function allImg = processTrials(trials,allImg,mask)
maskI = find(mask);
for kk = 1:length(trials)
    img = allImg{kk}.img;
    img = reshape(img,size(img,1)*size(img,2),size(img,3));
    mVals = img(maskI,:);
    dsig = findBleachingTrends(mVals');
    sig = mVals'./dsig;
    Fo = mean(sig(:));
    sig = (sig-Fo)./Fo;
    allImg{kk}.img = sig';
end

function writeToDisk (trials,allImg,root_folder)
for kk = 1:length(trials)
    fileName = makeName([trials{kk} '.mat'],root_folder);
    df = allImg{kk}.img;
    save(fileName,'df','-v7.3');
end


function df = findAverageOfTrials(trials,allImg)
df = zeros(size(allImg{1}.img));
for kk = 1:length(trials)
    df = df + allImg{kk}.img;
end
df = df/length(trials);


function dfSeq = readTiff(fileName,sp_frames)
% info = imfinfo(fileName);
% num_images = numel(info);
num_images = sp_frames(2)-sp_frames(1)+1;
dfSeq = zeros(128,128,num_images,'uint16');
remaind = mod(num_images,65536);
blocks = floor(num_images/65536);
for ii = 1:blocks
    shift = (ii - 1) * 65536;
    parfor k = 1:65536
        dfSeq(:,:,k + shift) = imread(fileName,'tif',k);
    end
end
shift = blocks * 65536;
parfor k = 1:remaind
    dfSeq(:,:,k + shift) = imread(fileName,'tif',(k+shift)-1+ sp_frames(1));
end


    
