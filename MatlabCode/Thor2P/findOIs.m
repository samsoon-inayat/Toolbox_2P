function ois = findOIs (ei1,ei2)
day1 = ei1.exp_date;
day2 = ei2.exp_date;
fileName1 = sprintf('overlap_indices_%s_%d_%s_%d.mat',day1,ei1.recording_number,day2,ei2.recording_number);
fileName2 = sprintf('overlap_indices_%s_%d_%s_%d.mat',day2,ei2.recording_number,day1,ei1.recording_number);
fileName1 = makeName(fileName1,ei1.folders.animalFolder);
fileName2 = makeName(fileName2,ei1.folders.animalFolder);
if exist(fileName1,'file')
    ois = load(fileName1);
    return;
end
% if exist(fileName2,'file')
%     ois = load(fileName2);
%     ois.OIs = ois.OIs';
%     temp = ois.idxs2;
%     ois.idxs1 = temp.idxs1;
%     idxs2 = temp.idxs2;
%     return;
% end

disp('Please wait ... finding overlap indices ... this will take a while');

yrange1 = ei1.ops1{1}.yrange;
xrange1 = ei1.ops1{1}.xrange;

yrange2 = ei2.ops1{1}.yrange;
xrange2 = ei2.ops1{1}.xrange;

tP1 = ei1.tP;
tP2 = ei2.tP;
mimg1 = ei1.ops1{1}.mimg1;
mimg2 = ei2.ops1{1}.mimg1;

fixed = mimg1;
moving = mimg2;

tformEstimate = imregcorr(moving,fixed);
Rfixed = imref2d(size(fixed));

cellz1 = zeros(1,length(tP1.stat));
for ii = 1:length(tP1.stat)
    cellz1(ii) = tP1.stat(ii).iscell;
end
[~,idxs1] = find(cellz1);
cellz2 = zeros(1,length(tP2.stat));
for ii = 1:length(tP2.stat)
    cellz2(ii) = tP2.stat(ii).iscell;
end
[~,idxs2] = find(cellz2);

hWaitBar = waitbar(0,sprintf('Processing'));
zImg = zeros(size(moving));
thisPixelIdxs2R = cell(1,length(idxs2));
for jjj = 1:length(idxs2)
    jj = idxs2(jjj);
    xs = tP2.stat(jj).xpix + min(xrange2);
    ys = tP2.stat(jj).ypix + min(yrange2);
    thisPixelIdxs2 = sub2ind(size(mimg2),ys,xs);
    moving = zImg;
    moving(thisPixelIdxs2) = ones(size(thisPixelIdxs2));
    movingReg = imwarp(moving,tformEstimate,'OutputView',Rfixed);
    [yss,xss] = find(movingReg > 0);
    thisPixelIdxs2R{jjj} = sub2ind(size(mimg2),yss,xss);
    waitbar(jjj/length(idxs2),hWaitBar,sprintf('Processing %d of %d',jjj,length(idxs2)));
end
close(hWaitBar);

OIs = zeros(length(idxs1),length(idxs2));
for iii = 1:length(idxs1)
    ii = idxs1(iii);
    xs = tP1.stat(ii).xpix + min(xrange1);
    ys = tP1.stat(ii).ypix + min(yrange1);
    thisPixelIdxs1 = sub2ind(size(mimg1),ys,xs);
    for jjj = 1:length(idxs2)
%         jj = idxs2(jjj);
%         xs = tP2.stat(jj).xpix + min(xrange2);
%         ys = tP2.stat(jj).ypix + min(yrange2);
%         thisPixelIdxs2 = sub2ind(size(mimg2),ys,xs);
%         moving = zImg;
%         moving(thisPixelIdxs2) = ones(size(thisPixelIdxs2));
%         movingReg = imwarp(moving,tformEstimate,'OutputView',Rfixed);
%         [yss xss] = find(movingReg > 0);
%         thisPixelIdxs2 = sub2ind(size(mimg2),yss,xss);
        sharedPixels = intersect(thisPixelIdxs1,thisPixelIdxs2R{jjj});
        OIs(iii,jjj) = (length(sharedPixels))./(length(thisPixelIdxs1)+length(thisPixelIdxs2R{jjj})-length(sharedPixels));
%         waitbar(sub2ind([length(idxs2) length(idxs1)],jjj,iii)/(length(idxs1)*length(idxs2)),hWaitBar,sprintf('Processing cells %d - %d',iii,jjj));
    end
%     waitbar(iii/length(idxs1),hWaitBar,sprintf('Processing %d of %d',iii,length(idxs1)));
end
% close(hWaitBar);
ois.idxs1 = idxs1;
ois.idxs2 = idxs2;
ois.OIs = OIs;
ois.tformEstimate = tformEstimate;
ois.Rfixed = Rfixed;

save(fileName1,'-struct','ois');
