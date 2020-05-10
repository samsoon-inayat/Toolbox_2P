% mimg1 = ei1.ops1{1}.mimg1;
% mimg2 = ei2.ops1{1}.mimg1;
% 
% 
% tP{1} = ei1.tP;
% tP{2} = ei2.tP;
% dday = 1;
% for ii = 1:length(tP{dday}.stat)
%     cellz1(ii) = tP{dday}.stat(ii).iscell;
% end
% 
% [vals1,idxs1] = find(cellz1);
% 
% 
% dday = 2;
% for ii = 1:length(tP{dday}.stat)
%     cellz2(ii) = tP{dday}.stat(ii).iscell;
% end
% 
% [vals2,idxs2] = find(cellz2);
% 
% 
% for iii = 1:length(idxs1)
%     ii = idxs1(iii);
%     xs = tP{1}.stat(ii).xpix;
%     ys = tP{1}.stat(ii).ypix;
%     thisPixelIdxs1 = sub2ind(size(mimg1),ys,xs);
%     for jjj = 1:length(idxs2)
%         jj = idxs2(jjj);
%         xs = tP{2}.stat(jj).xpix;
%         ys = tP{2}.stat(jj).ypix;
%         thisPixelIdxs2 = sub2ind(size(mimg2),ys,xs);
%         sharedPixels = intersect(thisPixelIdxs1,thisPixelIdxs2);
%         overlap_indices(iii,jjj) = (length(sharedPixels))./(length(thisPixelIdxs1)+length(thisPixelIdxs2)-length(sharedPixels));
%         n = 0;
%     end
% end
% 
% rt


[overlapRs,overlapCs] = find(overlap_indices > 0.3);


b1 = ei1.b;
b2 = ei2.b;


cc = 1;
while 1
    roiM1 = zeros(size(mimg1));
    roiM2 = zeros(size(mimg2));
    cn1 = idxs1(overlapRs(cc));
    xs = tP{1}.stat(cn1).xpix;
    ys = tP{1}.stat(cn1).ypix;
    thisPixelIdxs1 = sub2ind(size(mimg1),ys,xs);
    roiM1(thisPixelIdxs1) = ones(size(thisPixelIdxs1));
    cn2 = idxs2(overlapCs(cc));
    xs = tP{2}.stat(cn2).xpix;
    ys = tP{2}.stat(cn2).ypix;
    thisPixelIdxs2 = sub2ind(size(mimg2),ys,xs);
    roiM2(thisPixelIdxs2) = ones(size(thisPixelIdxs2));
    
    figure(1);clf;
    subplot 221;
    imagesc(mimg1);axis equal
    subplot 223;
    imagesc(roiM1);axis equal
    subplot 222;
    imagesc(mimg2);axis equal
    subplot 224;
    imagesc(roiM2);axis equal
    
    thisSignal = tP{1}.signals(cn1,:);
    b = b1;
    [xVals,ccSignal] = getTrialSignalsDist(thisSignal,b,b.air_onset,b.air_offset);
    figure(100);clf;
    subplot 221;
    plot(xVals,ccSignal');
    xlim([min(xVals) max(xVals)]);
    title(sprintf('%d - %d',cc,cn1));
    subplot 223;
    plot(xVals,normA(nanmean(ccSignal)));
    xlim([min(xVals) max(xVals)]);
    subplot 224;
    imagesc(normA(nanmean(ccSignal)));colorbar
    title(sprintf('%d - %d',cc,cn1));
    mCCsignal = max(ccSignal,[],2);
    norm_ccSignal = [];
    for mm = 1:length(mCCsignal)
        norm_ccSignal(mm,:) = ccSignal(mm,:);%./mCCsignal(mm);
    end
    subplot 222;
    imagesc(norm_ccSignal);
    colorbar
    plotRawTrials(xVals,ccSignal,10,150);

    
    
    thisSignal = tP{2}.signals(cn2,:);
    b = b2;
    [xVals,ccSignal] = getTrialSignalsDist(thisSignal,b,b.air_onset,b.air_offset);
    figure(200);clf;
    subplot 221;
    plot(xVals,ccSignal');
    xlim([min(xVals) max(xVals)]);
    title(sprintf('%d - %d',cc,cn1));
    subplot 223;
    plot(xVals,normA(nanmean(ccSignal)));
    xlim([min(xVals) max(xVals)]);
    subplot 224;
    imagesc(normA(nanmean(ccSignal)));colorbar
    title(sprintf('%d - %d',cc,cn1));
    mCCsignal = max(ccSignal,[],2);
    norm_ccSignal = [];
    for mm = 1:length(mCCsignal)
        norm_ccSignal(mm,:) = ccSignal(mm,:);%./mCCsignal(mm);
    end
    subplot 222;
    imagesc(norm_ccSignal);
    colorbar
    plotRawTrials(xVals,ccSignal,10,250);

    
    
    keyVal = getkey;%pause;
    if keyVal == 27
        break;
    end
    if keyVal == 30 % Up arrow
        if cc <= size(overlapRs,1)
            cc = cc + 1;
        end
    end
    if keyVal == 31 % down arrow
        if cc >= 1
            cc = cc - 1;
        end
    end
end
