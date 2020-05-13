function visualizeTrials_SP (ags,pags,spon_files_info,indices,maskFactor)
indx = indices(1);
gg = spon_files_info.list(indx,2);
aa = spon_files_info.list(indx,3);
rr = spon_files_info.list(indx,4);
% thisRecording = ags(gg).animals{aa}.eRecordings{rr};
pThisRecording = pags(gg).animals{aa}.sRecordings{rr};

coord = get_ROI_coordinates(ags,gg,aa);
roi_names = fieldnames(coord); roi_names(1) = [];
pixels = get_ROI_pixels(coord,[2 2]);
rois.names = roi_names;
rois.pixels = pixels;
mask = getMask(ags(gg).animals{aa}.root_folder,maskFactor);
masko = getMask(ags(gg).animals{aa}.root_folder,1);
if isempty(mask)
    display('Mask does not exist');
    return;
end

fileName = makeName('dfbyfo.mat',pThisRecording.root_folder);
try
    temp = load(fileName);
catch
    fileName = makeName([pThisRecording.name '.mat'],pThisRecording.root_folder);
    temp = load(fileName);
    temp.dfbyfo = temp.df;
end

% temp.dfbyfo = applyTemporalFilter(temp.df);
% temp.dfbyfo = applySpatialFilter(temp.dfbyfo);
code = visualizeImgSeq(temp.dfbyfo,mask,masko,pThisRecording,0,rois);


function mainQuit = visualizeImgSeq (dfbf1,mask,masko,info,trn,rois)
dfbf = zeros(size(mask,1),size(mask,2),size(dfbf1,2));
dfbf = reshape(dfbf,size(mask,1)*size(mask,2),size(dfbf1,2));
maskI = find(masko);
dfbf(maskI,:) = dfbf1;
dfbf = reshape(dfbf,size(mask,1),size(mask,2),size(dfbf1,2));
dfbf = applyMask(dfbf,mask);
dfbf = applySpatialFilter(dfbf);
minI = min(dfbf(:));
maxI = max(dfbf(:));
quit = 1;
frn = 21;
while quit
    figure(111);clf
    imagesc(dfbf(:,:,frn),[minI maxI]);
    axis equal; axis off;
    titleText = sprintf('%d - %d - %s - %d',trn,frn,info.name,info.props.sampling_frequency);
    title(titleText);
    colorbar;
    ch = getkey;
    if ch == 28 % left arrow
        if frn > 1
            frn = frn - 1;
        end
    end
    if ch == 31 % up arrow
        if (frn-50) > 1
            frn = frn - 50;
        end
    end
    if ch == 29 % right arrow
        if frn < size(dfbf,3)
            frn = frn + 1;
        end
    end
    if ch == 30 % right arrow
        if (frn+50) < size(dfbf,3)
            frn = frn + 50;
        end
    end
    if ch == 27 % escape
        mainQuit = 0;
        break;
    end
    if ch == 113 % q
        mainQuit = 1000;
        break;
    end
    if ch == 118 % v
        sel = generalGUIForSelection(rois.names);
        thisPixels = rois.pixels{sel};
        thisMask = makeMask(size(dfbf(:,:,1)),thisPixels);
        thisdfbyf = mean(getMaskedValues(dfbf,thisMask));
        figure(112);clf;
        plot(thisdfbyf);
%         plot(locdetrend(thisdfbyf,150,[0.2 0.05]));
    end
    if ch == 97 % a
        figure(113);clf
        mdf = mean(dfbf,3);
        imagesc(mdf,[minI maxI]);
        axis equal; axis off;
        titleText = sprintf('%d - %d - %s - %d',trn,frn,info.name,info.props.sampling_frequency);
        title(titleText);
        colorbar;
%         plot(locdetrend(thisdfbyf,150,[0.2 0.05]));
    end
end