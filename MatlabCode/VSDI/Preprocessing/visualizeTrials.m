function visualizeTrials (ags,pags,stim_folders_info,frames_info,indices,maskFactor)
indx = indices(1);
gg = stim_folders_info.list(indx,2);
aa = stim_folders_info.list(indx,3);
rr = stim_folders_info.list(indx,4);
% thisRecording = ags(gg).animals{aa}.eRecordings{rr};
pThisRecording = pags(gg).animals{aa}.eRecordings{rr};

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
    fileName = makeName('lo6_df.mat',pThisRecording.root_folder);
    temp = load(fileName);
    temp.dfbyfo = temp.df;
end

% temp.dfbyfo = applyTemporalFilter(temp.df);
% temp.dfbyfo = applySpatialFilter(temp.dfbyfo);
code = visualizeImgSeq(temp.dfbyfo,mask,masko,frames_info,pThisRecording,0,rois,indx+20);


function mainQuit = visualizeImgSeq (dfbf1,mask,masko,frames_info,info,trn,rois,indx)
% dfbf = zeros(size(mask,1),size(mask,2),size(dfbf1,2));
% dfbf = reshape(dfbf,size(mask,1)*size(mask,2),size(dfbf1,2));
% maskI = find(masko);
% dfbf(maskI,:) = dfbf1;
% dfbf = reshape(dfbf,size(mask,1),size(mask,2),size(dfbf1,2));
% dfbf = applyMask(dfbf,mask);
stimulusFrame = find((frames_info(4):frames_info(5)) == frames_info(1))+3;
dfbf = convertToImgSeq(dfbf1,masko);
dfbf = applySpatialFilter(dfbf);
dfbf = applyMask(dfbf,mask);
minI = min(dfbf(:));
maxI = max(dfbf(:));
quit = 1;
frn = stimulusFrame;
while quit
    figure(indx);
    imagesc(dfbf(:,:,frn),[minI maxI]);
    axis equal; axis off;
    titleText = sprintf('%d - %d - %s - %s',trn,frn,info.name,info.stimulus.stimulus_type);
    title(titleText);
    colorbar;
    if exist('RHR','var')
         RHR = thisPixels;
         rectangle('Position',[RHR(1,1),RHR(1,2),5,5],'LineWidth',1,'EdgeColor','k');
    end
    ch = getkey;
    if ch == 28 % left arrow
        if frn > 1
            frn = frn - 1;
        end
    end
    if ch == 31 % up arrow
        if (frn-5) > 1
            frn = frn - 5;
        end
    end
    if ch == 29 % right arrow
        if frn < size(dfbf,3)
            frn = frn + 1;
        end
    end
    if ch == 30 % right arrow
        if (frn+5) < size(dfbf,3)
            frn = frn + 5;
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
        figure(indx+1);clf;
        plot(thisdfbyf);hold on;
        plot(zeros(size(thisdfbyf)),'r');
        figure(indx);
        RHR = thisPixels;
        rectangle('Position',[RHR(1,1),RHR(1,2),5,5],'LineWidth',1,'EdgeColor','r');
        n = 0;
%         plot(locdetrend(thisdfbyf,150,[0.2 0.05]));
    end
    if ch == 97 % a
        figure(indx+2);clf
        mdf = mean(dfbf(:,:,stimulusFrame:(stimulusFrame+5)),3);
        imagesc(mdf,[minI maxI]);
        axis equal; axis off;
        titleText = sprintf('%d - %d - %s - %s',trn,frn,info.name,info.stimulus.stimulus_type);
        title(titleText);
        colorbar;
%         plot(locdetrend(thisdfbyf,150,[0.2 0.05]));
    end
end