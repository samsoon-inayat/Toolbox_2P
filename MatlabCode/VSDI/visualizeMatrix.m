function visualizeMatrix (dfbf1,mask,masko)
% dfbf = zeros(size(mask,1),size(mask,2),size(dfbf1,2));
% dfbf = reshape(dfbf,size(mask,1)*size(mask,2),size(dfbf1,2));
% maskI = find(masko);
% dfbf(maskI,:) = dfbf1;
% dfbf = reshape(dfbf,size(mask,1),size(mask,2),size(dfbf1,2));
% dfbf = applyMask(dfbf,mask);
dfbf = convertToImgSeq(dfbf1,masko);
dfbf = applySpatialFilter(dfbf);
dfbf = applyMask(dfbf,mask);
minI = min(dfbf(:));
maxI = max(dfbf(:));
quit = 1;
frn = 45;
while quit
    figure(111);clf
    imagesc(dfbf(:,:,22),[minI maxI]);
    axis equal; axis off;
%     titleText = sprintf('%d - %d - %s - %s',trn,frn,info.name,info.stimulus.stimulus_type);
%     title(titleText);
    colorbar;
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
    if ch == 27
        break;
    end
end