function markCenterOfResponse (animal,nums,stimType)

for an = nums%1:length(animal.data_folders)
    dataFolder = animal.data_folders{an};
    mask = getMask(dataFolder);
    if isempty(mask)
        display('Mask does not exist');
        break;
    end
    peDataFolder = animal.processed_data_folders{an}
    fileName = makeName('folders_data.mat',peDataFolder);
    temp = load(fileName);
    folders_data = temp.folders_data; clear temp;
    for ii = 1:length(folders_data)
        if ~strcmp(folders_data(ii).stimulus.stimulus_type,stimType)
            continue;
        end
        evFolder = makeName(folders_data(ii).name,peDataFolder);
        fileName = makeName('dfbyfo.mat',evFolder);
        temp = load(fileName);
        temp.dfbyfo = applySpatialFilter(temp.dfbyfo);
        [code,center,frn] = getCenter(temp.dfbyfo,mask,folders_data(ii),'Select Primary',47);
        if code ~= 1000
            [code,center,frn] = getCenter(temp.dfbyfo,mask,folders_data(ii),'Select Secondary',frn);
        end
        if code == 1000
            break;
        end
        folders_data(ii).response.center = center;
    end
    if code == 1000
        break;
    end
    fileName = makeName('folders_data.mat',peDataFolder);
    save(fileName,'folders_data');
end

function [mainQuit,center,frn] = getCenter(dfbf,mask,info,titleT,frns)
center = [0,0];
% dfbf = applyMask(dfbf,mask,'same');
minI = min(dfbf(:)); tminI = minI;
maxI = max(dfbf(:)); tmaxI = maxI;
quit = 1;
frn = frns;
mainQuit = 0;
while quit
    hf = figure(1);clf
    thisFrame = dfbf(:,:,frn);
    
    imagesc(thisFrame,[tminI tmaxI]);hold on;
    vals = getMaskedValues(thisFrame,mask);
    try
        BW = im2bw(thisFrame,mean(vals)+1*std(vals));
        s = regionprops(BW,'centroid');
        centroids = cat(1, s.Centroid);
        contour(BW);
        plot(centroids(:,1),centroids(:,2), 'b*')
    catch
    end
    axis equal; %axis off;
    titleText = sprintf('%d - %s - %s',frn,info.stimulus.stimulus_type,info.name);
    titleText = {titleT;titleText};
    title(titleText);
    colorbar;
    ch = getkey;
    if ch == 104 % h
        rect = getrect(hf);
        rect = round(rect);
        x1 = rect(1); y1 = rect(2);
        x2 = x1 + rect(3); y2 = y1 + rect(4);
        vals = thisFrame(y1:y2,x1:x2);
        tminI = min(vals(:));%mean(vals)-std(vals);
        tMaxI = max(vals(:));%mean(vals)+std(vals);
    end
    if ch == 106 % j
        tminI = minI; tmaxI = maxI;
    end
    if ch == 99 % c
        rect = getrect(hf);
        x1 = rect(1); y1 = rect(2);
        x2 = x1 + rect(3); y2 = y1 + rect(4);
        iis = [];
        for ii = 1:size(centroids,1)
            if centroids(ii,1) > x1 && centroids(ii,1) < x2
                iis = [iis ii];
            end
        end
        iiss = [];
        for ii = 1:length(iis)
            if centroids(iis(ii),2) > y1 && centroids(iis(ii),2) < y2
                iiss = [iiss iis(ii)];
            end
        end
        if length(iiss) == 1
            center = round([centroids(iiss,2),centroids(iiss,1)]);
        else
            mainQuit = 11;
        end
        break;
    end
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
        break;
    end
    if ch == 113 % q
        mainQuit = 1000;
        break;
    end
    if ch == 98 % b
        mainQuit = 10;
        break;
    end
end