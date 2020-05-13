function visualize_average_dfbyf (animal,fn,stimType,frames)

for an = fn%1:length(animal.data_folders)
    dataFolder = animal.data_folders{an};
    mask = getMask(dataFolder,0.9);
    if isempty(mask)
        display('Mask does not exist');
        break;
    end
    peDataFolder = animal.processed_data_folders{an}
    fileName = makeName('folders_data.mat',peDataFolder);
    temp = load(fileName);
    folders = temp.folders_data; clear temp;
    folders_data = folders;
    for ii = 1:length(folders_data)
        if strcmp(stimType,folders_data(ii).stimulus.stimulus_type)
            break;
        end
    end
    evFolder = makeName(folders(ii).name,peDataFolder);
    fileName = makeName('dfbyfo.mat',evFolder);
    temp = load(fileName);
%             temp.dfbyfo = applyTemporalFilter(temp.dfbyfo);
    temp.dfbyfo = applySpatialFilter(temp.dfbyfo);
    code = visualizeImgSeq(temp.dfbyfo,mask,folders(ii),frames);
end

function mainQuit = visualizeImgSeq (dfbf,mask,info,frames)
dfbf = applyMask(dfbf,mask,'shiftToZero');
minI = min(dfbf(:));
maxI = max(dfbf(:));
quit = 1;
frn = 45;
while quit
    figure(111111);clf
    imagesc(dfbf(:,:,frn),[minI maxI]);
    axis equal; axis off;
    titleText = sprintf('%d - %d - %s - %s',trn,frn,info.name,info.stimulus.stimulus_type);
    title(titleText);
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
    if ch == 112 % escape
        figure(111112);clf
        for ii = 1:size(dfbf,3)
            imagesc(dfbf(:,:,ii),[minI maxI]);
            axis equal; axis off;
            titleText = sprintf('%d - %s - %s',ii,info.name,info.stimulus.stimulus_type);
            title(titleText);
            colorbar;
            pause(0.05);
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
    if ch == 98 % b
        mainQuit = 10;
        break;
    end
end