function visualizeAverage (animal,fn,flagAvg)

for an = fn%1:length(animal.data_folders)
    dataFolder = animal.data_folders{an};
    mask = getMask(dataFolder);
    if isempty(mask)
        display('Mask does not exist');
        break;
    end
    peDataFolder = animal.processed_data_folders{an}
    fileName = makeName('folders_data.mat',peDataFolder);
    temp = load(fileName);
    folders = temp.folders_data; clear temp;
    for ii = 1:length(folders)
        trials = folders(ii).stimulus.actualTrials;
        evFolder = makeName(folders(ii).name,peDataFolder);
        if flagAvg
            fileName = makeName('dfbyfo.mat',evFolder);
            temp = load(fileName);
%             temp.dfbyfo = applyTemporalFilter(temp.dfbyfo);
            temp.dfbyfo = applySpatialFilter(temp.dfbyfo);
            code = visualizeImgSeq(temp.dfbyfo,mask,folders(ii),0);
        else
            badTrials = [];
            for jj = 1:trials
                trialName = sprintf('trial_%.2d',jj);
                fileName = makeName(trialName,evFolder);
                fileName = makeName('dfbyfo.mat',fileName);
                temp = load(fileName);
    %             temp.dfbyfo = applyTemporalFilter(temp.dfbyfo);
                temp.dfbyfo = applySpatialFilter(temp.dfbyfo);
                code = visualizeImgSeq(temp.dfbyfo,mask,folders(ii),jj);
                if code == 1000
                    break;
                end
                if code == 10
                    badTrials = [badTrials jj];
                end
            end
        end
        if code == 1000
            break;
        end
    end
    if code == 1000
        break;
    end
end

function mainQuit = visualizeImgSeq (dfbf,mask,info,trn)
dfbf = applyMask(dfbf,mask);
minI = min(dfbf(:));
maxI = max(dfbf(:));
quit = 1;
frn = 30;
while quit
    figure(1);clf
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