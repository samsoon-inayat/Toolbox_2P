function recordingFolder = getRecordingFolder(T,D)
nd = [];
for ii = 1:size(T,1)
    ii
    if ii == 15
        n = 0;
    end
    animal_id = (cell2mat(T{ii,1}));
    exp_date = datestr(cell2mat(T{ii,2}));
    ind = D.animalListNum == animal_id;
    if sum(ind) > 1
        tempD = D.animal(ind);
        for dd = 1:length(tempD)
            tempD_dateList = datestr(tempD(dd).dateList);
            ind_dateList = strcmp(cellstr(tempD_dateList),exp_date);
            if sum(ind_dateList) == 0
                continue;
            elseif sum(ind_dateList) == 1
                temp_ind = find(ind);
                ind = temp_ind(dd);
                break;
            else
                error;
            end
        end
        if length(ind) > 1
            error;
        end
    end
    animal_data = D.animal(ind);
    date_list = datestr(animal_data.dateList);
    ind = strcmp(cellstr(date_list),exp_date);
    if sum(ind) == 0 || isempty(ind)
        nd = [nd ii];
        continue;
    end
    animal_date_data = animal_data.date(ind);
    files = animal_date_data.files;
    if isempty(files)
        n = 0;
    end
    for ff = 1:length(files)
        sizeOfFile(ff) = files.bytes;
    end
    ind = sizeOfFile == max(sizeOfFile);
    recordingFolder{ii,1} = files(ind).folder;
end