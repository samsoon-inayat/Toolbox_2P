function [dS,T] = get_exp_info_from_folder(thisFolder,pFolder,animal_list,date_list)

if ~exist('animal_list','var')
    animal_list = [];
    date_list = [];
end

[dS,T] = get_exp_info_from_rfolder(thisFolder,pFolder,animal_list,date_list);


function [dS,T] = get_exp_info_from_rfolder(thisFolder,pFolder,animal_list,date_list)
ind = 1;
% for ii = 1:length(data_folder)
%     thisFolder = f.mainDataFolderList{ii};
    files = dir(thisFolder);
    for jj = 1:length(files)
        thisName = files(jj).name;
        if ~files(jj).isdir | strcmp(thisName,'.') | strcmp(thisName,'..')
            continue;
        end
        if ~isempty(str2num(thisName(1)))
            thisFolderD = fullfile(thisFolder,thisName);
            files_listD = dir(sprintf('%s\\**\\Image_0001_0001.raw',thisFolderD));
            filesD = dir(thisFolderD);
            if isempty(files_listD)
                files_listD = dir(sprintf('%s\\**\\Image_001_001.raw',thisFolderD));
                if isempty(files_listD)
                    continue;
                end
            end
            if ~any(strcmp(animal_list,thisName))
                continue;
            end
            if strcmp(thisName,'183762')
                n = 0;
            end
            animal(ind).ID = thisName;
            animalList{ind} = thisName;
            animalListNum(ind) = str2double(thisName);
            animal(ind).folder = thisFolderD;
%             for iii = 1:length(files_listD)
%                 thisFileFolder = files_listD(iii).folder;
%                 strfind(thisFileFolder,thisFolderD)
%             end
            indD = 1;
            for kk = 1:length(filesD)
                thisNameD = filesD(kk).name;
                if ~filesD(kk).isdir | strcmp(thisNameD,'.') | strcmp(thisNameD,'..')
                    continue;
                end
                if ~any(strcmp(date_list,thisNameD))
                    continue;
                end
                if ~isempty(str2num(thisNameD(1)))
                    thisFolderR = fullfile(thisFolderD,thisNameD);
                    files_list = dir(sprintf('%s\\**\\Image_0001_0001.raw',thisFolderR));
                    if isempty(files_list)
                        files_list = dir(sprintf('%s\\**\\Image_001_001.raw',thisFolderR));
                        if isempty(files_list)
                            continue;
                        end
                    end
                    posSlash = strfind(files_list.folder,'\');
                    animal(ind).date(indD).folder = thisFolderR;
                    animal(ind).dateList{indD} = thisNameD;
                    animal(ind).date(indD).files = files_list;
                    dS.exp_list_animal{ind,indD} = thisName;
                    dS.exp_list_date{ind,indD} = thisNameD;
                    dS.exp_list_expDir{ind,indD} = files_list.folder((posSlash(end)+1):end);
                    indD = indD + 1;
                end
            end
            ind = ind+1;
        end
    end
% end
dS.animal = animal;
dS.animalList = animalList;
dS.animalListNum = animalListNum;
dS.data_folder = thisFolder;
dS.s2p_processed_data_folder = pFolder{1};
dS.matlab_processed_data_folder = pFolder{2};
if length(pFolder) == 3
    dS = rmfield(dS,'matlab_processed_data_folder');
    dS.matlab_processed_data_folder_time = pFolder{2};
    dS.matlab_processed_data_folder_dist = pFolder{3};
end

T = get_table(dS);

function T = get_table(dS)
ind = 1;
for rr = 1:size(dS.exp_list_animal,1)
    for cc = 1:size(dS.exp_list_animal,2)
        if isempty(dS.exp_list_animal{rr,cc})
            continue;
        end
        T{ind,1} = dS.exp_list_animal{rr,cc};
        T{ind,2} = dS.exp_list_date{rr,cc};
        T{ind,3} = dS.exp_list_expDir{rr,cc};
        T{ind,6} = fullfile(dS.data_folder,dS.exp_list_animal{rr,cc},dS.exp_list_date{rr,cc},dS.exp_list_expDir{rr,cc});
        T{ind,7} = fullfile(dS.s2p_processed_data_folder,dS.exp_list_animal{rr,cc},dS.exp_list_date{rr,cc},dS.exp_list_expDir{rr,cc});
        if isfield(dS,'matlab_processed_data_folder')
            T{ind,8} = fullfile(dS.matlab_processed_data_folder,dS.exp_list_animal{rr,cc},dS.exp_list_date{rr,cc},dS.exp_list_expDir{rr,cc});
        else
            T{ind,8} = fullfile(dS.matlab_processed_data_folder_time,dS.exp_list_animal{rr,cc},dS.exp_list_date{rr,cc},dS.exp_list_expDir{rr,cc});
            T{ind,9} = fullfile(dS.matlab_processed_data_folder_dist,dS.exp_list_animal{rr,cc},dS.exp_list_date{rr,cc},dS.exp_list_expDir{rr,cc});
        end
        ind = ind + 1;
    end
end

T = array2table(T);

n = 0;