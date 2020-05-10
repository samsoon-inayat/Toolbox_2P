function [f,compName,varargout] = getFolders

f.mainDataFolder = '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data';
% f.mainDataFolder2 = '\\142.66.96.22\storage2\homes\samsoon.inayat\Data';
f.mainDataFolderList{1} = '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data';
f.mainDataFolderList{2} = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\Data';
f.nsDrive = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\S_Drive';
[~, name] = system('hostname');

if strcmp(strtrim(name),'Sam-WS')
    f.sDrive = 'S:\';
end

if strcmp(strtrim(name),'imaging3-psy16')
    f.sDrive = 'E:\Users\samsoon.inayat\S_Drive';
end

if strcmp(strtrim(name),'imaging4-psy16')
    f.sDrive = 'E:\Users\samsoon.inayat\S_Drive';
end

if strcmp(strtrim(name),'imaging1-psy16')
    f.sDrive = 'E:\Users\samsoon.inayat\S_Drive';
end

f.tifFolder = fullfile(f.sDrive,'Tif_Data');
f.pFolder = fullfile(f.nsDrive,'Processed_Data');
f.py_pFolder = fullfile(f.nsDrive,'pySuite2P_Processed_Data');

compName = strtrim(name);
if nargout == 3
    varargout{1} = getExpInfo(f);
end

function dS = getExpInfo(f)
ind = 1;
for ii = 1:length(f.mainDataFolderList)
    thisFolder = f.mainDataFolderList{ii};
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
                continue;
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
                if ~isempty(str2num(thisNameD(1)))
                    thisFolderR = fullfile(thisFolderD,thisNameD);
                    files_list = dir(sprintf('%s\\**\\Image_0001_0001.raw',thisFolderR));
                    if isempty(files_list)
                        continue;
                    end
                    animal(ind).date(indD).folder = thisFolderR;
                    animal(ind).dateList{indD} = thisNameD;
                    animal(ind).date(indD).files = files_list;
                    indD = indD + 1;
                end
            end
            ind = ind+1;
        end
    end
end
dS.animal = animal;
dS.animalList = animalList;
dS.animalListNum = animalListNum;
% animal = [];
% animalList = [];
% ind = 1;
% % for ii = 1:length(f.pFolder)
%     if iscell(f.pFolder)
%         thisFolder = f.pFolder{1};
%     else
%         thisFolder = f.pFolder;
%     end
%     files = dir(thisFolder);
%     for jj = 1:length(files)
%         thisName = files(jj).name;
%         if ~files(jj).isdir | strcmp(thisName,'.') | strcmp(thisName,'..')
%             continue;
%         end
%         if ~isempty(str2num(thisName(1)))
%             animal(ind).ID = thisName;
%             animalList{ind} = thisName;
%             thisFolderD = fullfile(thisFolder,thisName);
%             animal(ind).folder = thisFolderD;
%             filesD = dir(thisFolderD);
%             indD = 1;
%             for kk = 1:length(filesD)
%                 thisNameD = filesD(kk).name;
%                 if ~filesD(kk).isdir | strcmp(thisNameD,'.') | strcmp(thisNameD,'..')
%                     continue;
%                 end
%                 if ~isempty(str2num(thisNameD(1)))
%                     animal(ind).dateList{indD} = thisNameD;
%                     thisFolderR = fullfile(thisFolderD,thisNameD);
%                     animal(ind).date(indD).folder = thisFolderR;
%                     filesR = dir(thisFolderR);
%                     indR = 1;
%                     for mm = 1:length(filesR)
%                         thisNameR = filesR(mm).name;
%                         if ~filesR(mm).isdir | strcmp(thisNameR,'.') | strcmp(thisNameR,'..')
%                             continue;
%                         end
%                         animal(ind).date(indD).recording(indR).folder = fullfile(thisFolderR,thisNameR);
%                         animal(ind).date(indD).recordingList{indR} = thisNameR;
%                         indR = indR + 1;
%                     end
%                     indD = indD + 1;
%                 end
%             end
%             ind = ind+1;
%         end
%     end
% % end
% dS.p_animal = animal;
% dS.p_animalList = animalList;
