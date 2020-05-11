% function T = pre_processing
add_to_path
config = get_config;
% clear all
clc
fileName = 'Data_Info4.xlsx';
dPath{1} = '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data';
dPath{2} = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\Data';
pdPath = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\S_Drive\Processed_Data';
fileName = fullfile(dPath{1},fileName);

readAgain = 0;
% read excel file if not already read
if ~exist('raw','var') | readAgain
    [num,txt,raw] = xlsread(fileName,10,'A1:N303');
    [f,cName,D] = getFolders;
    system('C:\Users\samsoon.inayat\AppData\Local\Continuum\anaconda3\envs\suite2p\python.exe T:\GitHub\Toolbox_2P\Spatial\PreProcessing\readGoogleSheet.py','-echo')
    gS = load('temp.mat');
end
% 
% AD_Thy1_animals = [183224;183227;183228;183329];
% animals = AD_Thy1_animals;
% Thy1_animals = [173062;173511;173198;174374;173706;183633;183761;183745;183628;183762];
% Thy1_animals = [183628;183762];
Thy1_animals = [173062;173511;173198;174374;173706;183633;183761;183745;183628;183762];

animals = Thy1_animals;%
db = 0; bidi = 0; s2p = 0; list_r_folders = 0;
list = {}; pList = {};
T = [];
for ii = 1:length(animals)
    sel_id = animals(ii);
    TT = getTable(raw,sel_id,'recording','protocol_15_1'); 
    T = [T;TT];
end
indR = 1; nd = [];recordingFolder = [];
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
T(nd,:) = [];recordingFolder(nd,:) = []; 
if length(recordingFolder) == 1
    T(:,size(T,2)+1) = {recordingFolder};
else
    T(:,size(T,2)+1) = recordingFolder;
end
T.Properties.VariableNames{6} = 'RecordingFolder';
[pdPaths,py_pdPaths] = make_db_and_pdPaths(T);
if length(recordingFolder) == 1
    T(:,(size(T,2)+1)) = {py_pdPaths};
else
    T(:,(size(T,2)+1)) = py_pdPaths;
end
T.Properties.VariableNames{7} = 'Py_ProcessedData';
if length(recordingFolder) == 1
    T(:,(size(T,2)+1)) = {pdPaths(:,1)};
else
    T(:,(size(T,2)+1)) = pdPaths(:,1);
end
T.Properties.VariableNames{8} = 'ProcessedData_Plane1';
if size(pdPaths,2) == 2
    T(:,(size(T,2)+1)) = pdPaths(:,2);
else
    T(:,(size(T,2)+1)) = cell(size(pdPaths,1),1);
end
T.Properties.VariableNames{9} = 'ProcessedData_Plane2';
[pdbidi,py_pdbidi] = check_bidiShift(T);

abf_file_status = check_abf_file(T);

[s2p,nm_s2p,py_s2p] = check_s2p_run(T,f);

list = T{:,6};
fileName = 'recording_list.txt';
fid = fopen(fileName,'w');
for ii = 1:length(list)
    try
        thisLine = list{ii};
    catch
        thisLine = list;
    end
    fprintf(fid, '%s\n',thisLine);
end
fclose(fid);

fileName = 'T_15_1_Thy1.mat';
save(fileName,'T');


function [s2p,nm_s2p,py_s2p] = check_s2p_run(T,f)
    for ii = 1:size(T,1)
        try
            rFolder  = cell2mat(T{ii,6});
        catch
            rFolder  = (T{ii,6});
        end
        eip = thorGetExperimentInfo(rFolder);
        nP = getNumberOfPlanes(eip);
        bs = [];
        for pp = 1:nP
            try
                pd_path = cell2mat(T{ii,7+pp});
            catch
                pd_path = (T{ii,7+pp});
            end
            files = dir(sprintf('%s\\F_*plane1.mat',pd_path));
            if isempty(files)
                f.pd_path = pd_path;
                try
                    f.recordingFolder = cell2mat(T{ii,6});
                catch
                    f.recordingFolder = (T{ii,6});
                end
%                 runSuite2P(f,num2str(cell2mat(T{ii,1})),datestr(cell2mat(T{ii,2}),'yyyy-mm-dd'),'',pp,0,NaN);
                s2p(ii,pp) = 0;
            else
                s2p(ii,pp) = 1;
            end
            files = dir(sprintf('%s\\F_*_proc.mat',pd_path));
            if isempty(files)
                nm_s2p(ii,pp) = 0;
            else
                nm_s2p(ii,pp) = 1;
            end

        end
        try
            pd_path = cell2mat(T{ii,7});
        catch
            pd_path = (T{ii,7});
        end
        files = dir(sprintf('%s\\**\\Fall.mat',pd_path));
        if isempty(files)
            py_s2p(ii,1) = 0;
        else
            py_s2p(ii,1) = 1;
        end
    end
end

function out = check_abf_file(T)
for ii = 1:size(T,1)
    try
        rFolder  = cell2mat(T{ii,6});
    catch
        rFolder  = (T{ii,6});
    end
    files = dir(sprintf('%s\\*.abf',rFolder));
    if isempty(files)
        out(ii,1) = 0;
    else
        out(ii,1) = 1;
    end
end
end

function [pd,py_pd] = check_bidiShift(T)
for ii = 1:size(T,1)
    try
        rFolder  = cell2mat(T{ii,6});
    catch
        rFolder  = (T{ii,6});
    end
    eip = thorGetExperimentInfo(rFolder);
    nP = getNumberOfPlanes(eip);
    fileName = 'bidishift.mat';
    bs = [];
    for pp = 1:nP
        try
            pd_path = cell2mat(T{ii,7+pp});
        catch
            pd_path = (T{ii,7+pp});
        end
        fn = fullfile(pd_path,fileName);
        if exist(fn,'file')
            pd(ii,pp) = 1;
            temp = load(fn);
            bs(pp) = temp.bidishift;
        else
            if isempty(pd_path)
                pd(ii,pp) = NaN;
            else
                bidishift = get_bidi_shift(eip,pp,500);
                if ~exist(pd_path,'dir')
                    mkdir(pd_path);
                end
                save(fn,'bidishift');
                bs(pp) = bidishift;
                pd(ii,pp) = 1;
            end
        end
        
    end
    try
        pd_path = cell2mat(T{ii,7});
    catch
        pd_path = (T{ii,7});
    end
    fn = fullfile(pd_path,fileName);
%     if exist(fn,'file')
%         py_pd(ii,1) = 1;
%     else
        bidishift = bs;
        if ~exist(pd_path,'dir')
            mkdir(pd_path);
        end
        save(fn,'bidishift');
        py_pd(ii,1) = 1;
%     end
end
end

function [pdPaths,py_pdPaths] = make_db_and_pdPaths(T)
pdPath = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\S_Drive\Processed_Data';
py_pdPath = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\S_Drive\pySuite2p_Processed_Data';
for ii = 1:size(T,1)
    try
        rFolder  = cell2mat(T{ii,6});
    catch
        rFolder  = (T{ii,6});
    end
    if ~isempty(strfind(cell2mat(T{ii,4}),'10'))
        copyfile('metaD15.xml',rFolder);
    else
        copyfile('metaD16.xml',rFolder);
    end
    eip = thorGetExperimentInfo(rFolder);
    nP = getNumberOfPlanes(eip);
    dbFile = fullfile(rFolder,'make_db.m');
%     if exist(dbFile,'file')
%         continue;
%     end
    pos = strfind(rFolder,'\Data\');
    py_pdPaths{ii,1} = fullfile(py_pdPath,rFolder(pos+6:end));
    itf = 1;
    dbFileStr{itf,1} = sprintf('i = 0;');itf = itf+1;
    for kk = 1:nP
        if cell2mat(T{ii,1}) == 173511
            n = 0;
        end
        dbFileStr{itf,1} = sprintf('i = i+1;');itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).mouse_name    	= ''%d'';',cell2mat(T{ii,1}));itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).date                = ''%s'';',datestr(cell2mat(T{ii,2}),'yyyy-mm-dd'));itf = itf+1;
%         dbFileStr{itf,1} = sprintf('db(i).expText             = ''%s'';',thisName);itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).expts               = [%d];',kk);itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).selectedPlane    	= %d;',kk);itf = itf+1;
        if ~isempty(strfind(cell2mat(T{ii,4}),'10'))
            dbFileStr{itf,1} = sprintf('db(i).metaD           = ''%s'';','metaD15.xml');itf = itf+1;
        else
            dbFileStr{itf,1} = sprintf('db(i).metaD           = ''%s'';','metaD16.xml');itf = itf+1;
        end
        dbFileStr{itf,1} = sprintf('db(i).nplanes       	= 1;');itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).BiDiPhase     	= 0;');itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).diameter      	= 30;');itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).downsamplespace 	= 1;');itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).downsampletime 	= 1;');itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).downsampletime 	= 1;');itf = itf+1;
        substr = rFolder((pos+6):end);
        posls = strfind(substr,'\');
        pdPaths{ii,kk} = fullfile(pdPath,sprintf('%s\\%d',substr(1:(posls(end)-1)),kk));
    end
    fileID = fopen(dbFile,'w');
    for kk = 1:(itf-1)
        fprintf(fileID,sprintf('%s\n',dbFileStr{kk}));
    end
    fclose(fileID);
end
end

