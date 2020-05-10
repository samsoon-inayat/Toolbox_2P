% function experiments_XL_record
add_to_path

% clear all
% clc
fileName = 'Data_Info3.xlsx';
dPath{1} = '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data';
dPath{2} = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\Data';
pdPath = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\S_Drive\pySuite2p_Processed_Data';
fileName = fullfile(dPath{1},fileName);

readAgain = 0;
% read excel file if not already read
if ~exist('raw','var') | readAgain
    [num,txt,raw] = xlsread(fileName,10,'A1:N293');
    [f,cName] = getFolders;
end

AD_Thy1_animals = [183224;183227;183228;183329];
Thy1_animals = [173062;173511;173198;174374;173706;183633;183761;183745;183628;183762];
db = 0; bidi = 0; s2p = 0;
for ii = 1:length(Thy1_animals)
    sel_id = Thy1_animals(ii);
    preProcess(f,raw,sel_id,[db bidi s2p]);
end

function preProcess(f,raw,sel_id,opts)
dPath{1} = '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data';
dPath{2} = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\Data';
pdPath = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\S_Drive\Processed_Data';
% pdPath = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\S_Drive\pySuite2p_Processed_Data';




% columns in excel file
IDs = getSelectedCol(raw,'A');
RecordingDate = getSelectedCol(raw,'B');
TrainingOrRecording = getSelectedCol(raw,'E');
Protocol = getSelectedCol(raw,'I');
PerformanceNotes = getSelectedCol(raw,'N');

selection = [];



% specify min and max age
ids = getIDs(raw,'A');
selection{1} = find(ids == sel_id);
IndexC = strfind(lower(TrainingOrRecording), 'recording');
selection{2} = find(not(cellfun('isempty', IndexC)));
IndexC = strfind(lower(Protocol), 'protocol 10');
selection{3} = find(not(cellfun('isempty', IndexC)));


% common = intersect(selection{1},selection{2});
common = selection{1};
if length(selection) > 1
    for ii = 2:length(selection)
        common = intersect(common,selection{ii});
    end
end
% 

% varNames = {'ID','RecordingDate','TrainingOrRecording','Protocol','PerformanceNotes'};
% T = table(IDs(common),RecordingDate(common),TrainingOrRecording(common),Protocol(common),PerformanceNotes(common),'VariableNames',varNames)

dataFolder = [];
bad = [];
for ii = 1:length(common)
    if ~isempty(strfind(lower(PerformanceNotes{common(ii)}),'euthanize')) || ~isempty(strfind(lower(PerformanceNotes{common(ii)}),'bad'))
        bad = [bad ii];
        continue;
    end
    thisDate = RecordingDate(common(ii));
    thisDateF = datestr(thisDate);
    dateFolder = datestr(thisDateF,'yyyy-mm-dd');
    for jj = 1:length(dPath)
        thisdPath = dPath{jj};
        folderName = fullfile(thisdPath,num2str(sel_id),dateFolder);
        if exist(folderName,'dir')
            dataFolder{ii,1} = folderName;
        else
            n = 0;
        end
    end
end
try
    common(bad) = [];
catch
end
try
    dataFolder(bad) = [];
catch
end
% varNames = {'ID','RecordingDate','TrainingOrRecording','Protocol','dataFolder'};
% T = table(IDs(common),RecordingDate(common),TrainingOrRecording(common),Protocol(common),dataFolder,'VariableNames',varNames)


pti = 1;

    for ii = 1:length(common)
        ci = common(ii);
        thisDataFolder = dataFolder{ii};
        if isempty(thisDataFolder)
            continue;
        end 
        thisDate = RecordingDate(common(ii));
        thisDateF = datestr(thisDate);
        dateFolder = datestr(thisDateF,'yyyy-mm-dd');
        files = dir(thisDataFolder);
        dbFileStr = [];itf = 1;
        folderList = [];ind = 1;
        for jj = 1:length(files)
            thisName = files(jj).name;
            if ~files(jj).isdir | strcmp(thisName,'.') | strcmp(thisName,'..')
                continue;
            end
            rFolder = fullfile(thisDataFolder,thisName);
            try
                eip = thorGetExperimentInfo(rFolder);
            catch
                continue;
            end
            fileName = fullfile(rFolder,'Image_0001_0001.raw');
            if ~exist(fileName,'file')
                continue;
            end
            nP = getNumberOfPlanes(eip);
            for mm = 1:nP
                pID{pti,1} = IDs{ci};
                pD{pti,1} = RecordingDate{ci};
                pTR{pti,1} = TrainingOrRecording{ci};
                pPr{pti,1} = Protocol{ci};
                pPlane{pti,1} = mm;
                p_planePath = fullfile(pdPath,num2str(sel_id),dateFolder,num2str(mm));
                fName = sprintf('bidishift.mat');
                fileName = fullfile(p_planePath,fName);
                if ~exist(fileName,'file')
                    pBiDi{pti,1} = 'No';
                    pS2p{pti,1} = 'No';
                    pNM{pti,1} = 'No';
                    continue;
                else
                    pBiDi{pti,1} = 'Yes';
                end
                p_planePath = fullfile(pdPath,num2str(sel_id),dateFolder,num2str(mm));
                fName = sprintf('F_%s_%s_plane1.mat',num2str(sel_id),datestr(thisDateF,'yyyy-mm-dd'));
                fileName = fullfile(p_planePath,fName);
                if ~exist(fileName,'file')
                    pS2p{pti,1} = 'No';
                    pNM{pti,1} = 'No';
                    continue;
                else
                    pS2p{pti,1} = 'Yes';
                end
                
                p_planePath = fullfile(pdPath,num2str(sel_id),dateFolder,num2str(mm));
                fName = sprintf('F_%s_%s_plane1_proc.mat',num2str(sel_id),datestr(thisDateF,'yyyy-mm-dd'));
                fileName = fullfile(p_planePath,fName);
                if ~exist(fileName,'file')
                    pNM{pti,1} = 'No';
                else
                    pNM{pti,1} = 'Yes';
                end
                pti = pti + 1;
            end
        end
    end
    if exist('pID','var')
    varNames = {'ID','RecordingDate','Plane','TrainingOrRecording','Protocol','BiDi','S2P','New_Main'};
    T = table(pID,pD,pPlane,pTR,pPr,pBiDi,pS2p,pNM,'VariableNames',varNames)
    end

end

function ids = getIDs(raw,colA,selRows)
    Alphabets = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    col = strfind(Alphabets,colA);
    if ~exist('selRows','var')
        selRows = 1:size(raw,1);
    end
    for ii = 1:length(selRows)
        try
            ids(ii) = (raw{selRows(ii),col});
        catch
            ids(ii) = NaN;
        end
    end
end

function rows = getSelectedCol(raw,colA,selRows)
    Alphabets = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    col = strfind(Alphabets,colA);
    if ~exist('selRows','var')
        selRows = 1:size(raw,1);
    end
    for ii = 1:length(selRows)
        if isnan(raw{selRows(ii),col})
%             if strcmp(colA,'A')
%                 rows(ii,1) = ' ';
%             else
                rows{ii,1} = ' ';
%             end
        else
%             if strcmp(colA,'A')
%                 rows(ii,1) = raw{selRows(ii),col};
%             else
                rows{ii,1} = raw{selRows(ii),col};
%             end
        end
    end
end

