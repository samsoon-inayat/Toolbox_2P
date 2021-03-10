function training_data = behaviorProcessor

% [f,cName,D] = getFolders;
temp = load('mohajerani_nas_drive_data_list.mat');
f = temp.f;
cName = temp.cName;
D = temp.D;

T15_1 = load('T_15_1_Thy1.mat');
T15 = load('T15.mat');
% T16 = load('T16.mat');
% T10 = load('T.mat');

T = T15_1.T;

T = [T;T15.T];
% T = [T;T16.T];
% T = [T;T10.T];

animalIDs = [];
for ii = 1:size(T,1)
    animalIDs = [animalIDs;cell2mat(T{ii,1})];
end

uaids = unique(animalIDs);
uaids = [183633;183761;183745;183628;183762];



fileName = 'Data_Info5.xlsx';
dPath{1} = '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data';
dPath{2} = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\Data';
pdPath = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\S_Drive\Processed_Data';
fileName = fullfile(dPath{1},fileName);

[num,txt,raw] = xlsread(fileName,10,'A1:N333');
Tf = [];
emptyRow = {NaN,'','','',''};
for ii = 1:length(uaids)
    ii
    animalID = uaids(ii);
    TT = getTable(raw,animalID,'','');
    Tf = [Tf;TT];
    Tf = [Tf;emptyRow];
end
n = 0;
% writetable(Tf,'allData.xls');
uaids1 = uaids;%([1 3:length(uaids)]);
% uaids1(uaids1 == 183745) = [];
Tr = [];
for ii = 1:length(uaids1)
    if ii == 2
        n = 0;
    end
    animalID = uaids1(ii);
    subT = T(animalID == animalIDs,:);
    datesHere = convertToDate(subT{:,2});
    datesHere = sort(datesHere);
    firstRecordingDate = datesHere(1);
    recordingFolders = subT{:,6};
    pdFolders = subT{:,7};
    TT = getTable(raw,animalID,'training',[]);
    trainingDates = convertToDate(TT{:,2});
    inds = find(trainingDates < firstRecordingDate);
    inds = sort(inds);
    N(ii) = length(inds);
    folderName = recordingFolders{1};
    pos = strfind(folderName,'\');
    folderName = folderName(1:(pos(end-1)-1));
    pfolderName = pdFolders{1};
    pos = strfind(pfolderName,'\');
    pfolderName = pfolderName(1:(pos(end-1)-1));
    for jj = 1:length(inds)
        thisTrainingDate = trainingDates(inds(jj));
        files = dir(sprintf('%s\\**\\*.abf',folderName));
        found = 0;
        for kk = 1:length(files)
            thisFileName = files(kk).name;
            yearstr = str2double(thisFileName(1:4));
            monthstr = str2double(thisFileName(6:7));
            daystr = str2double(thisFileName(9:10));
            if yearstr == year(thisTrainingDate) && monthstr == month(thisTrainingDate) && daystr == day(thisTrainingDate)
                found = 1;
                break;
            end
        end
        if found == 0
            continue;
        end
        fileName = fullfile(files(kk).folder,files(kk).name);
        disp(fileName);
        pos = strfind(thisFileName,'.');
        pfileName = fullfile(folderName,sprintf('%s.mat',thisFileName(1:(pos-1))));
        if exist(pfileName,'file')
            temp = load(pfileName);
            b = calcBehav(temp.b);
            bs{ii,jj} = b;
            try
                td(ii,jj) = day(thisTrainingDate - trainingDates(1));
            catch
                td(ii,jj) = (minutes(thisTrainingDate - trainingDates(1))/60)/24;
            end
            t_date{ii,jj} = thisTrainingDate;
            continue;
        else
            try
                b = process_abf_here(fileName);
                save(pfileName,'b');
                b = calcBehav(b);
                bs{ii,jj} = b;
                td(ii,jj) = day(thisTrainingDate - trainingDates(1));
                t_date{ii,jj} = thisTrainingDate;
            catch
                bs{ii,jj} = [];
                td(ii,jj) = NaN;
                t_date{ii,jj} = thisTrainingDate;
                continue;
            end
        end
    end
end
training_data.bs = bs;
training_data.training_dates = t_date;
training_data.training_days = td;
training_data.animalIDs = uaids1;

n = 0;
function b = process_abf_here(fileName)
% db channel names defined for abf files
db.channel{1} = 'frames';
db.channel{2} = 'ch_a';
db.channel{3} = 'ch_b';
db.channel{4} = 'photo_sensor';
db.channel{5} = 'air_puff';
[d,si] = abf2load(fileName);
channel = identify_abf_channels(d,si);
if ~isequal(db.channel,channel)
    db.channel = channel;
    disp('Channels inconsistency was found in abf');
end

b = process_abf_data(d,si,db);