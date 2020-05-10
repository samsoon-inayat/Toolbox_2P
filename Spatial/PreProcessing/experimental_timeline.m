function bs = experimental_timeline

% [f,cName,D] = getFolders;
temp = load('mohajerani_nas_drive_data_list.mat');
f = temp.f;
cName = temp.cName;
D = temp.D;

T15_1 = load('T_15_1_Thy1.mat');
T15 = load('T15.mat');
T16 = load('T16.mat');
T10 = load('T.mat');

T = T15_1.T;

T = [T;T15.T];
T = [T;T16.T];
T = [T;T10.T];

animalIDs = [];
for ii = 1:size(T,1)
    animalIDs = [animalIDs;cell2mat(T{ii,1})];
end

uaids = unique(animalIDs);



fileName = 'Data_Info5.xlsx';
dPath{1} = '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data';
dPath{2} = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\Data';
pdPath = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\S_Drive\Processed_Data';
fileName = fullfile(dPath{1},fileName);

[num,txt,raw] = xlsread(fileName,10,'A1:N303');
Tf = [];
emptyRow = {'','','','',''};
for ii = 1:length(uaids)
    animalID = uaids(ii);
    TT = getTable(raw,animalID,'',''); 
    Tf = [Tf;TT];
    Tf = [Tf;emptyRow];
end
n = 0;
writetable(Tf,'allData.xls');

% T15_1 = load('T_15_1_Thy1.mat');
% T15 = load('T15.mat');
% T16 = load('T16.mat');
% T10 = load('T.mat');
% 
% T = T15_1.T;
% 
% T = [T;T15.T];
% T = [T;T16.T];
% T = [T;T10.T];
% 
% for ii = 1:size(T,1)
%     thisDate = T{ii,2};
%     pos = strfind(thisDate,'/')
%     if isempty(pos{1})
%         dateC(ii,1) = datetime(thisDate,'Format','yyyy-MM-dd');
%     else
%         dateC(ii,1) = datetime(thisDate,'Format','MM/dd/yyyy');
%     end
% end
% 
% [I,inds] = sort(dateC);
% 
% Tm = T(inds,:);
% for ii = 1:size(T,1)
%     Tm{ii,2} = cellstr(I(ii));
% end
% 
% animalIDs = [];
% for ii = 1:size(T,1)
%     animalIDs = [animalIDs;cell2mat(Tm{ii,1})];
% end
% 
% uaids = unique(animalIDs);
% 
% for ii = 1:length(uaids)
%     inds = find(animalIDs == uaids(ii));
%     subT = Tm{inds,1:4};
%     animals{ii} = subT;
% end
% n = 0;