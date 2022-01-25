fileName = 'Data_Info5.xlsx';
dPath{1} = '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data';
fileName = fullfile(dPath{1},fileName);
[num,txt,raw] = xlsread(fileName,6,'A1:P72');

%%

%%
% T = getTable(raw,{'171034'},'recording',''); 

IDs = getSelectedCol(raw,'A');
ids = getIDs(raw,'A');
uids = unique(ids');
RecordingDate = getSelectedCol(raw,'B');

tbl = [];
for ii = 1:length(uids)
    if isnan(uids(ii))
        continue;
    end
    inds = ids == uids(ii);
    tbl = [tbl;raw(inds,:)];
    tbl = [tbl;cell(2,size(raw,2))];
end

tbl_xl = cell2table(tbl);
writetable(tbl_xl,'data_RSC_2.xls');


