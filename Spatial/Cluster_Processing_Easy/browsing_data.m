% browsing_data
fileName = 'Data_Info.xlsx';
dPath{1} = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Data\RSEG_PSEG_more_data';
fileName = fullfile(dPath{1},fileName);
[num,txt,raw] = xlsread(fileName,1,'A1:T56');

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
writetable(tbl_xl,'data_basic.xls');

