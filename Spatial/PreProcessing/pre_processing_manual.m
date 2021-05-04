% '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data\183633\2019-06-04\1_001'
% '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data\183761\2019-06-06\1_002'
% '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data\183745\2019-06-07\1_001'
% '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data\183628\2019-06-11\1_001'
% '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data\183762\2019-06-11\1_001'
rFolder_prefix = '\\mohajerani-nas.uleth.ca\storage\homes\samsoon.inayat\Data';
pdFolder_prefix = 'Z:\homes\brendan.mcallister\2P\Processed_Data_1';
folder_postfix = '\183633\2019-06-04\1_001';
rFolder = [rFolder_prefix folder_postfix];
pd_path = [pdFolder_prefix folder_postfix];

eip = thorGetExperimentInfo(rFolder);

if ~exist(pd_path,'dir')
    mkdir(pd_path);
end

fileName = 'bidishift.mat';
fn = fullfile(pd_path,fileName);
bidishift = get_bidi_shift(eip,1,100);
save(fn,'bidishift');

ii = 1;
% pd_path = d15{ii}.