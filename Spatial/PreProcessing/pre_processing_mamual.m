rFolder = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Data\RSC_HPC\003351\2020-08-09\1'
pd_path = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Processed_Data\RSC_HPC\003351\2020-08-09\1'

eip = thorGetExperimentInfo(rFolder);

fileName = 'bidishift.mat';
fn = fullfile(pd_path,fileName);
bidishift = get_bidi_shift(eip,1,100);
save(fn,'bidishift');
