sourceFile = '\\mohajerani-nas.uleth.ca\storage\homes\brendan.mcallister\2P\Data\RSEG_PSEG_more_data\3328\2021-05-01\1\metaD.xml';
for ii = 1:size(T,1)
    thisDate = cell2mat(T{ii,2});
    if strcmp(thisDate,'2021-04-30')
        continue;
    end
    rFolder = cell2mat(T{ii,6});
    destFile = fullfile(rFolder,'metaD.xml');
    if strcmp(sourceFile,destFile)
        continue;
    end
    copyfile(sourceFile,destFile);
end