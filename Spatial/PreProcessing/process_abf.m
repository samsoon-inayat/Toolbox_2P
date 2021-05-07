function process_abf(T,owr)

for ii = 1:size(T,1)
    f.recordingFolder = cell2mat(T{ii,6});
    disp(f.recordingFolder);
    if ~isempty(strfind(f.recordingFolder,'Missing'))
        continue;
    end
    db = get_db(f);
    ei = thorGetExperimentInfo(f.recordingFolder);
    ei.recordingFolder = f.recordingFolder;
    for jj = 1:length(db)
        ei.db = db(jj);
        abf2behavior_1(ei,cell2mat(T{ii,7}),'overwrite_behavior',owr);
    end
end