function out = get_define_contexts_file(ei)

fileName = fullfile(ei.recordingFolder,'define_contexts.m');

if ~exist(fileName,'File')
    disp('file does not exist')
end

out = get_mfile_vars(fileName);

