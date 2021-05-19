function out = get_define_contexts_file(ei,dcfilename)

if ~exist('dcfilename','var')
    dcfilename = 'define_contexts.m';
end

if isempty(strfind(dcfilename,'.m'))
    dcfilename = [dcfilename '.m'];
end

fileName = fullfile(ei.recordingFolder,dcfilename);

if ~exist(fileName,'File')
    disp('file does not exist')
    out = NaN;
end

out = get_mfile_vars(fileName);
n = 0;

