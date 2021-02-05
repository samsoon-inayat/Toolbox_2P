function edit_define_contexts_file(ei,dcfilename,bp)
if ~exist('dcfilename','var')
	dcfilename = 'define_contexts.m';
end

if ~exist('bp','var')
    bp = 0;
end
if bp
    behaviorPlot({ei});
end
fileName = fullfile(ei.recordingFolder,dcfilename);

if ~exist(fileName,'File')
    if strcmp(dcfilename,'define_contexts.m')
        copyfile('define_contexts.m',fileName);
    else
        fileName1 = fullfile(ei.recordingFolder,'define_contexts.m');
        copyfile(fileName1,fileName);
    end
end

edit(fileName);
