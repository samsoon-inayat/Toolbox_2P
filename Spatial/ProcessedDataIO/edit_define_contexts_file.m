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
fileName = fullfile(ei.recordingFolder,'define_contexts.m');

if ~exist(fileName,'File')
    copyfile('define_contexts.m',fileName);
end

edit(fileName);
