function edit_define_contexts_file(ei,bp)
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
