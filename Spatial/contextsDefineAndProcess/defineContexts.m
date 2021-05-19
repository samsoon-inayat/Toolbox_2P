function contexts = defineContexts(ei,dcfilename)


if ~exist('dcfilename','var')
    dcfilename = 'define_contexts.m';
end

allContexts = contextDefinitions;

cvs = get_define_contexts_file(ei,dcfilename);
fields = fieldnames(cvs);
for ii = 1:length(fields)
    cmdTxt = sprintf('%s = cvs.%s;',fields{ii},fields{ii}); eval(cmdTxt);
end
table(contextIDs,allContexts.names(contextIDs))

table(contextIDs,allContexts.names(contextIDs),contextTrials)

% what do you want to look at in the trials
typesOfMarkers = allContexts.typesOfMarkers;
table((1:length(typesOfMarkers))',typesOfMarkers)

for ii = 1:length(contextIDs)
    cid = contextIDs(ii);
    thisContext = allContexts.names{cid};
    thisTrials = contextTrials{ii};
    thisMarkers = typesOfMarkers(contextTypeOfMarkers{ii});
    cmdTxt = sprintf('contexts(ii).name = thisContext;');
    eval(cmdTxt);
    cmdTxt = sprintf('contexts(ii).trials = thisTrials;');
    eval(cmdTxt);
    cmdTxt = sprintf('contexts(ii).stimMarkers = thisMarkers;');
    eval(cmdTxt);
end
n = 0;
% 
% for ii = 1:length(ei)
%     for pp = 1:length(ei{ii}.plane)
%         thispFolder = ei{ii}.plane{pp}.folder;
%         disp(sprintf('Saving contexts %s',thispFolder));
%         fileName = fullfile(thispFolder,'contexts.mat');
%         save(fileName,'contexts','belt_length','protocol');
%     end
% end

good = 1;