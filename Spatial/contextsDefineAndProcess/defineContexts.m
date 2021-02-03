function contexts = defineContexts(ei,ow)

% if ~exist('ei','var')
%     ei = evalin('base','ei([1])');
% end
% 
% if ~exist('ow','var')
%     ow = 0;
% end
% 
% 
% thispFolder = ei{1}.plane{1}.folder;
% fileName = fullfile(thispFolder,'contexts.mat');
% if exist(fileName,'file') && ow == 0
%     disp('Found existing file ... if you want to run edit code to disable this check');
%     good = 0;
%     return;
% end

allContexts = contextDefinitions;

cvs = get_define_contexts_file(ei);
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