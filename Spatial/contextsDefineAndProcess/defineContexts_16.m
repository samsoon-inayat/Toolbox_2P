function good = defineContexts_16(ei,ow)

if ~exist('ei','var')
    ei = evalin('base','ei([1])');
end

if ~exist('ow','var')
    ow = 0;
end

thispFolder = ei{1}.plane{1}.folder;
fileName = fullfile(thispFolder,'contexts.mat');
if exist(fileName,'file') && ow == 0
    disp('Found existing file ... if you want to run edit code to disable this check');
    good = 0;
    return;
end

allContexts = contextDefinitions;



contextIDs = [7 8 7 9 7 6 7 10 1 2]'; % for protocol 16
table(contextIDs,allContexts.names(contextIDs))

protocol = 'Protocol 16';
belt_length = 150; % in cm


% contextTrials = {
%     [1:11];
%     [12:22];
%     [23:33];
%     [34:43];
%     [44:54];
%     [55:65];
%     [66:76];
%     [81:90];
%     [91:100];
%     [77:86];
%     };
% 
contextTrials = {
    [1:11];
    [12:22];
    [23:33];
    [34:44];
    [45:55];
    [57:67];
    [68:78];
    [81:90];
    [91:100];
    [79:88];
    };


table(contextIDs,allContexts.names(contextIDs),contextTrials)

% what do you want to look at in the trials
typesOfMarkers = allContexts.typesOfMarkers;
table((1:length(typesOfMarkers))',typesOfMarkers)

contextTypeOfMarkers = {
    [2];
    [2];
    [2];
    [2];
    [2];
    [2];
    [2];
    [12];
    [11];
    [1];
    };

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

for ii = 1:length(ei)
    for pp = 1:length(ei{ii}.plane)
        thispFolder = ei{ii}.plane{pp}.folder;
        disp(sprintf('Saving contexts %s',thispFolder));
        fileName = fullfile(thispFolder,'contexts.mat');
        save(fileName,'contexts','belt_length','protocol');
    end
end


good = 1;