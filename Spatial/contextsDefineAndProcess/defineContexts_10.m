

function good = defineContexts_10(ei,ow)

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



% contextIDs = [7 8 7 9 7 6 7 10 1 2]'; perhaps for protocol 16
contextIDs = [3 5 6 7]'; % for protocol 10


table(contextIDs,allContexts.names(contextIDs))

% contextTrials = {
%     [1:10];
%     [15:24];
%     [25:35];
%     [36:46];
%     };

protocol = 'Protocol 10';
belt_length = 150; % in cm

contextTrials = {
    [1:10];
    [12:21];
    {[23:32];[22:32];[23:32];[23:32];[23:32];[23:32];};
    {[34:43];[33:43];[34:43];[34:43];[34:43];[34:43];};
    };

table(contextIDs,allContexts.names(contextIDs),contextTrials)

% what do you want to look at in the trials
typesOfMarkers = allContexts.typesOfMarkers;
table((1:length(typesOfMarkers))',typesOfMarkers)

contextTypeOfMarkers = {
    [1 2 9];%13 14 17 18];
    [1 2 9];%13 14 17 18];
    [1 2 9];%13 14 17 18];
    [1 2 9];%13 14 17 18];
    };

% contextTypeOfMarkers = {
%     [1 2 13 14 17 18];
%     [1 2 13 14 17 18];
%     [1 2 13 14 17 18];
%     [1 2 13 14 17 18];
%     };

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