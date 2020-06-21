function good = defineContexts_protocol15(ei,ow)

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
contextIDs = [1 2 3 4 3 1 2]'; % for protocol 10
table(contextIDs,allContexts.names(contextIDs))

% for first 3 animals before 174374 protocol 10
% contextTrials = {
%     [1:10];
%     [15:24];
%     [25:35];
%     [36:46];
%     };
protocol = 'Protocol 15';
belt_length = 150; % in cm
% for later 6 animals starting 174374 protocol 10
% contextTrials = {
%     [1:10];
%     [12:21];
%     [22:32];
%     [33:43];
%     };
contextTrials = {
    [1:10];
    [1:10];
    [11:20];
    {[21:30];[21:30];[21:30];[11:20];};
    [31:40];
    [21:30];
    [41:50];
    };
table(contextIDs,allContexts.names(contextIDs),contextTrials)
% what do you want to look at in the trials
typesOfMarkers = allContexts.typesOfMarkers;
table((1:length(typesOfMarkers))',typesOfMarkers)

contextTypeOfMarkers = {
    [23];
    [24];
    [1 2 9];
    [1 2 9 23];
    [1 2 9];
    [23];
    [24];
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