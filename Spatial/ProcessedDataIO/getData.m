function tei = getData (f,T)

for ii = 1:size(T,1)
    f.recordingFolder = cell2mat(T{ii,6});
    disp(f.recordingFolder);
    db = get_db(f);
    tei{ii}.recordingFolder = f.recordingFolder;
    tei{ii}.thorExp = thorGetExperimentInfo(f.recordingFolder);
    tei{ii}.db = db;
    plane{1} = cell2mat(T{ii,8});
    plane{2} = cell2mat(T{ii,9});
    for pp = 1:length(plane)
        if isempty(plane{pp})
            continue;
        end
        sel_plane = plane{pp};
        fileName = makeName('behavior.mat',sel_plane);
        disp(sprintf('Loading 2P plane %d behavior',pp));
        b = load(fileName);
        tei{ii}.plane{pp}.b.frames_f = b.frames_f;
        if pp == 1
            if length(b.frames_f) < tei{ii}.thorExp.timePoints
                tei{ii}.thorExp.timePoints = length(b.frames_f)-1;
                tei{ii}.thorExp.totalFrames = length(b.frames_f)-1;
                tei{ii}.thorExp.streaming_frames = (length(b.frames_f)-1) * (tei{ii}.thorExp.zSteps + 1);
            end
        end
        tei{ii}.plane{pp}.folder = plane{pp};
        files = dir(sprintf('%s\\F_*_proc.mat',sel_plane));
        fileName = fullfile(files(1).folder,files(1).name);
        disp(sprintf('Loading 2P plane %d %s',pp,sel_plane));
        temp = load(fileName);
        tP = temp.dat;
        fileName = sprintf('%s\\Fcell_baseline.mat',sel_plane);
        if exist(fileName,'file')
            disp(sprintf('Loading 2P plane %d baselines',pp));
            temp = load(fileName);
            tP.Fcell_baseline{1} = temp.Fcell_baseline;
        else
            disp(sprintf('Finding 2P plane %d baselines',pp));
            Fcell_baseline = findBleachingTrend(double(tP.Fcell{1}'),4)';
            save(fileName,'Fcell_baseline');
            tP.Fcell_baseline{1} = Fcell_baseline;
        end
        tP.signals = 100*(double(tP.Fcell{1})-tP.Fcell_baseline{1})./tP.Fcell_baseline{1};
        areCells = zeros(length(tP.stat),1);
        for jj = 1:length(tP.stat)
            if tP.stat(jj).iscell && sum(isnan(tP.signals(jj,:)))==0
                areCells(jj,1) = 1;
            end
        end
        tP.areCells = areCells;
        disp(sprintf('Loading 2P plane %d spikes',pp));
        [tP.deconv.caSigAll,tP.deconv.spSigAll] = getSpikes(tei{ii}.plane{pp},tP);
        tei{ii}.plane{pp}.tP = tP;
        fileCheck = sprintf('%s\\reg*.mat',sel_plane);
        files = dir(fileCheck);
        fileName = fullfile(files(1).folder,files(1).name);
        temp1 = load(fileName);
        tei{ii}.plane{pp}.ops1 = temp1.ops1;
    end
    b = rmfield(b,'frames_f');
    b = rmfield(b,'frames_r');
    b = calcBehav(b);
    tei{ii}.b = b;
end
% return;
% ei.folders = f;
% db = get_db(f);
% ei.db = db(recording_number);
% ii = recording_number;
% posstr = strfind(f.pd_path,'Processed_Data\');
% subFolder = f.pd_path((posstr+15):end);
% ei.folders.rawDataFolder = f.recordingFolder;
% ei.folders.thisTifFolder = fullfile(ei.folders.tifFolder,subFolder);
% ei.folders.thispFolder = f.pd_path
% ei.folders.animalFolder = makeName(db(ii).mouse_name,ei.folders.pFolder);
% if ~exist(ei.folders.thispFolder,'dir')
%     mkdir(ei.folders.thispFolder);
% end
% if ~exist(ei.folders.thisTifFolder,'dir')
%     mkdir(ei.folders.thisTifFolder);
% end
%   
% ei = thorGetExperimentInfo(ei);
% ei.animal_id = animal_id;
% ei.exp_date = exp_date;
% ei.recording_number = recording_number;
% 
% 
% if p.Results.twoP == 1
%     display('Loading 2P data');
%     fileCheck = sprintf('%s\\F_*_proc.mat',ei.folders.thispFolder);
%     files = dir(fileCheck);
%     fileName = makeName(files(1).name,ei.folders.thispFolder);
%     temp = load(fileName);
%     tP = temp.dat;
%     fileName = sprintf('%s\\Fcell_baseline.mat',ei.folders.thispFolder);
%     if exist(fileName,'file')
%         temp = load(fileName);
%         tP.Fcell_baseline{1} = temp.Fcell_baseline;
%     else
%         Fcell_baseline = findBleachingTrend(double(tP.Fcell{1}'),4)';
%         save(fileName,'Fcell_baseline');
%         tP.Fcell_baseline{1} = Fcell_baseline;
%     end
%     tP.signals = 100*(double(tP.Fcell{1})-tP.Fcell_baseline{1})./tP.Fcell_baseline{1};
%     ei.tP = tP;
% 
%     fileCheck = sprintf('%s\\reg*.mat',ei.folders.thispFolder);
%     files = dir(fileCheck);
%     fileName = makeName(files(1).name,ei.folders.thispFolder);
%     temp1 = load(fileName);
%     ei.ops1 = temp1.ops1;
% 
%     areCells = [];
%     for ii = 1:length(tP.stat)
%         if tP.stat(ii).iscell && sum(isnan(tP.signals(ii,:)))==0
%             areCells = [areCells ii];
%         end
%     end
%     ei.areCells = areCells';
%     ei.signals = tP.signals(areCells,:);
% end
% 
% if p.Results.twoP == 2
%     display('Loading 2P data');
%     fileCheck = sprintf('%s\\F_*_proc.mat',ei.folders.thispFolder);
%     files = dir(fileCheck);
%     fileName = makeName(files(1).name,ei.folders.thispFolder);
%     temp = load(fileName);
%     tP = temp.dat;
%     fileName = sprintf('%s\\Fcell_baseline.mat',ei.folders.thispFolder);
%     if exist(fileName,'file')
%         temp = load(fileName);
%         tP.Fcell_baseline{1} = temp.Fcell_baseline;
%     else
%         Fcell_baseline = findBleachingTrend(double(tP.Fcell{1}'),4)';
%         save(fileName,'Fcell_baseline');
%         tP.Fcell_baseline{1} = Fcell_baseline;
%     end
%     tP.signals = 100*(double(tP.Fcell{1})-tP.Fcell_baseline{1})./tP.Fcell_baseline{1};
%     ei.tP = tP;
% 
%     fileCheck = sprintf('%s\\reg*.mat',ei.folders.thispFolder);
%     files = dir(fileCheck);
%     fileName = makeName(files(1).name,ei.folders.thispFolder);
%     temp1 = load(fileName);
%     ei.ops1 = temp1.ops1;
% 
%     areCells = [];
%     for ii = 1:length(tP.stat)
%         if tP.stat(ii).iscell && sum(isnan(tP.signals(ii,:)))==0
%             areCells = [areCells ii];
%         end
%     end
%     ei.areCells = areCells';
%     ei.signals = tP.signals(areCells,:);
% end
% 
% 
% function trials = selectTrials(b)
% thresh = 0.07; % in in micro secs
% for ii = 1:length(b.trials)
%     trialStart = b.air_puff_r(ii);
%     trialEnd = b.air_puff_f(ii);
%     speeds = b.speed(trialStart:trialEnd);
%     aSpeeds(ii) = mean(speeds);
% end
% trials = find(aSpeeds>thresh);
% n=0;