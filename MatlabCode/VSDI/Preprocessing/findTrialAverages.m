function findTrialAverages (animal)
for an = 1:length(animal.data_folders)
    peDataFolder = animal.processed_data_folders{an}
    fileName = makeName('folders_data.mat',peDataFolder);
    temp = load(fileName);
    folders = temp.folders_data; clear temp;
    for ii = 1:length(folders)
        trials = folders(ii).stimulus.actualTrials;
        evFolder = makeName(folders(ii).name,peDataFolder);
        dfbyfo = [];
        for jj = 1:trials
            trialName = sprintf('trial_%.2d',jj);
            fileName = makeName(trialName,evFolder);
            fileName = makeName('dfbyfo.mat',fileName);
            temp = load(fileName);
            if jj == 1
                dfbyfo = temp.dfbyfo;
            else
                dfbyfo = dfbyfo + temp.dfbyfo;
            end
        end
        dfbyfo = dfbyfo/trials;
        fileName = makeName('dfbyfo.mat',evFolder);
        save(fileName,'dfbyfo');
    end
end
