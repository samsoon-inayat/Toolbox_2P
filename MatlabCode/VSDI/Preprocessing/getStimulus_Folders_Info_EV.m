function stim_info = getStimulus_Folders_Info_EV (ags,stimulus_types,pags)
stim_info.listLabels = {'stimulus (ss)','Group (gg)','AnimalNumber (aa)','RecordingNumber (rr)'};
stim_info.list = [];
for ss = 1:length(stimulus_types)
    N_groups = length(ags);
    for gg = 1:N_groups
        thisGroup = ags(gg);
        N_animals = length(thisGroup.animals);
        for aa = 1:N_animals
            thisAnimal = thisGroup.animals{aa};
            N_recordings = length(thisAnimal.eRecordings);
            for rr = 1:N_recordings
                thisRecording = thisAnimal.eRecordings{rr};
                if strcmp(thisRecording.stimulus.stimulus_type,stimulus_types{ss})
                    stim_info.list = [stim_info.list;[ss gg aa rr]];
                    stim_info.root_folders{size(stim_info.list,1),1} = thisRecording.root_folder;
                    stim_info.p_root_folders{size(stim_info.list,1),1} = pags(gg).animals{aa}.eRecordings{rr}.root_folder;
                end
            end
        end
    end
end




