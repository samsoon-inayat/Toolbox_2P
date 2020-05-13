function stim_info = get_Files_Info_SP (ags,pags)
stim_info.listLabels = {'stimulus (ss)','Group (gg)','AnimalNumber (aa)','RecordingNumber (rr)'};
stim_info.list = [];
N_groups = length(ags);
for gg = 1:N_groups
    thisGroup = ags(gg);
    N_animals = length(thisGroup.animals);
    for aa = 1:N_animals
        thisAnimal = thisGroup.animals{aa};
        N_recordings = length(thisAnimal.sRecordings);
        for rr = 1:N_recordings
            thisRecording = thisAnimal.sRecordings{rr};
            stim_info.list = [stim_info.list;[0 gg aa rr]];
            stim_info.root_folders{size(stim_info.list,1),1} = thisRecording.root_folder;
            stim_info.p_root_folders{size(stim_info.list,1),1} = pags(gg).animals{aa}.sRecordings{rr}.root_folder;
        end
    end
end


