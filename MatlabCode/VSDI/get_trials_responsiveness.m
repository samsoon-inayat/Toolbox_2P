function presp = get_trials_responsiveness (ags,pags,stim_folders_info,indices,sel_roi)
presp = [];
for ii = 1:length(indices)
    indx = indices(ii)
    presp = [presp processFolder(ags,pags,stim_folders_info,indx,sel_roi)];
end


function presp = processFolder(ags,pags,stim_folders_info,indx,sel_roi)
gg = stim_folders_info.list(indx,2);
aa = stim_folders_info.list(indx,3);
rr = stim_folders_info.list(indx,4);
thisRecording = ags(gg).animals{aa}.eRecordings{rr};
pThisRecording = pags(gg).animals{aa}.eRecordings{rr};
fileName = makeName('trials_responsiveness.mat',pThisRecording.root_folder);
resp = load(fileName);
tresp = [resp.FRVhijb(sel_roi,:) resp.FRVlojb(sel_roi,:)];
presp = 100*sum(tresp == 1)/length(tresp);





