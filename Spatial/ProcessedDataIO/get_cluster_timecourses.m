function deconv = get_cluster_timecourses(recording)
filename = 'timecourses.mat';
filename = fullfile(recording.processed_data_folder,filename);
temp = load(filename);
deconv = temp.tcs;

