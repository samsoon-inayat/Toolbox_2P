function deconv = get_cluster_timecourses(recording)
if isstruct(recording)
filename = 'timecourses.mat';
filename = fullfile(recording.processed_data_folder,filename);
temp = load(filename);
deconv = temp.tcs;
end

if ischar(recording)
    filename = 'timecourses.mat';
    filename = fullfile(recording,filename);
    temp = load(filename);
    deconv = temp.tcs;
end

