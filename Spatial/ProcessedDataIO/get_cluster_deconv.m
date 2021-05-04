function deconv = get_cluster_deconv(recording)
filename = 'deconv.mat';
filename = fullfile(recording.processed_data_folder,filename);
temp = load(filename);
deconv = temp.deconv;

