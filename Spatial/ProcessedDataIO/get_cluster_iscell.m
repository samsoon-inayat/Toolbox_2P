function deconv = get_cluster_iscell(recording)
filename = 'Fall.mat';
filename = fullfile(recording.processed_data_folder,'suite2p','plane0',filename);
mf = matfile(filename);
deconv = mf.iscell;

