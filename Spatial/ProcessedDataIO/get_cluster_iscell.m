function deconv = get_cluster_iscell(recording,plane)
filename = 'Fall.mat';
filename = fullfile(recording.s2p_processed_data_folder,'suite2p',plane,filename);
mf = matfile(filename);
deconv = mf.iscell;

