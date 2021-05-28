function deconv = get_cluster_deconv(recording,plane)
try
    filename = 'deconv.mat';
    filename = fullfile(recording.s2p_processed_data_folder,filename);
    temp = load(filename);
    deconv = temp.deconv;
    return;
catch
    filename = 'deconv.mat';
    filename = fullfile(recording.s2p_processed_data_folder,'suite2p',plane,filename);
    temp = load(filename);
    deconv = temp.deconv;
end


