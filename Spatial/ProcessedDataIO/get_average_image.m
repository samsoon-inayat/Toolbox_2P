function deconv = get_cluster_deconv(recording,plane)
if isstruct(recording)
    try
        filename = 'avgstack.mat';
        filename = fullfile(recording.s2p_processed_data_folder,filename);
        temp = load(filename);
        deconv = temp.deconv;
        return;
    catch
        filename = 'avgstack.mat';
        filename = fullfile(recording.s2p_processed_data_folder,'suite2p',plane,filename);
        temp = load(filename);
        deconv = temp.deconv;
    end
end

if ischar(recording)
    try
        filename = 'avgstack.mat';
        filename = fullfile(recording,filename);
        temp = load(filename);
        deconv = temp.avg_stack;
        return;
    catch
        filename = 'avgstack.mat';
        filename = fullfile(recording,'suite2p',plane,filename);
        temp = load(filename);
        deconv = temp.avg_stack;
    end
end

