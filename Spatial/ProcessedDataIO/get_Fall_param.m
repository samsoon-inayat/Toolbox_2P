function deconv = get_Fall_param(recording,plane)
filename = 'Fall.mat';
if isstruct(recording)
    filename = fullfile(recording.s2p_processed_data_folder,'suite2p',plane,filename);
else
    filename = fullfile(recording,filename);
end
mf = matfile(filename);
deconv = mf.ops;

