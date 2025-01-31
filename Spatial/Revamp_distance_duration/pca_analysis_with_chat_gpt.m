
tR = Rs_C{1,3};

firing_rate = tR.sp_rasters1;
firing_rate = fillmissing(firing_rate, 'linear', 1);  % Replace NaNs along the 1st dimension (trials)
%%

% Step 1: Center the Data by Subtracting the Mean
% firing_rate is an MxNxP matrix (M=trials, N=bins, P=cells)
mean_firing_rate = mean(firing_rate, [1, 2], 'omitnan');  % Mean across trials and bins for each cell

% Center the data by subtracting the mean
centered_firing_rate = firing_rate - repmat(mean_firing_rate, size(firing_rate, 1), size(firing_rate, 2), 1);

% Step 2: Reshape the Data for PCA
% Reshape the centered data to a 2D matrix for PCA (collapse trials and bins into one dimension)
reshaped_data = reshape(centered_firing_rate, [], size(firing_rate, 3));  % Collapse trials and bins (MxN)

% Step 3: Perform PCA
[coeff, score, latent] = pca(reshaped_data');  % Perform PCA, score gives the projection

% Step 4: Reduce the Dimensionality
num_components = 50;  % Choose the number of components to retain
reduced_data = score(:, 1:num_components);  % Project onto the first num_components PCs

% Step 5: Reconstruct the Data
% Reconstruct the data with the reduced components
reconstructed_data = reduced_data * coeff(:, 1:num_components)';  % Re-project into original space

% Reshape back to the original dimensions (MxNxP)
reconstructed_data = reshape(reconstructed_data', size(firing_rate));  % Make sure to transpose to match dimensions

% Step 6: Re-add the Mean
reconstructed_data = reconstructed_data + repmat(mean_firing_rate, size(firing_rate, 1), size(firing_rate, 2), 1);

% Optional: Visualize the Reconstructed Data
figure(1000);clf
imagesc(squeeze(mean(reconstructed_data, 1)));  % Example visualization of the mean over trials
colorbar;
title('Reconstructed Data with Reduced Dimensions');

%%
% Compute the average firing rate across trials for each time bin (after noise reduction)
avg_firing_rate = mean(reconstructed_data, 1);

% Plot tuning curve for a specific neuron
neuron_idx = 8;
figure(1000);clf
plot(squeeze(avg_firing_rate(:, :, neuron_idx)));
title(['Tuning Curve for Neuron ', num2str(neuron_idx)]);
xlabel('Time Bins');
ylabel('Firing Rate');

%%
% Example: Compute coherence between two neurons
neuron1 = squeeze(firing_rate(:, :, 2));  % Firing rates for neuron 1
neuron2 = squeeze(firing_rate(:, :, 3));  % Firing rates for neuron 2

% Compute coherence
[Cxy, f] = mscohere(neuron1(:), neuron2(:));
figure(1000);clf
plot(f, Cxy);
xlabel('Frequency (Hz)');
ylabel('Coherence');
title('Coherence between Neuron 1 and Neuron 2');