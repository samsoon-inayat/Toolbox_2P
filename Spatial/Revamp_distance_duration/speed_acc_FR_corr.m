function speed_acc_FR_corr

%% Load Data
tic
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei = evalin('base','ei');
udata = evalin('base','udata');
udata1 = evalin('base','udata1');
udataT = evalin('base','udataT');
udataD = evalin('base','udataD');

n  = 0;

%%
% Define number of animals
num_animals = length(udataT);

% Initialize storage for correlation results
correlation_results = cell(num_animals, 1); % Store per animal

% Loop through each animal
for a = 1:num_animals
    % Extract firing rate and speed
    firing_rate = udataT{a}.firing_rate; % [Neurons x Time Bins]
    speed = udataT{a}.speed; % Time-binned speed [1 x Time Bins]
    ts = udataT{a}.ts; % Time vector

    % Compute derivatives
    acceleration = diff(speed) ./ diff(ts); % Acceleration = d(Speed)/dt
    acceleration = [acceleration NaN]; % Append NaN to match time points

    delta_FR = diff(firing_rate, 1, 2) ./ repmat(diff(ts), size(firing_rate, 1), 1); % Delta FR = d(FR)/dt
    delta_FR = [delta_FR NaN(size(firing_rate, 1), 1)]; % Append NaN for alignment

    % Store derivatives in udataT
    udataT{a}.acceleration = acceleration;
    udataT{a}.delta_FR = delta_FR;

    % Initialize storage for correlations
    num_neurons = size(firing_rate, 1);
    corr_results = NaN(num_neurons, 4); % 4 correlation types

    % Compute correlations for each neuron
    for n = 1:num_neurons
        % Correlation 1: Speed vs. Firing Rate
        corr_results(n, 1) = corr(speed', firing_rate(n, :)', 'Type', 'Pearson','Rows', 'complete');

        % Correlation 2: Speed vs. Delta FR
        corr_results(n, 2) = corr(speed', delta_FR(n, :)', 'Type', 'Pearson','Rows', 'complete');

        % Correlation 3: Acceleration vs. Firing Rate
        corr_results(n, 3) = corr(acceleration', firing_rate(n, :)', 'Type', 'Pearson','Rows', 'complete');

        % Correlation 4: Acceleration vs. Delta FR
        corr_results(n, 4) = corr(acceleration', delta_FR(n, :)', 'Type', 'Pearson','Rows', 'complete');
    end

    % Store correlation results
    correlation_results{a} = corr_results;
end

% Display example results
disp('Correlation Results (Example from Animal 1, First 10 Neurons):');
disp(array2table(correlation_results{1}(1:10, :), ...
    'VariableNames', {'Speed-FR', 'Speed-DeltaFR', 'Accel-FR', 'Accel-DeltaFR'}));

%%
% Extract correlation values for all neurons (Example: First Animal)
corr_values = correlation_results{1}; % First animal's correlation data

% Compute correlation between Speed-FR and Acceleration-Delta FR
[r, p] = corr(corr_values(:, 1), corr_values(:, 4), 'Rows', 'complete');

% Display result
fprintf('Correlation between Speed-FR and Acceleration-Delta FR: r = %.3f, p = %.3f\n', r, p);

%%
figure(1000);clf;
subplot(2,3,1);scatter(corr_values(:,1),corr_values(:,2))
subplot(2,3,2);scatter(corr_values(:,1),corr_values(:,3))
subplot(2,3,3);scatter(corr_values(:,1),corr_values(:,4))
subplot(2,3,3);scatter(corr_values(:,2),corr_values(:,3))
subplot(2,3,4);scatter(corr_values(:,2),corr_values(:,4))
subplot(2,3,5);scatter(corr_values(:,3),corr_values(:,4))
%%
corr_res = cell(1,5);
for ii = 1:5
    tcorR = correlation_results{ii};
    corr_res{ii}.speed_FR = tcorR(:,1)';
    corr_res{ii}.speed_dFR = tcorR(:,2)';
    corr_res{ii}.accel_FR = tcorR(:,3)';
    corr_res{ii}.accel_dFR = tcorR(:,4)';
end
save('correlation_results.mat', 'corr_res', '-v7.3');

