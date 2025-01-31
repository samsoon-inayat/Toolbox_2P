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
    acceleration = [acceleration; NaN]; % Append NaN to match time points

    delta_FR = diff(firing_rate, 1, 2) ./ diff(ts)'; % Delta FR = d(FR)/dt
    delta_FR = [NaN(size(firing_rate, 1), 1), delta_FR]; % Append NaN for alignment

    % Store derivatives in udataT
    udataT{a}.acceleration = acceleration;
    udataT{a}.delta_FR = delta_FR;

    % Initialize storage for correlations
    num_neurons = size(firing_rate, 1);
    corr_results = NaN(num_neurons, 4); % 4 correlation types

    % Compute correlations for each neuron
    for n = 1:num_neurons
        % Correlation 1: Speed vs. Firing Rate
        corr_results(n, 1) = corr(speed', firing_rate(n, :)', 'Type', 'Pearson');

        % Correlation 2: Speed vs. Delta FR
        corr_results(n, 2) = corr(speed', delta_FR(n, :)', 'Type', 'Pearson');

        % Correlation 3: Acceleration vs. Firing Rate
        corr_results(n, 3) = corr(acceleration, firing_rate(n, :)', 'Type', 'Pearson');

        % Correlation 4: Acceleration vs. Delta FR
        corr_results(n, 4) = corr(acceleration, delta_FR(n, :)', 'Type', 'Pearson');
    end

    % Store correlation results
    correlation_results{a} = corr_results;
end

% Display example results
disp('Correlation Results (Example from Animal 1, First 10 Neurons):');
disp(array2table(correlation_results{1}(1:10, :), ...
    'VariableNames', {'Speed-FR', 'Speed-DeltaFR', 'Accel-FR', 'Accel-DeltaFR'}));
