function data = get_unlv_analysis_data_time_bin (udata,binwidth)
%%
for ani = 1:5
    disp(ani)
   data{ani} = return_data_animal(udata{ani},binwidth);
end




function o = return_data_animal(tudata,bw)

field_names = fieldnames(tudata);
for ii = 1:length(field_names)
    varname = field_names{ii};
    cmdTxt = sprintf('%s = tudata.%s;',varname,varname);
    eval(cmdTxt);
end
% firing_rate = tudata.firing_rate;
% Define time bins (every 0.3 seconds)
time_bin_size = bw;  % 0.3 seconds
time_bins = 0:time_bin_size:max(ds);  % Time bins from 0 to max time in ds variable
frf = tudata.frf;
firing_rate = tudata.firing_rate;
% initialize variables in memory
for ii = 1:length(field_names)
    varname = field_names{ii};
    cmdTxt = sprintf('%s_binned = NaN(1, length(time_bins)-1);',varname);
    eval(cmdTxt);
end


res = false;

% Discretize the time points into the defined bins
bin_indices = discretize(ds, time_bins);  % bin_indices maps each ds value to a time bin
% Ensure no NaN or 0 indices in bin_indices
bin_indices(isnan(bin_indices)) = length(time_bins)-1;  % Assign NaNs to the last valid bin
bin_indices(bin_indices == 0) = 1;    % Replace zero indices with first bin

% Initialize matrices to store binned data
num_neurons = size(firing_rate, 1);  % Number of neurons
firing_rate_binned = NaN(num_neurons, length(time_bins)-1);  % Preallocate binned firing rate

% Use frf to restrict firing rate binning to frames (assuming frf is a logical index)
frame_times = ds(frf);  % Get the actual frame times corresponding to 'frf'
frame_bin_indices = discretize(frame_times, time_bins);  % Bin the frame times

% Ensure frame_bin_indices don't exceed valid bounds
frame_bin_indices(frame_bin_indices > length(time_bins)-1) = length(time_bins)-1;
frame_bin_indices(frame_bin_indices == 0) = 1;  % Ensure no 0 index
frame_bin_indices(isnan(frame_bin_indices)) = length(time_bins)-1;  % Assign NaNs to the last valid bin


% Parallelized loop to calculate the binned firing rate for each neuron
parfor neuron_idx = 1:num_neurons
    % Calculate the binned firing rate for each neuron using the frame-based bin indices
    firing_rate_binned(neuron_idx, :) = accumarray(frame_bin_indices', firing_rate(neuron_idx, :), ...
        [length(time_bins)-1, 1], @mean, NaN)';
end

% Use accumarray to bin the speed
% speed_binned = accumarray(bin_indices', speed, [length(time_bins)-1, 1], @mean, NaN)';

% Use accumarray to bin the bnb (brake status)
% bnb_binned = accumarray(bin_indices', bnb, [length(time_bins)-1, 1], @(x) any(x), NaN)';


for ii = 1:length(field_names)
    varname = field_names{ii};
    if strcmp(varname,'firing_rate')
        continue;
    end
    cmdTxt = sprintf('res = islogical(%s);',varname);eval(cmdTxt)
    if res
        cmdTxt = sprintf('%s_binned = accumarray(bin_indices'', %s, [length(time_bins)-1, 1], @(x) any(x), false)'';',varname,varname); eval(cmdTxt);
    else
        cmdTxt = sprintf('%s_binned = accumarray(bin_indices'', %s, [length(time_bins)-1, 1], @mean, NaN)'';',varname,varname); eval(cmdTxt);
        % speed_binned = accumarray(bin_indices', speed, [length(time_bins)-1, 1], @mean, NaN)';
    end
end

for ii = 1:length(field_names)
    varname = field_names{ii};
    cmdTxt = sprintf('o.%s = %s_binned;',varname,varname);
    eval(cmdTxt);
end

% % Display the binned results for verification
% disp('Binned firing rates, speed, and logical variables (bnb, C1, C2):');
% disp(firing_rate_binned);
% disp(speed_binned);
% disp(bnb_binned);
% disp(C1_binned);
% disp(C2_binned);
