function data = get_unlv_analysis_data_bin (udata,binwidth,type)
%%
for ani = 1:5
    disp(ani)
   data{ani} = return_data_animal(udata{ani},binwidth,type);
end




function o = return_data_animal(tudata,bw,type)

field_names = fieldnames(tudata);
for ii = 1:length(field_names)
    varname = field_names{ii};
    cmdTxt = sprintf('%s = tudata.%s;',varname,varname);
    eval(cmdTxt);
end
% firing_rate = tudata.firing_rate;
% Define time bins (every 0.3 seconds)
if strcmp(type,'time')
    var = ts;
    disp(sprintf('Time - min %.3f, max %.3f',min(var),max(var)));
    bins = 0:bw:(1.5*max(ts));  % Time bins from 0 to max time in ts variable
end
if strcmp(type,'distance')
    var = ds; var(var<0) = NaN; var = fillmissing(var, 'spline');
    dvar = diff(var); spike_idx = find(abs(dvar) > 50) + 1;
    var(spike_idx) = NaN; var(spike_idx-1) = NaN; var(spike_idx+1) = NaN; var = fillmissing(var, 'spline');
    % figure(100);clf;plot(dvar);hold on;plot(var);
    disp(sprintf('Distance - min %.3f, max %.3f',min(var),max(var)));
    bins = 0:bw:(1.5*max(ds)); % Time bins from 0 to max time in ts variable
end

frf = tudata.frf;
firing_rate = tudata.firing_rate;
% initialize variables in memory
for ii = 1:length(field_names)
    varname = field_names{ii};
    cmdTxt = sprintf('%s_binned = NaN(1, length(bins)-1);',varname);
    eval(cmdTxt);
end


res = false;

% Discretize the time points into the defined bins
[bin_indices,edges] = discretize(var, bins);  % bin_indices maps each var value to a time bin
% Ensure no NaN or 0 indices in bin_indices
% bin_indices(isnan(bin_indices)) = length(bins)-1;  % Assign NaNs to the last valid bin
% % % % % % % % bin_indices(bin_indices == 0) = 1;    % Replace zero indices with first bin
% bin_indices = fillmissing(bin_indices,"next");


% Initialize matrices to store binned data
num_neurons = size(firing_rate, 1);  % Number of neurons
firing_rate_binned = NaN(num_neurons, length(bins)-1);  % Preallocate binned firing rate

% Use frf to restrict firing rate binning to frames (assuming frf is a logical index)
frame_times = var(frf);  % Get the actual frame times corresponding to 'frf'
frame_bin_indices = discretize(frame_times, bins);  % Bin the frame times

% % % % % % % Ensure frame_bin_indices don't exceed valid bounds
% % % % % % frame_bin_indices(frame_bin_indices > length(bins)-1) = length(bins)-1;
% % % % % % frame_bin_indices(frame_bin_indices == 0) = 1;  % Ensure no 0 index
% % % % % % frame_bin_indices(isnan(frame_bin_indices)) = length(bins)-1;  % Assign NaNs to the last valid bin
% frame_bin_indices(isnan(frame_bin_indices)) = 0;  % Assign NaNs to the last valid bin
% frame_bin_indices = fillmissing(frame_bin_indices,"next");


% Parallelized loop to calculate the binned firing rate for each neuron
parfor neuron_idx = 1:num_neurons
    % Calculate the binned firing rate for each neuron using the frame-based bin indices
    firing_rate_binned(neuron_idx, :) = accumarray(frame_bin_indices', firing_rate(neuron_idx, :),[length(bins)-1, 1], @mean, NaN)';
end

% Use accumarray to bin the speed
% speed_binned = accumarray(bin_indices', speed, [length(bins)-1, 1], @mean, NaN)';

% Use accumarray to bin the bnb (brake status)
% bnb_binned = accumarray(bin_indices', bnb, [length(bins)-1, 1], @(x) any(x), NaN)';


for ii = 1:length(field_names)
    varname = field_names{ii};
    if strcmp(varname,'firing_rate')
        continue;
    end
    cmdTxt = sprintf('res = islogical(%s);',varname);eval(cmdTxt)
    if res
        cmdTxt = sprintf('%s_binned = accumarray(bin_indices'', %s, [length(bins)-1, 1], @(x) any(x), false)'';',varname,varname); eval(cmdTxt);
    else
        cmdTxt = sprintf('%s_binned = accumarray(bin_indices'', %s, [length(bins)-1, 1], @mean, NaN)'';',varname,varname); eval(cmdTxt);
        % speed_binned = accumarray(bin_indices', speed, [length(bins)-1, 1], @mean, NaN)';
    end
end

vars_to_consider = {'bnb','air','light','C1','C2','C3','C4','C5','C6','C7','air_trials_on','air_trials_off'};
if strcmp(type,'distance')
    for ii = 1:length(vars_to_consider)
        if ii == 11
            n = 0;
        end
        cmdTxt = sprintf('temp = %s_binned;',vars_to_consider{ii}); eval(cmdTxt);
        if temp(1) > 0
            temp(1) = 0;
            cmdTxt = sprintf('%s_binned = temp;',vars_to_consider{ii}); eval(cmdTxt);
        end
    end
end


for ii = 1:length(field_names)
    varname = field_names{ii};
    cmdTxt = sprintf('o.%s = %s_binned;',varname,varname);
    eval(cmdTxt);
end

% Compute bin duration (difference of ts within each bin)
o.bindur = (accumarray(bin_indices', ts', [], @(x) max(x) - min(x)))';
 % Compute bin distance (difference of ds within each bin)
o.bindist = (accumarray(bin_indices', ds', [], @(x) max(x) - min(x)))';


n = 0;
% % Display the binned results for verification
% disp('Binned firing rates, speed, and logical variables (bnb, C1, C2):');
% disp(firing_rate_binned);
% disp(speed_binned);
% disp(bnb_binned);
% disp(C1_binned);
% disp(C2_binned);
