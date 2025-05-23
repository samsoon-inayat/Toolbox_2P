function out = get_trial_metrics(data_an,type,bw,air_phase,conf)

% from the data_an variable, get all the fields as independent variables
field_names = fieldnames(data_an);
for ii = 1:length(field_names)
    varname = field_names{ii};
    cmdTxt = sprintf('%s = data_an.%s;',varname,varname);
    eval(cmdTxt);
    if strcmp(varname,'firing_rate')
        continue;
    end
    cmdTxt = sprintf('a_%s = [];',varname);    eval(cmdTxt);
end

% select the appropriate air phase and configuration number
% cmdTxt = sprintf('csel = air & %s;',conf); eval(cmdTxt);
num_frames = sum(frf==1);
frf_n(frf==1) = 1:num_frames;
% select the air_trials
if strcmp(air_phase,'ON')
    cmdTxt = sprintf('air_trials = air_trials_on .* %s;',conf); eval(cmdTxt);
end
if strcmp(air_phase,'OFF')
    cmdTxt = sprintf('air_trials = air_trials_off .* %s;',conf); eval(cmdTxt);
end

bins = 0:bw:1000; % set a large number of bins
a_tts = []; a_tds = []; % the first one will accumulate the trial times where the first time from the start of trial will be 0, same thing for the second which will collect trial distance
a_bins = []; a_trials = []; a_bindur = []; a_bindist = [];
for ii = 1:10
    trials = (air_trials == ii); % get the trials as a big vector containing 1s corresponding to the trial in question defined by ii (iterated)
    trial = find(frf & trials); % get the indices of the trials
    tts = ts(trial)-ts(trial(1)); tds = ds(trial)-ds(trial(1)); % for the tds and tts variables we want their first value to be zero i.e., start of the trial
    %     frf_ti = frf_n(frf & trials);
    %     FR = firing_rate(:,frf_ti);
    %     frf_t = ts(frf & trials)-ts(trial(1));
    % a_tts = [a_tts tts]; a_tds = [a_tds tds];
    if strcmp(type,'time')
        bin_indices = discretize(tts, bins);
    end
    if strcmp(type,'distance')
        tds(tds<0)=0;
        bin_indices = discretize(tds, bins);
    end
    tts_binned = accumarray(bin_indices',tts',[],@mean);    tts_binned = tts_binned - tts_binned(1);    a_tts = [a_tts tts_binned'];
    tds_binned = accumarray(bin_indices',tds',[],@mean);    tds_binned = tds_binned - tds_binned(1);    a_tds = [a_tds tds_binned'];
    a_bins = [a_bins 1:length(tts_binned)];    
    a_trials = [a_trials (ones(size(tts_binned))*ii)'];
    % a_bindur = [a_bindur accumarray(bin_indices',ones(size(bin_indices))*bw,[],@sum)'];
    a_bindur = [a_bindur accumarray(bin_indices',tts,[],@(x) max(x)-min(x))'];
    a_bindist = [a_bindist accumarray(bin_indices',tds',[],@(x) max(x)-min(x))'];
    % in the following loop we want to bin the values based on the
    % bin_indices. For the continuous variables the binned value will be
    % mean over bins
    for ii = 1:length(field_names)
        varname = field_names{ii};
        if strcmp(varname,'firing_rate')
            % if strcmp(type,'time')
            %     frf_ti = frf_n(frf & trials);
            %     FR = firing_rate(:,frf_ti);
            %     frf_t = ts(frf & trials)-ts(trial(1));
            %     frame_bin_indices = discretize(frf_t, bins);  % Bin the frame times
            % end
            % if strcmp(type,'distance')
            %     frf_t = ds(frf & trials)-ds(trial(1));
            %     frame_bin_indices = discretize(frf_t, bins);  % Bin the frame times
            % end
            % % Parallelized loop to calculate the binned firing rate for each neuron
            % parfor neuron_idx = 1:size(firing_rate,1)
            %     % Calculate the binned firing rate for each neuron using the frame-based bin indices
            %     firing_rate_binned(neuron_idx, :) = accumarray(frame_bin_indices', FR(neuron_idx, :),[], @mean, NaN)';
            % end
            continue;
        end
        cmdTxt = sprintf('res = islogical(%s);',varname);eval(cmdTxt)
        if res
            cmdTxt = sprintf('%s_binned = accumarray(bin_indices'', %s(trial), [], @(x) any(x), false)'';',varname,varname); eval(cmdTxt);
        else
            cmdTxt = sprintf('%s_binned = accumarray(bin_indices'', %s(trial), [], @mean, NaN)'';',varname,varname); eval(cmdTxt);
        end
        cmdTxt = sprintf('a_%s = [a_%s %s_binned];',varname,varname,varname); eval(cmdTxt);
    end
end

for ii = 1:length(field_names)
    varname = field_names{ii};
    if strcmp(varname,'firing_rate')
        continue;
    end
    cmdTxt = sprintf('out.%s = a_%s;',varname,varname); eval(cmdTxt);
end

out.dist = a_tds;
out.time = a_tts;
out.binnum = a_bins;
out.trialnum = a_trials;
out.bindur = a_bindur;
out.bindist = a_bindist;


% Creating the model for firing rate as a function of speed, distance, time, and motion
% X = [atts', atsp', atds'];%, atac', attm'];    % Predictor matrix (continuous variables)
% out.time = atts';
% out.dist = atds';
% out.speed = atsp';
% out.FR = aFR';
% out.trialnum = attn';
% out.binnum = atbn';
% out.bindur = atbdur;
% out.bindist = atbdist;
[out.trial_metrics_code,out.trial_metrics] = get_metrics(out);

% 
function [trial_metrics_code,trial_metrics] = get_metrics(out)
% Define functions for different variance metrics
mean_fun = @(x) mean(x);
std_fun = @(x) std(x, 0);  % Standard deviation (unbiased)
cv_fun = @(x) std(x, 0) ./ mean(x);  % Coefficient of Variation (CV)
skew_fun = @(x) skewness(x);  % Skewness
kurt_fun = @(x) kurtosis(x);  % Kurtosis
max_fun = @(x) max(x); 
latency_fun = @(x, t) t(find(x > 0, 1, 'first'));  
rlatency_fun = @(x, t) t(find(x > 0, 1, 'last'));  

trial_nums = out.trialnum; speed = out.speed; time = out.time; distance = out.dist; bindur = out.bindur; bindist = out.bindist;
% Compute standard deviation, CV, skewness, and kurtosis for each trial
trial_mean = accumarray(trial_nums', speed, [], mean_fun); 
trial_std = accumarray(trial_nums', speed, [], std_fun); 
trial_cv = accumarray(trial_nums', speed, [], cv_fun);
trial_skew = accumarray(trial_nums', speed, [], skew_fun);
trial_kurt = accumarray(trial_nums', speed, [], kurt_fun);
% Compute total time and total distance for each trial
last_bin_time = accumarray(trial_nums', time, [], max_fun);
total_time = accumarray(trial_nums', bindur, [], @sum);
total_distance = accumarray(trial_nums', bindist, [], @sum);
last_bin_distance = accumarray(trial_nums', distance, [], max_fun);
% Compute movement latency for each trial
movement_latency = accumarray(trial_nums', speed, [], @(x) latency_fun(x, time));
rest_latency = accumarray(trial_nums', speed, [], @(x) latency_fun(x, time));

trial_metrics_code = {'time','distance','movement_latency','rest_latency','mean_speed','std_speed','cov_speed','ske_speed','kurt_speed'};
trial_metrics = [total_time,total_distance,rest_latency,movement_latency,trial_mean,trial_std,trial_cv,trial_skew,trial_kurt];