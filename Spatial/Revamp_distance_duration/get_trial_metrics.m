function out = get_trial_metrics(data_an,type,bw,air_phase,conf)

% from the data_an variable, get all the fields as independent variables
field_names = fieldnames(data_an);
for ii = 1:length(field_names)
    varname = field_names{ii};
    cmdTxt = sprintf('%s = data_an.%s;',varname,varname);
    eval(cmdTxt);
    cmdTxt = sprintf('a_%s = [];',varname);    eval(cmdTxt);
end

% select the appropriate air phase and configuration number
% cmdTxt = sprintf('csel = air & %s;',conf); eval(cmdTxt);
frf_n = double(frf);
num_frames = sum(frf==1);
frf_n(frf==1) = 1:num_frames;
% select the air_trials
if strcmp(air_phase,'ON')
    cmdTxt = sprintf('air_trials = air_trials_on .* %s;',conf); eval(cmdTxt);
end
if strcmp(air_phase,'OFF')
    cmdTxt = sprintf('air_trials = air_trials_off .* %s;',conf); eval(cmdTxt);
end


ds(ds<0) = NaN; var = fillmissing(ds, 'spline');
dds = diff(ds); spike_idx = find(abs(dds) > 50) + 1;
ds(spike_idx) = NaN; ds(spike_idx-1) = NaN; ds(spike_idx+1) = NaN; ds = fillmissing(ds, 'spline');

bins = 0:bw:1000; % set a large number of bins
a_tts = []; a_tds = []; % the first one will accumulate the trial times where the first time from the start of trial will be 0, same thing for the second which will collect trial distance
a_bins = []; a_trials = []; a_bindur = []; a_bindist = [];
firing_rate = data_an.firing_rate;
for tn = 1:10
    trials = (air_trials == tn); % get the trials as a big vector containing 1s corresponding to the trial in question defined by ii (iterated)
    trialH = find(trials); % get the indices of the trials
    trial = find(frf & trials); % get the indices of the trials
    tts = ts(trial)-ts(trialH(1)); tds = ds(trial)-ds(trialH(1)); % for the tds and tts variables we want their first value to be zero i.e., start of the trial
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
            frame_nums = frf_n(frf & trials);
            FR = firing_rate(:,frame_nums);
            % Parallelized loop to calculate the binned firing rate for each neuron
            firing_rate_binned = [];
            parfor neuron_idx = 1:size(firing_rate,1)
                % Calculate the binned firing rate for each neuron using the frame-based bin indices
                firing_rate_binned(neuron_idx, :) = accumarray(bin_indices', FR(neuron_idx, :),[], @mean, NaN)';
            end
            a_firing_rate = [a_firing_rate firing_rate_binned];
            continue;
        else
            cmdTxt = sprintf('res = islogical(%s);',varname);eval(cmdTxt)
            if res
                cmdTxt = sprintf('%s_binned = accumarray(bin_indices'', %s(trial), [], @(x) any(x), false)'';',varname,varname); eval(cmdTxt);
            else
                cmdTxt = sprintf('%s_binned = accumarray(bin_indices'', %s(trial), [], @mean, NaN)'';',varname,varname); eval(cmdTxt);
            end
        end
        cmdTxt = sprintf('a_%s = [a_%s %s_binned];',varname,varname,varname); eval(cmdTxt);
    end
    total_time(tn,1) = tts_binned(end);
    total_distance(tn,1) = tds_binned(end);
    movement_latency(tn,1) = tts_binned(find(speed_binned > 0, 1, 'first'));
    rest_latency(tn,1) = tts_binned(find(speed_binned > 0, 1, 'last'));
    mean_speed(tn,1) = mean(speed_binned); 
    std_speed(tn,1) = std(speed_binned); 
    ske_speed(tn,1) = skewness(speed_binned); 
    kurt_speed(tn,1) = kurtosis(speed_binned); 
    cov_speed(tn,1) = mean(speed_binned)/std(speed_binned);
    n = 0;
end

for ii = 1:length(field_names)
    varname = field_names{ii};
    cmdTxt = sprintf('out.%s = a_%s;',varname,varname); eval(cmdTxt);
end

out.dist = a_tds;
out.time = a_tts;
out.binnum = a_bins;
out.trialnum = a_trials;
out.bindur = a_bindur;
out.bindist = a_bindist;
out.trial_metrics_code = {'time','distance','movement_latency','rest_latency','mean_speed','std_speed','cov_speed','ske_speed','kurt_speed'};
out.trial_metrics = [total_time,total_distance,rest_latency,movement_latency,mean_speed,std_speed,cov_speed,ske_speed,kurt_speed];
