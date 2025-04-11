function time_trial_distance_speed

udata = evalin('base','udata1');
ei = evalin('base','ei');
configurations = {'C3','C4','C5'};
air_phases = {'ON','OFF'};

time_bins = 0:0.15:1000; % set a large number of bins
atimecc = {}; adistcc = {}; aspeedcc = {}; aFRcc = {}; atrialcc = {};  o_atimecc = {}; o_adistcc = {}; o_aspeedcc = {}; o_atrialcc = {}; 
for an = 1:5
    data_an = udata{an};
    field_names = fieldnames(data_an);
    for ii = 1:length(field_names)
        varname = field_names{ii};
        cmdTxt = sprintf('clear %s;',varname);eval(cmdTxt);
        cmdTxt = sprintf('%s = data_an.%s;',varname,varname);eval(cmdTxt);
    end
    firing_rate = data_an.firing_rate;
    % firing_rate = data_an.ca_signal;
    for cn = 1:length(configurations)
        for ap = 1:2
            timecc = []; distcc = []; speedcc = []; FRcc = []; trialcc = [];
            otimecc = []; odistcc = []; ospeedcc = []; otrialcc = [];
            for tn = 1:10
                [an cn ap tn]
                % get time to complete trial
                if ap == 1
                    cmdTxt = sprintf('air_trials = air_trials_on .* %s;',configurations{cn}); eval(cmdTxt);
                else
                    cmdTxt = sprintf('air_trials = air_trials_off .* %s;',configurations{cn}); eval(cmdTxt);
                end
                idx = find(air_trials == tn & frf);
                tts = ts(idx)-ts(idx(1)); tds = ds(idx)-ds(idx(1)); tds(tds<0) = 0; 
                tsp = speed(idx);

                otimecc = [otimecc;tts']; odistcc = [odistcc;tds']; ospeedcc = [ospeedcc;tsp']; otrialcc = [otrialcc;tn*ones(size(tts'))];

                idx_fr = frf_n(idx);
                FR = firing_rate(:,idx_fr);
                bin_indices = discretize(tts,time_bins);
                time_binned = accumarray(bin_indices',tts,[],@mean); time_binned = time_binned - time_binned(1);
                dist_binned = accumarray(bin_indices',tds,[],@mean); dist_binned = dist_binned - dist_binned(1);
                speed_binned = accumarray(bin_indices',tsp,[],@mean);
                timecc = [timecc;time_binned]; distcc = [distcc;dist_binned]; speedcc = [speedcc;speed_binned];
                trialcc = [trialcc;tn*ones(size(time_binned))];
                FR_binned = [];
                parfor neuron_idx = 1:size(firing_rate,1)
                    % Calculate the binned firing rate for each neuron using the frame-based bin indices
                    FR_binned(neuron_idx, :) = accumarray(bin_indices', FR(neuron_idx, :),[], @mean, NaN)';
                end
                FR_binned = FR_binned';
                FRcc = [FRcc;FR_binned];

            end
            atimecc{an,cn,ap} = timecc; adistcc{an,cn,ap} = distcc; aspeedcc{an,cn,ap} = speedcc; aFRcc{an,cn,ap} = FRcc; atrialcc{an,cn,ap} = trialcc;
            o_atimecc{an,cn,ap} = otimecc; o_adistcc{an,cn,ap} = odistcc; o_aspeedcc{an,cn,ap} = ospeedcc; o_atrialcc{an,cn,ap} = otrialcc;
            n = 0;
        end
    end
end
n = 0;
%%
noofbMI = 10; nshuffles = 0;
MI_fun = @(x,y,noofbMI,nshuffles) calc_metric_MI(x,y,noofbMI,nshuffles);
PC_fun = @(x,y,nshuffles) calc_metric_PC(x,y,nshuffles);
variable_combs = {'FR_time','FR_dist','FR_speed'};
avar = []; Mvals={};
for an = 1:5
    anvar = [];
    for cn = 1:3
        for ap = 1:2
            timecc = atimecc{an,cn,ap}; distcc = adistcc{an,cn,ap}; speedcc = aspeedcc{an,cn,ap}; FRcc = aFRcc{an,cn,ap}; trialcc = atrialcc{an,cn,ap};
            for vn = 1:length(variable_combs)
                [an cn ap vn]
                var_name = variable_combs{vn};
                idx_us = strfind(var_name,'_');
                var1 = var_name(1:(idx_us-1)); var2 = var_name((idx_us+1):end);
                cmdTxt = sprintf('var1v = %scc;',var1);eval(cmdTxt); cmdTxt = sprintf('var2v = %scc;',var2);eval(cmdTxt);
                grouped_var2v = accumarray(trialcc, var2v, [], @(x) {x});
                num_neurons = size(var1v,2);
                PC_vals = []; MI_vals = []; PC_valscc = []; MI_valscc = [];
                clear metric_vals
                for nn = 1:num_neurons
                    thisFR = var1v(:,nn);
                    grouped_var1v = accumarray(trialcc, var1v(:,nn), [], @(x) {x});
                    thisvar = arrayfun(@(i) MI_fun(grouped_var1v{i}, grouped_var2v {i},noofbMI,nshuffles), 1:length(grouped_var1v));
                    MI_vals(:,nn) = thisvar;
                    thisvar = arrayfun(@(i) PC_fun(grouped_var1v{i}, grouped_var2v {i},nshuffles), 1:length(grouped_var1v));
                    PC_vals(:,nn) = thisvar';
                    MI_valscc(nn) = MI_fun(var1v(:,nn), var2v,noofbMI,nshuffles);
                    PC_valscc(nn) = PC_fun(var1v(:,nn), var2v,nshuffles);
                end
                metric_vals.PC_vals = PC_vals; metric_vals.PC_valscc = PC_valscc; metric_vals.MI_vals = MI_vals; metric_vals.MI_valscc = MI_valscc;
                Mvals{an,cn,ap,vn} = metric_vals;
            end
        end
    end
end
n = 0;
%% average speed and other speed characteristics
% for total time, distance use the variable_combs = {'time'} or distance
% and then use the max_fun. you may also change ap to 1 or 2 depending upon
% whether we want air on phase or air off phase

mean_fun = @(x) mean(x);
std_fun = @(x) std(x, 0);  % Standard deviation (unbiased)
cv_fun = @(x) std(x, 0) ./ mean(x);  % Coefficient of Variation (CV)
skew_fun = @(x) skewness(x);  % Skewness
kurt_fun = @(x) kurtosis(x);  % Kurtosis
max_fun = @(x) max(x); 
min_fun = @(x) min(x);
% latency_fun = @(x, t) t(find(x > 0, 1, 'first'));  
% rlatency_fun = @(x, t) t(find(x > 0, 1, 'last'));  

variable_combs = {'speed'};
% variable_combs = {'time'};
avar = [];
for an = 1:5
    anvar = [];
    for cn = 1:3
        for ap = 1:2
            timecc = o_atimecc{an,cn,ap}; distcc = o_adistcc{an,cn,ap}; speedcc = o_aspeedcc{an,cn,ap}; trialcc = o_atrialcc{an,cn,ap};
            for vn = 1:length(variable_combs)
                [an cn ap vn]
                var_name = variable_combs{vn};
                idx_us = strfind(var_name,'_');
                if length(variable_combs) == 1
                    cmdTxt = sprintf('var1v = %scc;',var_name);eval(cmdTxt);
                end
                thisvar = (accumarray(trialcc,var1v,[],max_fun))';
                anvar = [anvar thisvar];% outD.trial_metrics(:,idx)'];
                n = 0;
            end
        end
    end
    avar = [avar;anvar];
end
%
n = 0;
clc
[within,dvn,xlabels,awithinD] = make_within_table({'CN','AP','TN'},[3,2,10]);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
%%
% the following or the previous one to this will give error depending on
% the value of ap if it is only air on phase or air off phase
[within,dvn,xlabels,awithinD] = make_within_table({'CN','TN'},[3,10]);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
%% total time, distance and other things 

mean_fun = @(x) mean(x);
std_fun = @(x) std(x, 0);  % Standard deviation (unbiased)
cv_fun = @(x) std(x, 0) ./ mean(x);  % Coefficient of Variation (CV)
skew_fun = @(x) skewness(x);  % Skewness
kurt_fun = @(x) kurtosis(x);  % Kurtosis
max_fun = @(x) max(x); 
min_fun = @(x) min(x);
latency_fun = @(x, t) t(find(x > 0, 1, 'first'));  % movement_latency = accumarray(trialcc, speed, [], @(x) latency_fun(x, time));
rlatency_fun = @(x, t) t(find(x > 0, 1, 'last'));  % rest_latency = accumarray(trialcc, speed, [], @(x) latency_fun(x, time));
sp_thr = 0;
mon_fun = @(x) length(find_rising_edge(x > sp_thr,0.1,500));
moff_fun = @(x) length(find_falling_edge(x > sp_thr,-0.1,500));
latency_fun1 = @(x, t) t(find_first_rising_edge(x > sp_thr,0.1,0));  % movement_latency = accumarray(trialcc, speed, [], @(x) latency_fun(x, time));
rlatency_fun1 = @(x, t) t(find_first_falling_edge(x > sp_thr,-0.1,0));

variable_combs = {'speed'};
avar = [];
for an = 1:5
    anvar = [];
    for cn = 1:3
        for ap = 2
            timecc = o_atimecc{an,cn,ap}; distcc = o_adistcc{an,cn,ap}; speedcc = o_aspeedcc{an,cn,ap}; trialcc = o_atrialcc{an,cn,ap};
            for vn = 1:length(variable_combs)
                [an cn ap vn]
                var_name = variable_combs{vn};
                idx_us = strfind(var_name,'_');
                if length(variable_combs) == 1
                    cmdTxt = sprintf('var1v = %scc;',var_name);eval(cmdTxt);
                end
                grouped_var1v = accumarray(trialcc, var1v, [], @(x) {x}); grouped_timecc = accumarray(trialcc, timecc, [], @(x) {x});
                thisvar = arrayfun(@(i) rlatency_fun1(grouped_var1v{i}, grouped_timecc{i}), 1:length(grouped_var1v));
                % thisvar = (accumarray(trialcc,var1v,[],moff_fun))';
                anvar = [anvar thisvar];
                n = 0;
            end
        end
    end
    avar = [avar;anvar];
end
descriptiveStatistics(avar(:));
%%
n = 0;
clc
[within,dvn,xlabels,awithinD] = make_within_table({'CN','AP','TN'},[3,2,10]);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
%%
clc
raR = RMA_bonferroni(ra,1);
%%
clc
[within,dvn,xlabels,awithinD] = make_within_table({'CN','TN'},[3,10]);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)

%%
noofbMI = 10; nshuffles = 0;
MI_fun = @(x,y,noofbMI,nshuffles) calc_metric_MI(x,y,noofbMI,nshuffles);
PC_fun = @(x,y,nshuffles) calc_metric_PC(x,y,nshuffles);
variable_combs = {'time_dist','time_speed','dist_speed'};
variable_combs = {'time_speed','dist_speed'};
avar = [];
for an = 1:5
    anvar = [];
    for cn = 1:3
        for ap = 1:2
            timecc = atimecc{an,cn,ap}; distcc = adistcc{an,cn,ap}; speedcc = aspeedcc{an,cn,ap}; FRcc = aFRcc{an,cn,ap}; trialcc = atrialcc{an,cn,ap};
            for vn = 1:length(variable_combs)
                [an cn ap vn]
                var_name = variable_combs{vn};
                idx_us = strfind(var_name,'_');
                var1 = var_name(1:(idx_us-1)); var2 = var_name((idx_us+1):end);
                cmdTxt = sprintf('var1v = %scc;',var1);eval(cmdTxt); cmdTxt = sprintf('var2v = %scc;',var2);eval(cmdTxt);
                grouped_var1v = accumarray(trialcc, var1v, [], @(x) {x}); grouped_var2v = accumarray(trialcc, var2v, [], @(x) {x});
                thisvar = arrayfun(@(i) MI_fun(grouped_var1v{i}, grouped_var2v {i},noofbMI,nshuffles), 1:length(grouped_var1v));
                % thisvar = arrayfun(@(i) PC_fun(grouped_var1v{i}, grouped_var2v {i},nshuffles), 1:length(grouped_var1v));

                % thisvar = MI_fun(var1v, var2v,noofbMI,nshuffles);
                % thisvar = PC_fun(var1v, var2v,nshuffles);

                anvar = [anvar thisvar];% outD.trial_metrics(:,idx)'];
                n = 0;
            end
        end
    end
    avar = [avar;anvar];
end
descriptiveStatistics(avar(:));
%%
n = 0;
clc
fac_names = {'CN','TN'}; fac_levels = [3,10];
fac_names = {'CN','PT','TN'}; fac_levels = [3,3,10];
% fac_names = {'CN','AP','TN'}; fac_levels = [3,2,10];
fac_names = {'CN','AP','PT','TN'}; fac_levels = [3,2,2,10];
% fac_names = {'CN','AP','PT'}; fac_levels = [3,2,3];
% fac_names = {'CN','PT'}; fac_levels = [3,3];
[within,dvn,xlabels,awithinD] = make_within_table(fac_names,fac_levels);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
%%
ra1 = ra;
ra2 = ra;
%%
clc
raR = RMA_bonferroni(ra,2);
%%
clc
raR = RMA_subset(ra,'AP');
%%
raRR = RMA_bonferroni(raR{1},'PT');
% [xdata,mVar,semVar,combs,p,h,nB] = RMA_get_multcompare(raR{1},{'TN','hsd',0.05},[1 1]);
%%
raRB = RMA_bonferroni(raR.ras{2},3);
[xdata,mVar,semVar,combs,p,h,nB] = RMA_get_multcompare(raRB{1},{'TN','hsd',0.05},[1 1]);
%%
n = 0;
clc
[within,dvn,xlabels,awithinD] = make_within_table({'CN','TN','PT'},[3,10,3]);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
%% speed time and dist glm
avar = [];
for an = 1:5
    anvar = [];
    for cn = 1:3
        for ap = 1:2
            timecc = o_atimecc{an,cn,ap}; distcc = o_adistcc{an,cn,ap}; speedcc = o_aspeedcc{an,cn,ap}; trialcc = o_atrialcc{an,cn,ap};
            Y = speedcc; % Firing rate for current neuron (timepoints x 1)
            X = [timecc,distcc]; % Predictors (time, distance, speed)
            
            % Fit the GLM (here we use a linear model, but you can modify it for other GLM families)
            mdl = fitglm(X, Y, 'linear');
            aglm{an,cn,ap} = mdl;
            anvar = [anvar mdl.Coefficients.pValue(1) mdl.Coefficients.pValue(2)];
        end
    end
    avar = [avar;anvar];
end
clc
[within,dvn,xlabels,awithinD] = make_within_table({'CN','AP','DT'},[3,2,2]);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
%%
aglm = {};
for an = 1:5
    for cn = 1:3
        for ap = 1:2
            [an cn ap]
            timecc = atimecc{an,cn,ap}; distcc = adistcc{an,cn,ap}; speedcc = aspeedcc{an,cn,ap}; FRcc = aFRcc{an,cn,ap};
            % Assuming `timecc`, `distcc`, `speedcc`, and `FRcc` are already defined

            numNeurons = size(FRcc, 2); % Number of neurons
            results = cell(numNeurons, 1); % Store GLM results for each neuron
            parfor neuron_idx = 1:numNeurons
                % Prepare the data for the GLM
                Y = FRcc(:, neuron_idx); % Firing rate for current neuron (timepoints x 1)
                X = [timecc, distcc, speedcc]; % Predictors (time, distance, speed)
                
                % Fit the GLM (here we use a linear model, but you can modify it for other GLM families)
                mdl = fitglm(X, Y, 'linear'); % You can change 'linear' to other distributions as needed
                
                % Store the model results for the current neuron
                results{neuron_idx} = mdl;
            end
            aglm{an,cn,ap} = results;
        end
    end
end
%%
% Assuming aglm{an, cn, ap} contains the GLM results for each animal, configuration, and phase
numAnimals = size(aglm, 1); % Number of animals
numConfigurations = size(aglm, 2); % Number of configurations
numAirPhases = size(aglm, 3); % Number of air-phases
% numNeurons = size(FRcc{1,1,1}, 2); % Number of neurons (assuming FRcc has neurons in columns)

% Initialize a cell to store the tuning results (for each animal, configuration, and phase)
tuningResults = cell(numAnimals, numConfigurations, numAirPhases);

% Loop through each animal, configuration, and air-phase
for an = 1:numAnimals
    for cn = 1:numConfigurations
        for ap = 1:numAirPhases
            % Get the GLM results for this animal, configuration, and air-phase
            results = aglm{an, cn, ap};
            numNeurons = length(results);
            % Initialize an array to store the tuning category for each neuron
            neuronTuning = cell(numNeurons, 1);
            
            % Loop through each neuron
            for neuron_idx = 1:numNeurons
                % Get the GLM model for this neuron
                mdl = results{neuron_idx};
                
                % Extract the coefficients for time, distance, and speed
                coefficients = mdl.Coefficients.Estimate(2:4); % Time, Distance, Speed (assuming first coefficient is intercept)
                pvals = mdl.Coefficients.pValue(2:4); % Time, Distance, Speed (assuming first coefficient is intercept)
                sigpval = pvals' < 0.05;
                sigpvaldec = (sigpval(1) * 4 + sigpval(2) * 2 + sigpval(3) * 1);
                tuningTypes = {'Untuned','Speed','Distance','Dist-Speed','Time','Time-Speed','Time-Dist','Time-Dist-Speed'};
                thistuning = tuningTypes{sigpvaldec + 1};
                neuronTuning{neuron_idx} = thistuning;
            end
            % Store the tuning results for this animal, configuration, and air-phase
            tuningResults{an, cn, ap} = neuronTuning;
        end
    end
end
%%
% Initialize arrays to store percentages (or fractions) for time, distance, and speed tuning
timeTunedPercent = zeros(numAnimals, numConfigurations, numAirPhases);
distanceTunedPercent = zeros(numAnimals, numConfigurations, numAirPhases);
speedTunedPercent = zeros(numAnimals, numConfigurations, numAirPhases);

% Loop through the tuning results and calculate the percentages (or fractions)
avar = [];
for an = 1:numAnimals
    anvar = [];
    for cn = 1:numConfigurations
        for ap = 1:numAirPhases
            % Get the tuning results for this animal, configuration, and phase
            neuronTuning = tuningResults{an, cn, ap};
            numNeurons = length(neuronTuning);
            % tuningTypes = {'Untuned','Speed','Distance','Dist-Speed','Time','Time-Speed','Time-Dist','Time-Dist-Speed'};
            numtuned = [];
            for tii = 1:length(tuningTypes)
                fractuned(1,tii) = sum(strcmp(neuronTuning, tuningTypes{tii}))/numNeurons;
            end
            anvar = [anvar fractuned([2 3 5])];
        end
    end
    avar = [avar;anvar];
end

%%
n = 0;
clc
[within,dvn,xlabels,awithinD] = make_within_table({'CN','AP','TT'},[3,2,3]);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
%%
clc
raR = RMA_bonferroni(ra,3);


%%
noofbMI = 10; nshuffles = 0;
MI_fun = @(x,y,noofbMI,nshuffles) calc_metric_MI(x,y,noofbMI,nshuffles);
PC_fun = @(x,y,nshuffles) calc_metric_PC(x,y,nshuffles);
variable_combs = {'FR_time','FR_dist','FR_speed'};
avar = []; FRvals={};
for an = 1:5
    anvar = [];
    for cn = 1:3
        for ap = 1:2
            timecc = atimecc{an,cn,ap}; distcc = adistcc{an,cn,ap}; speedcc = aspeedcc{an,cn,ap}; FRcc = aFRcc{an,cn,ap}; trialcc = atrialcc{an,cn,ap};
            for vn = 1:length(variable_combs)
                [an cn ap vn]
                var_name = variable_combs{vn};
                idx_us = strfind(var_name,'_');
                var1 = var_name(1:(idx_us-1)); var2 = var_name((idx_us+1):end);
                cmdTxt = sprintf('var1v = %scc;',var1);eval(cmdTxt); cmdTxt = sprintf('var2v = %scc;',var2);eval(cmdTxt);
                grouped_var2v = accumarray(trialcc, var2v, [], @(x) {x});
                num_neurons = size(var1v,2);
                metric_vals = [];
                for nn = 1:num_neurons
                    thisFR = var1v(:,nn);
                    grouped_var1v = accumarray(trialcc, var1v(:,nn), [], @(x) {x});
                    thisvar = arrayfun(@(i) MI_fun(grouped_var1v{i}, grouped_var2v {i},noofbMI,nshuffles), 1:length(grouped_var1v));
                    thisvar = arrayfun(@(i) PC_fun(grouped_var1v{i}, grouped_var2v {i},nshuffles), 1:length(grouped_var1v));
                    metric_vals(:,nn) = thisvar';
                    % thisvar = MI_fun(var1v, var2v,noofbMI,nshuffles);
                    % thisvar = PC_fun(var1v, var2v,nshuffles);
                end
                FRvals{an,cn,ap,vn} = metric_vals;
                thisvar = nanmean(metric_vals,2);
                anvar = [anvar thisvar'];% outD.trial_metrics(:,idx)'];
                n = 0;
            end
        end
    end
    avar = [avar;anvar];
end
descriptiveStatistics(avar(:));
%%
n = 0;
clc
[within,dvn,xlabels,awithinD] = make_within_table({'CN','AP','TN','PT'},[3,2,10,3]);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)
%%
clc
raR = RMA_bonferroni(ra,4);
%%
n = 0;
clc
[within,dvn,xlabels,awithinD] = make_within_table({'CN','AP','PT'},[3,2,3]);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
print_for_manuscript(ra)