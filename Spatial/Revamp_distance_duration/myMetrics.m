function out = myMetrics(x,y,params)

% params = {'no_of_bins_for_MI',10,'no_of_shuffles_for_norm',500,'animal_info',ei,'overwrite_processing',[0 0]};

param = 'no_of_bins_for_MI'; idx = find(strcmp(params,param)); cmdTxt = sprintf('%s = params{idx+1};',param);eval(cmdTxt);
param = 'no_of_shuffles_for_norm'; idx = find(strcmp(params,param)); cmdTxt = sprintf('%s = params{idx+1};',param);eval(cmdTxt);
param = 'animal_info'; idx = find(strcmp(params,param)); cmdTxt = sprintf('%s = params{idx+1};',param);eval(cmdTxt);
param = 'overwrite_processing'; idx = find(strcmp(params,param)); cmdTxt = sprintf('%s = params{idx+1};',param);eval(cmdTxt);
param = 'air_phase'; idx = find(strcmp(params,param)); cmdTxt = sprintf('%s = params{idx+1};',param);eval(cmdTxt);
param = 'configuration'; idx = find(strcmp(params,param)); cmdTxt = sprintf('%s = params{idx+1};',param);eval(cmdTxt);
param = 'variables'; idx = find(strcmp(params,param)); cmdTxt = sprintf('%s = params{idx+1};',param);eval(cmdTxt);
param = 'bin_type'; idx = find(strcmp(params,param)); cmdTxt = sprintf('%s = params{idx+1};',param);eval(cmdTxt);
param = 'trial_type'; idx = find(strcmp(params,param)); cmdTxt = sprintf('%s = params{idx+1};',param);eval(cmdTxt);


file_name = sprintf('metrics_%s_%s_%s_%s_%s.mat',variables,bin_type,trial_type,configuration,air_phase);
folder = fullfile(animal_info.matlab_folder,"unlv"); 
if ~exist(folder,'dir')
    mkdir(folder)
end
file_name = fullfile(folder,file_name);
if overwrite_processing(1) < 0
    if exist(file_name)
        disp(file_name);
        delete(file_name);
        out.MI = NaN; out.PC = NaN;
        return;
    end
end
if overwrite_processing(1) >= 1 || ~exist(file_name)
    % disp(file_name);
    MI = NaN(size(x,2),3); PC = MI; nobmi = no_of_bins_for_MI; noshffl = no_of_shuffles_for_norm;
    parfor nn = 1:size(x,2)
        [MI(nn,:),PC(nn,:)] = calc_metrics(x(:,nn),y,nobmi,noshffl);
    end
    out.MI = MI; out.PC = PC;
    if overwrite_processing < 2 
        save(file_name,'-struct','out','-v7.3');
    end
else
    out = load(file_name);
end


% rng(3,'twister');
% [output ~] = info_metrics_S(x, y, no_of_bins_for_MI, [], no_of_shuffles_for_norm);
n = 0;


