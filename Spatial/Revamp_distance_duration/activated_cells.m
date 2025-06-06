function activated_cells_trial

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
% 
% an = 1;
% tudata = udataT{an}
% field_names = fieldnames(tudata);
% for ii = 1:length(field_names)
%     varname = field_names{ii};
%     cmdTxt = sprintf('%s = tudata.%s;',varname,varname);
%     eval(cmdTxt);
% end
%%
% Define parameters
configs = {'C2', 'C3', 'C4', 'C5', 'C7'};
num_animals = length(udataT); % Number of animals
num_configs = length(configs); % Number of configurations
num_trials = 10; % Trials per configuration
num_before_after = 2; % Before/After (1 = Before, 2 = After)
num_time_windows = 4; % 4 time windows per side
num_time_points = 3; % 3 bins per window
bin_size = mean(diff(udataT{1}.ts)); % Bin size from timestamps
num_bins_per_window = round(1 / bin_size); % 1-second window = num_bins_per_window

% Initialize multidimensional percent_active
percent_active = NaN(num_animals, num_configs, num_trials, num_before_after, num_time_windows, num_time_points);

% Loop through each animal
for a = 1:num_animals
    % Get animal-specific data
    firing_rate = udataT{a}.firing_rate;
    air = udataT{a}.air;
    ts = udataT{a}.ts;
    
    % Initialize structure for air-on events
    air_events = struct();
    
    % Loop through configurations
    for c = 1:num_configs
        config = configs{c};
        
        % Get configuration-specific indices
        config_idx = udataT{a}.(config) == 1;
        
        % Find air-onset indices within this configuration
        air_onsets = find(diff(air) == 1 & config_idx(2:end));
        
        % Store air-onset windows (4 before, 4 after, each 1 second wide)
        air_events.(config).onset_windows = arrayfun(@(x) ...
            [  (x-4*num_bins_per_window:x-3*num_bins_per_window-1); ...
               (x-3*num_bins_per_window:x-2*num_bins_per_window-1); ...
               (x-2*num_bins_per_window:x-1*num_bins_per_window-1); ...
               (x-1*num_bins_per_window:x+0*num_bins_per_window-1); ... % 4 BEFORE
               (x+0*num_bins_per_window:x+1*num_bins_per_window-1); ...
               (x+1*num_bins_per_window:x+2*num_bins_per_window-1); ...
               (x+2*num_bins_per_window:x+3*num_bins_per_window-1); ...
               (x+3*num_bins_per_window:x+4*num_bins_per_window-1) ], ... % 4 AFTER
               air_onsets, 'UniformOutput', false);
        
        % Compute percentage of active cells for each air-onset event
        for t = 1:length(air_onsets)
            for b = 1:2 % 1 = Before, 2 = After
                for w = 1:num_time_windows
                    % Determine window index
                    window_idx = (b - 1) * num_time_windows + w;
                    
                    % Get bin indices for the window
                    bin_idx = air_events.(config).onset_windows{t}(window_idx, :);
                    
                    % Calculate percentage of active cells for each time point
                    for p = 1:num_time_points
                        active_cells = firing_rate(:, bin_idx(p)) > 0; % Active = firing > 0
                        percent_active(a, c, t, b, w, p) = sum(active_cells) / size(firing_rate, 1) * 100;
                    end
                end
            end
        end
    end
end

% Display summary of percent_active
disp('Multidimensional percent_active created:');
disp(size(percent_active)); % Should be [num_animals, num_configs, num_trials, 2, 4, 3]
%%
% Initialize parameters
num_animals = size(percent_active, 1);
num_configs = size(percent_active, 2);
num_trials = size(percent_active, 3);
num_before_after = size(percent_active, 4);
num_time_windows = size(percent_active, 5);
num_time_points = size(percent_active, 6);

% Total columns for ANOVA matrix
total_columns = num_configs * num_trials * num_before_after * num_time_windows * num_time_points;

% Initialize ANOVA-style matrix
anova_matrix = NaN(num_animals, total_columns);

% Column index for filling the matrix
col_idx = 1;

% Loop through configurations
for c = 1:num_configs
    % Loop through trials
    for t = 1:num_trials
        % Loop through Before/After
        for b = 1:num_before_after
            % Loop through time windows
            for w = 1:num_time_windows
                % Loop through time points
                for p = 1:num_time_points
                    % Extract data for all animals
                    anova_matrix(:, col_idx) = percent_active(:, c, t, b, w, p);
                    
                    % Increment column index
                    col_idx = col_idx + 1;
                end
            end
        end
    end
end

% Display summary
disp('ANOVA-style matrix created:');
disp(size(anova_matrix)); % Should be [5, 1200]


%%
[within,dvn,xlabels,awithinD] = make_within_table({'Conf','TR','BA','TW','TP'},[5,10,2,4,3]);
dataT = make_between_table({anova_matrix},dvn);
ra = RMA(dataT,within,{0.05,{''}});
ra.ranova
print_for_manuscript(ra)
%%
% redF = [1,2]; redV = {[1,2],1};
% [dataTR,withinR] = reduce_within_between(dataT,within,redF,redV);
% raR = RMA(dataTR,withinR,{0.025,{'hsd'}});
% raR.ranova
% print_for_manuscript(raR)
% 
% Average over Time Points (TP) in the percent_active variable
% Input dimensions of percent_active: [Animals, Configurations, Trials, BA, TW, TP]

% Average along the TP dimension (6th dimension)
percent_active_avg = mean(percent_active, 6, 'omitnan'); % Resulting size: [Animals, Configurations, Trials, BA, TW]

% Display size of the averaged percent_active matrix
disp('Dimensions of percent_active after averaging over TP:');
disp(size(percent_active_avg)); % Should be [Animals, Configurations, Trials, BA, TW]

% Prepare an ANOVA-style matrix
% Initialize parameters
num_animals = size(percent_active_avg, 1);
num_configs = size(percent_active_avg, 2);
num_trials = size(percent_active_avg, 3);
num_before_after = size(percent_active_avg, 4);
num_time_windows = size(percent_active_avg, 5);

% Total columns for ANOVA matrix
total_columns = num_configs * num_trials * num_before_after * num_time_windows;

% Initialize ANOVA-style matrix
anova_matrix_avg = NaN(num_animals, total_columns);

% Column index for filling the matrix
col_idx = 1;

% Loop through configurations
for c = 1:num_configs
    % Loop through trials
    for t = 1:num_trials
        % Loop through Before/After
        for b = 1:num_before_after
            % Loop through Time Windows
            for w = 1:num_time_windows
                % Extract data for all animals
                anova_matrix_avg(:, col_idx) = percent_active_avg(:, c, t, b, w);
                
                % Increment column index
                col_idx = col_idx + 1;
            end
        end
    end
end

% Display summary of the ANOVA matrix
disp('ANOVA-style matrix after averaging over TP created:');
disp(size(anova_matrix_avg)); % Should be [num_animals, total_columns]
%%
[within,dvn,xlabels,awithinD] = make_within_table({'Conf','TR','BA','TW'},[5,10,2,4]);
dataT = make_between_table({anova_matrix_avg},dvn);
ra = RMA(dataT,within,{0.05,{''}});
ra.ranova
print_for_manuscript(ra)