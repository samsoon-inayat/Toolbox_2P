%%
for ani = 1:5
    tudata = udataT{ani};
    mibnb = get_unlv_analysis_multiFunc(tudata,'MI_BNB');
end
%%
[within,dvn,xlabels] = make_within_table({'SC'},[7]);
dataT = make_between_table({mAM},dvn);
ra = RMA(dataT,within,{0.05,{'hsd'}});
%     ra.ranova
print_for_manuscript(ra);

magfac = mData.magfac;
ff = makeFigureRowsCols(107,[3 5 2 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.3],...
    'widthHeightAdjustment',[10 -450]);
MY = 40; ysp = 9; mY = 0; ystf = 7; ysigf = 1;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.33*magfac; widths = [1.1 0.35 2.85 1]*magfac+0.061; gap = 0.09105*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 1.5];

tcolors = repmat(mData.colors(5:end),1,6);
[hbs,xdata] = view_results_rmanova(ff.h_axes(1,1),ra,{'SC','hsd',0.05},xs_gaps,tcolors,[mY MY ysp ystf ysigf],mData);
xticklabels = {'Comb1','Comb2','Comb3','Comb4','Comb5','Comb6','Comb7'}; set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
ylabel('Mean Firing Rate');

% tcolors = repmat(mData.colors(1:2),1,1);
% [hbs,xdata] = view_results_rmanova(ff.h_axes(1,2),raCCCC,{'Ph','hsd',0.05},xs_gaps,tcolors,[mY MY ysp ystf ysigf],mData);
% xticklabels = {'AOn','AOff'};  set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
ht = axes_title(ff,{1},{'All Cells'},axes_title_shifts_line,axes_title_shifts_text,'no');
set(ht,'FontWeight','Bold');


save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);
%%

ani = 3;
tudata = udataT{ani};
field_names = fieldnames(tudata);
for ii = 1:length(field_names)
    varname = field_names{ii};
    cmdTxt = sprintf('%s = tudata.%s;',varname,varname);
    eval(cmdTxt);
end

whos firing_rate speed air belt light ts ds bnb C1 C2 C3 C4 C5 C6 C7 animal_motion frf

xs = ts;
figure(1000);clf;
plot(xs,light,'r');hold on;
plot(xs,air,'b');
plot(xs,bnb*1.2,'g');
plot(xs,C1*1.1,'k');
plot(xs,C2*1.1,'k');
plot(xs,C3*1.1,'k');
plot(xs,C4*1.1,'k');
plot(xs,C5*1.1,'k');
plot(xs,C6*1.1,'k');
plot(xs,C7*1.1,'k');
% plot(frf * 0.5,'c')

%%
% Find the peak (max) firing rate for each neuron
peak_values = max(firing_rate, [], 2);  % Find max along the second dimension (time)

% Replicate the peak values across the time dimension for element-wise division
peak_values_rep = repmat(peak_values, 1, size(firing_rate, 2));

% Normalize firing rate by the peak values
firing_rate_normalized = firing_rate ./ peak_values_rep;

% Handle cases where peak values are zero to avoid division by zero
firing_rate_normalized(isnan(firing_rate_normalized)) = 0;


% Sort neurons by peak value
[~, sort_idx] = sort(peak_values, 'descend');  % Get indices of neurons sorted by peak values in descending order

% Select the top 100 neurons
top_100_idx = sort_idx(1:100);  % Indices of the top 100 neurons

% Extract the firing rates of the top 100 neurons
firing_rate_top_100 = firing_rate_normalized(top_100_idx, :);


% Reorder the firing_rate_normalized matrix based on the sorted neuron indices
firing_rate_normalized_sorted = firing_rate_normalized(sort_idx, :);


figure(1000);clf
imagesc(firing_rate_top_100);



%%
% Define binary labels for each combination of Brake, Air, and Light, including no stimulus
comb_Brake_Light = (bnb == 1) & (air == 0) & (light == 1);  % Brake + Light
comb_Brake_Air = (bnb == 1) & (air == 1) & (light == 0);    % Brake + Air
comb_Air = (bnb == 0) & (air == 1) & (light == 0);          % Air only
comb_Air_Light = (bnb == 0) & (air == 1) & (light == 1);    % Air + Light
comb_No_Stimulus = (bnb == 0) & (air == 0) & (light == 0);  % No stimulus

% Define a threshold to identify locomotion (e.g., speed > 1 cm/s)
locomotion_threshold = 1;  % Define based on data (adjust as needed)

% Define binary labels for locomotion and non-locomotion periods
locomotion = speed > locomotion_threshold;  % Locomotion periods (speed above threshold)
non_locomotion = speed <= locomotion_threshold;  % Non-locomotion periods (speed below threshold)

% Initialize result matrix for mean firing rates (now with 7 conditions)
mean_firing_rates_comb = zeros(size(firing_rate, 1), 7);  % Include 7 conditions

% Calculate mean firing rate for each condition, including locomotion
mean_firing_rates_comb(:, 1) = mean(firing_rate(:, comb_Brake_Light), 2);   % Brake + Light
mean_firing_rates_comb(:, 2) = mean(firing_rate(:, comb_Brake_Air), 2);     % Brake + Air
mean_firing_rates_comb(:, 3) = mean(firing_rate(:, comb_Air), 2);           % Air only
mean_firing_rates_comb(:, 4) = mean(firing_rate(:, comb_Air_Light), 2);     % Air + Light
mean_firing_rates_comb(:, 5) = mean(firing_rate(:, comb_No_Stimulus), 2);   % No stimulus
mean_firing_rates_comb(:, 6) = mean(firing_rate(:, locomotion), 2);         % Locomotion
mean_firing_rates_comb(:, 7) = mean(firing_rate(:, non_locomotion), 2);     % Non-locomotion

% Plot mean firing rates for all conditions, including locomotion
figure;
bar(mean_firing_rates_comb);
xlabel('Neuron');
ylabel('Mean Firing Rate');
legend({'Brake + Light', 'Brake + Air', 'Air', 'Air + Light', 'No Stimulus', 'Locomotion', 'Non-Locomotion'}, 'Location', 'northeastoutside');
title('Mean Firing Rate for Different Stimulus Combinations and Locomotion');

%%
% Example: Perform peri-event analysis for each neuron in a specific window
pre_event_window = 50;  % bins before the event
post_event_window = 100; % bins after the event

% Align firing rates to a specific stimulus onset (e.g., Air Onset)
aligned_firing_rates_air = zeros(length(find(air == 1)), pre_event_window + post_event_window, size(firing_rate, 1));
air_onset_times = find(air == 1);

for i = 1:length(air_onset_times)
    if air_onset_times(i) > pre_event_window && air_onset_times(i) + post_event_window <= length(ts)
        aligned_firing_rates_air(i, :, :) = firing_rate(:, air_onset_times(i) - pre_event_window : air_onset_times(i) + post_event_window - 1)';
    end
end

% Compute average firing rate across all air-on events for each neuron
mean_aligned_firing_rate_air = squeeze(mean(aligned_firing_rates_air, 1));

% Plot the peri-event activity for a neuron of interest
neuron_idx = 1;  % Example neuron
figure;
plot(-pre_event_window:post_event_window-1, mean_aligned_firing_rate_air(:, neuron_idx));
xlabel('Time from Air On (bins)');
ylabel('Average Firing Rate');
title(['Neuron ', num2str(neuron_idx), ' Peri-Event Average (Air)']);

%%
% Collect speeds for each combination into one array
all_speeds = [speed(comb_Brake_Light)', speed(comb_Brake_Air)', speed(comb_Air)', speed(comb_Air_Light)', speed(comb_No_Stimulus)'];

% Create a grouping variable for the combinations
group_labels = [repmat({'Brake + Light'}, sum(comb_Brake_Light), 1); ...
                repmat({'Brake + Air'}, sum(comb_Brake_Air), 1); ...
                repmat({'Air'}, sum(comb_Air), 1); ...
                repmat({'Air + Light'}, sum(comb_Air_Light), 1); ...
                repmat({'No Stimulus'}, sum(comb_No_Stimulus), 1)];

% Perform one-way ANOVA
[p, tbl, stats] = anova1(all_speeds, group_labels);

% Display post-hoc comparisons
multcompare(stats);

%%
% Define binary labels for Brake and No-Brake conditions
brake_conditions = {C1, C2, C6, C7};      % Configurations under Brake condition
no_brake_conditions = {C3, C4, C5};       % Configurations under No-Brake condition

% Initialize result arrays for mean speeds across Brake and No-Brake configurations
mean_speed_brake = zeros(1, 4);   % For C1, C2, C6, C7 (Brake)
mean_speed_no_brake = zeros(1, 3); % For C3, C4, C5 (No-Brake)

% Calculate mean speed for each configuration under Brake
mean_speed_brake(1) = mean(speed(C1));  % Mean speed for C1 (Brake)
mean_speed_brake(2) = mean(speed(C2));  % Mean speed for C2 (Brake)
mean_speed_brake(3) = mean(speed(C6));  % Mean speed for C6 (Brake)
mean_speed_brake(4) = mean(speed(C7));  % Mean speed for C7 (Brake)

% Calculate mean speed for each configuration under No-Brake
mean_speed_no_brake(1) = mean(speed(C3));  % Mean speed for C3 (No-Brake)
mean_speed_no_brake(2) = mean(speed(C4));  % Mean speed for C4 (No-Brake)
mean_speed_no_brake(3) = mean(speed(C5));  % Mean speed for C5 (No-Brake)

% Plot the mean speed for Brake and No-Brake conditions
figure;

% Plot Brake conditions
subplot(1,2,1);
bar(mean_speed_brake);
xlabel('Brake Configurations');
ylabel('Mean Speed (cm/s)');
ylim([0 10]);
% xticklabels({'C1', 'C2', 'C6', 'C7'});
title('Mean Speed for Brake Conditions');

% Plot No-Brake conditions
subplot(1,2,2);
bar(mean_speed_no_brake);
xlabel('No-Brake Configurations');
ylabel('Mean Speed (cm/s)');
ylim([0 10]);
% xticklabels({'C3', 'C4', 'C5'});
title('Mean Speed for No-Brake Conditions');

%%
msaonoff = [];
for ani = 1:5
    out = get_unlv_analysis_multiFunc(udata{ani},'C345_AON_AOFF_avgspeed');
    mean_speed_air_on = out.msaon;
    mean_speed_air_off = out.msaoff
    msaonoff = [msaonoff;[mean_speed_air_on mean_speed_air_off]];

    
    % Plot the mean speeds for Air On and Air Off periods
    figure(1000);clf
    
    
    % Plot for Air On periods
    subplot(1,2,1);
    bar(mean_speed_air_on);
    xlabel('No-Brake Configurations');
    ylabel('Mean Speed (cm/s)');ylim([0 25]);
    set(gca, 'XTickLabel', {'C3', 'C4', 'C5'});  % Alternative method for xticklabels
    title('Mean Speed During Air On (No-Brake)');
    
    % Plot for Air Off periods
    subplot(1,2,2);
    bar(mean_speed_air_off);
    xlabel('No-Brake Configurations'); ylim([0 25]);
    ylabel('Mean Speed (cm/s)');
    set(gca, 'XTickLabel', {'C3', 'C4', 'C5'});  % Alternative method for xticklabels
    title('Mean Speed During Air Off (No-Brake)');
    pause(1)
end

%%
[within,dvn,xlabels] = make_within_table({'APhase','Conf'},[2,3]);
dataT = make_between_table({msaonoff},dvn);
ra = RMA(dataT,within,{0.05,{'hsd'}});
%     ra.ranova
print_for_manuscript(ra);
%%
magfac = mData.magfac;
ff = makeFigureRowsCols(107,[3 5 2 1],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.3],...
    'widthHeightAdjustment',[10 -450]);
MY = 40; ysp = 9; mY = 0; ystf = 7; ysigf = 1;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.33*magfac; widths = [1.1 0.35 2.85 1]*magfac+0.061; gap = 0.09105*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 1.5];

tcolors = repmat(mData.colors(5:end),1,6);
[hbs,xdata] = view_results_rmanova(ff.h_axes(1,1),ra,{'APhase:Conf','hsd',0.05},xs_gaps,tcolors,[mY MY ysp ystf ysigf],mData);
xticklabels = {'Comb1','Comb2','Comb3','Comb4','Comb5','Comb6','Comb7'}; set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
ylabel('Mean Firing Rate');

% tcolors = repmat(mData.colors(1:2),1,1);
% [hbs,xdata] = view_results_rmanova(ff.h_axes(1,2),raCCCC,{'Ph','hsd',0.05},xs_gaps,tcolors,[mY MY ysp ystf ysigf],mData);
% xticklabels = {'AOn','AOff'};  set(gca,'xtick',xdata,'xticklabels',xticklabels); xtickangle(30);
ht = axes_title(ff,{1},{'All Cells'},axes_title_shifts_line,axes_title_shifts_text,'no');
set(ht,'FontWeight','Bold');


save_pdf(ff.hf,mData.pdf_folder,'bar_graph.pdf',600);

%%
% Define Air On and Air Off periods
air_on = (air == 1);    % Air On periods
air_off = (air == 0);   % Air Off periods

% Define the No-Brake configurations (C3, C4, C5)
no_brake_C3 = C3;  
no_brake_C4 = C4;  
no_brake_C5 = C5;  

% --- Time Binning ---

% Define time bins (every 0.3 seconds)
time_bin_size = 0.3;  % 0.3 seconds
time_bins = 0:time_bin_size:max(ts);  % Assuming max time in ts defines the trial length

% Initialize firing rate matrix for time bins (rows: neurons, columns: bins)
firing_rate_time_on = zeros(size(firing_rate, 1), length(time_bins) - 1);  % Air On
firing_rate_time_off = zeros(size(firing_rate, 1), length(time_bins) - 1); % Air Off

% Initialize logical arrays for air on/off within each time bin
air_on_time_bins = false(1, length(time_bins) - 1);
air_off_time_bins = false(1, length(time_bins) - 1);

% Calculate firing rates and logical variables for Air On and Air Off periods (time tuning)
for i = 1:length(time_bins)-1
    % Air On period (variable running time)
    bin_idx_on = (ts >= time_bins(i)) & (ts < time_bins(i+1)) & air_on;
    firing_rate_time_on(:, i) = mean(firing_rate(:, bin_idx_on), 2);
    air_on_time_bins(i) = any(bin_idx_on);  % Update logical array
    
    % Air Off period (fixed 15-second interval)
    bin_idx_off = (ts >= time_bins(i)) & (ts < time_bins(i+1)) & air_off;
    firing_rate_time_off(:, i) = mean(firing_rate(:, bin_idx_off), 2);
    air_off_time_bins(i) = any(bin_idx_off);  % Update logical array
end

% --- Distance Binning ---

% Define the distance bins (every 3 cm)
distance_bin_size = 3;  % 3 cm
distance_bins = 0:distance_bin_size:max(ds);

% Initialize firing rate matrix for distance bins (rows: neurons, columns: bins)
firing_rate_dist_on = zeros(size(firing_rate, 1), length(distance_bins) - 1);  % Air On
firing_rate_dist_off = zeros(size(firing_rate, 1), length(distance_bins) - 1); % Air Off

% Initialize logical arrays for air on/off within each distance bin
air_on_distance_bins = false(1, length(distance_bins) - 1);
air_off_distance_bins = false(1, length(distance_bins) - 1);

% Calculate firing rates and logical variables for Air On and Air Off periods (distance tuning)
for i = 1:length(distance_bins)-1
    % Air On period
    bin_idx_on = (ds >= distance_bins(i)) & (ds < distance_bins(i+1)) & air_on;
    firing_rate_dist_on(:, i) = mean(firing_rate(:, bin_idx_on), 2);
    air_on_distance_bins(i) = any(bin_idx_on);  % Update logical array
    
    % Air Off period
    bin_idx_off = (ds >= distance_bins(i)) & (ds < distance_bins(i+1)) & air_off;
    firing_rate_dist_off(:, i) = mean(firing_rate(:, bin_idx_off), 2);
    air_off_distance_bins(i) = any(bin_idx_off);  % Update logical array
end

% --- Plotting Example Neuron ---

neuron_idx = 1;  % Example neuron to plot

% Plot Time Tuning
figure;
subplot(1,2,1);
plot(time_bins(1:end-1), firing_rate_time_on(neuron_idx, :));
xlabel('Time (s)');
ylabel('Firing Rate (Hz)');
title('Time Tuning (Air On)');

subplot(1,2,2);
plot(time_bins(1:end-1), firing_rate_time_off(neuron_idx, :));
xlabel('Time (s)');
ylabel('Firing Rate (Hz)');
title('Time Tuning (Air Off)');

% Plot Distance Tuning
figure;
subplot(1,2,1);
plot(distance_bins(1:end-1), firing_rate_dist_on(neuron_idx, :));
xlabel('Distance (cm)');
ylabel('Firing Rate (Hz)');
title('Distance Tuning (Air On)');

subplot(1,2,2);
plot(distance_bins(1:end-1), firing_rate_dist_off(neuron_idx, :));
xlabel('Distance (cm)');
ylabel('Firing Rate (Hz)');
title('Distance Tuning (Air Off)');


%%
% Assuming you have:
% - firing_rate: Original firing rate matrix (neurons x time or distance bins)
% - firing_rate_binned_time: Firing rates binned by time (neurons x time bins)
% - firing_rate_binned_distance: Firing rates binned by distance (neurons x distance bins)
ani = 1;

firing_rate = udata{ani}.firing_rate;
firing_rate_binned_time = udataT{ani}.firing_rate;
firing_rate_binned_distance = udataD{ani}.firing_rate;
% Plot raster plot for individual neurons (e.g., neuron 1)
neuron_idx = 2;  % Index of the neuron to plot

% Option 1: Raster plot for time-binned data
figure;
imagesc(firing_rate_binned_time(neuron_idx, :));  % Plot firing rate over time
colormap(gray);
colorbar;
xlabel('Time Bins');
ylabel('Trials');
title(['Raster Plot (Time Bins) for Neuron ' num2str(neuron_idx)]);

% Option 2: Raster plot for distance-binned data
figure;
imagesc(firing_rate_binned_distance(neuron_idx, :));  % Plot firing rate over distance
colormap(gray);
colorbar;
xlabel('Distance Bins');
ylabel('Trials');
title(['Raster Plot (Distance Bins) for Neuron ' num2str(neuron_idx)]);

% Optional: Loop through multiple neurons and plot raster plots
for neuron_idx = 1:5  % Example: Loop through the first 5 neurons
    figure;
    
    % Raster plot for time-binned data
    subplot(2,1,1);
    imagesc(firing_rate_binned_time(neuron_idx, :));
    colormap(gray);
    colorbar;
    xlabel('Time Bins');
    ylabel('Trials');
    title(['Raster Plot (Time Bins) for Neuron ' num2str(neuron_idx)]);
    
    % Raster plot for distance-binned data
    subplot(2,1,2);
    imagesc(firing_rate_binned_distance(neuron_idx, :));
    colormap(gray);
    colorbar;
    xlabel('Distance Bins');
    ylabel('Trials');
    title(['Raster Plot (Distance Bins) for Neuron ' num2str(neuron_idx)]);
end
