function out = get_unlv_analysis_multiFunc(tudata,anname)
% ani = 1;
% tudata = udata{ani};
firing_rate = tudata.firing_rate;
bnb = tudata.bnb;
air = tudata.air;
light = tudata.light;
belt = tudata.belt;
ts = tudata.ts;
tsm = tudata.tsm;
ds = tudata.ds;
speed = tudata.speed;
animal_motion = tudata.animal_motion;
C1 = tudata.C1; C2 = tudata.C2; C3 = tudata.C3; 
C4 = tudata.C4; C5 = tudata.C5; C6 = tudata.C6; 
C7 = tudata.C7;

if strcmp(anname,'C345_AON_AOFF_avgspeed')
    % Define binary labels for No-Brake conditions
    no_brake_C3 = C3;  % No-Brake C3
    no_brake_C4 = C4;  % No-Brake C4
    no_brake_C5 = C5;  % No-Brake C5
    
    % Air On and Air Off periods
    air_on = (air == 1);    % Air On periods
    air_off = (air == 0);   % Air Off periods
    
    % Initialize result arrays for mean speeds during Air On and Air Off for No-Brake conditions
    mean_speed_air_on = zeros(1, 3);  % For C3, C4, C5 (No-Brake, Air On)
    mean_speed_air_off = zeros(1, 3); % For C3, C4, C5 (No-Brake, Air Off)
    
    % Calculate mean speed for Air On periods within No-Brake configurations
    mean_speed_air_on(1) = mean(speed(no_brake_C3 & air_on));  % Air On for C3
    mean_speed_air_on(2) = mean(speed(no_brake_C4 & air_on));  % Air On for C4
    mean_speed_air_on(3) = mean(speed(no_brake_C5 & air_on));  % Air On for C5
    
    % Calculate mean speed for Air Off periods within No-Brake configurations
    mean_speed_air_off(1) = mean(speed(no_brake_C3 & air_off));  % Air Off for C3
    mean_speed_air_off(2) = mean(speed(no_brake_C4 & air_off));  % Air Off for C4
    mean_speed_air_off(3) = mean(speed(no_brake_C5 & air_off));  % Air Off for C5
    out.msaon = mean_speed_air_on;
    out.msaoff = mean_speed_air_off;
    return;
end

if strcmp(anname,'MI_BNB')
    % Assuming you have:
    % - firing_rate: Firing rates binned by time (neurons x time bins)
    % - time_bins: Time bin indices (same length as time dimension of firing_rate)
    % - bnb: Binary vector indicating brake condition (1 = brake on, 0 = no brake)
    
    % Find the indices for brake on and brake off
    brake_on_idx = find(bnb == 1);   % Indices where brake is on
    brake_off_idx = find(bnb == 0);  % Indices where brake is off
    
    % Extract the firing rates for brake on and brake off
    firing_rate_brake_on = firing_rate(:, brake_on_idx);   % Firing rates during brake on
    firing_rate_brake_off = firing_rate(:, brake_off_idx); % Firing rates during brake off
    
    % Extract the corresponding time bins for brake on and brake off
    time_bins_brake_on = ts(brake_on_idx);   % Time bins during brake on
    time_bins_brake_off = ts(brake_off_idx); % Time bins during brake off
    
    % Reshape the data to match the format expected by MutualInformation function
    firing_rate_brake_on_reshaped = reshape(firing_rate_brake_on, [], 1);   % Reshape firing rates for brake on
    firing_rate_brake_off_reshaped = reshape(firing_rate_brake_off, [], 1); % Reshape firing rates for brake off
    
    % Reshape the time bins to match the size of the firing rates
    time_bins_brake_on_reshaped = repmat(time_bins_brake_on, size(firing_rate_brake_on, 1), 1);   % Replicate time bins
    time_bins_brake_off_reshaped = repmat(time_bins_brake_off, size(firing_rate_brake_off, 1), 1); % Replicate time bins
    
    % Compute mutual information between firing rate and time for brake on condition
    mutual_info_brake_on = MutualInformation(firing_rate_brake_on_reshaped, time_bins_brake_on_reshaped(:));
    
    % Compute mutual information between firing rate and time for brake off condition
    mutual_info_brake_off = MutualInformation(firing_rate_brake_off_reshaped, time_bins_brake_off_reshaped(:));
    
    % Display the results
    disp(['Mutual Information (Brake On, Firing Rate vs Time): ', num2str(mutual_info_brake_on)]);
    disp(['Mutual Information (Brake Off, Firing Rate vs Time): ', num2str(mutual_info_brake_off)]);
    
    % Plot the comparison
    figure;
    bar([mutual_info_brake_on, mutual_info_brake_off]);
    set(gca, 'XTickLabel', {'Brake On', 'Brake Off'});
    ylabel('Mutual Information (bits)');
    title('Mutual Information: Firing Rate vs Time (Brake On vs Brake Off)');

 
end