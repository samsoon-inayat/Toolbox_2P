function mFR = get_unlv_analysis_mFR_stim_combs(tudata)
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

% whos firing_rate speed air belt light ts ds bnb
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

% Calculate mean speed for each stimulus combination
mean_speed_comb(1, 1) = mean(speed(comb_Brake_Light));  % Mean speed for Brake + Light
mean_speed_comb(1, 2) = mean(speed(comb_Brake_Air));    % Mean speed for Brake + Air
mean_speed_comb(1, 3) = mean(speed(comb_Air));          % Mean speed for Air only
mean_speed_comb(1, 4) = mean(speed(comb_Air_Light));    % Mean speed for Air + Light
mean_speed_comb(1, 5) = mean(speed(comb_No_Stimulus));  % Mean speed for No Stimulus
mean_speed_comb(1, 6) = mean(speed(locomotion));         % Locomotion
mean_speed_comb(1, 7) = mean(speed(non_locomotion));     % Non-locomotion


% Calculate mean speed for each stimulus combination
mean_speed_comb1(1, 1) = mean(animal_motion(comb_Brake_Light));  % Mean speed for Brake + Light
mean_speed_comb1(1, 2) = mean(animal_motion(comb_Brake_Air));    % Mean speed for Brake + Air
mean_speed_comb1(1, 3) = mean(animal_motion(comb_Air));          % Mean speed for Air only
mean_speed_comb1(1, 4) = mean(animal_motion(comb_Air_Light));    % Mean speed for Air + Light
mean_speed_comb1(1, 5) = mean(animal_motion(comb_No_Stimulus));  % Mean speed for No Stimulus
mean_speed_comb1(1, 6) = mean(animal_motion(locomotion));         % Locomotion
mean_speed_comb1(1, 7) = mean(animal_motion(non_locomotion));     % Non-locomotion

% % Plot mean firing rates across all combinations
% figure;
% bar(mean_firing_rates_comb);
% xlabel('Neuron');
% ylabel('Mean Firing Rate');
% legend({'Brake + Light', 'Brake + Air', 'Air', 'Air + Light'}, 'Location', 'northeastoutside');
% title('Mean Firing Rate for Different Combinations of Brake, Air, and Light');
mFR.FR = mean(mean_firing_rates_comb);
mFR.SP = mean_speed_comb;
mFR.AM = mean_speed_comb1;