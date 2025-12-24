function fig_dlc

vp = evalin('base','vp');

vp_l = evalin('base','vp_labeled');
vf = evalin('base','vf');
v = evalin('base','v');
mD = evalin('base','mData'); colors = mD.colors; sigColor = mD.sigColor; axes_font_size = mD.axes_font_size;
mData = mD;
animal = evalin('base','animal');
n = 0;
%%


%% ---------- Load and Clean DLC Data ----------
% Using the path structure you provided
file_path = fullfile(animal(1).pdir,'video_20251216_165824DLC_resnet50_gcamp16declimbDec18shuffle1_185000_filtered.csv');

% Set up import options to handle the triple-header (scorer, bodyparts, coords)
opts = detectImportOptions(file_path);
% Start reading from row 4 where the numeric data begins
opts.DataLines = [4, Inf]; 
opts.VariableNamingRule = 'preserve';
tbl = readtable(file_path, opts);

% Mapping likelihood columns based on your CSV image:
% Column D=4 (Front Right), G=7 (Front Left), J=10 (Hind Right), 
% M=13 (Hind Left), P=16 (Tail Base), S=19 (Nose)
lik_indices = [4, 7, 10, 13, 16, 19];
bodyparts = {'Front Right', 'Front Left', 'Hind Right', 'Hind Left', 'Tail Base', 'Nose'};

% Extract likelihood matrix [40748 frames x 6 bodyparts]
likelihoods = table2array(tbl(:, lik_indices));

%% ---------- Calculate Quality Metrics ----------
p_thresh = 0.9;
% Percentage of frames where p > 0.9 for each body part
percent_high_conf = sum(likelihoods > p_thresh, 1) / size(likelihoods, 1) * 100;

%% ---------- Updated Bodypart Colors (RGB) ----------
% These values are sampled directly from the markers in your images:
% Front Right (Blue/Indigo), Front Left (Sky Blue), Hind Right (Teal), 
% Hind Left (Bright Green), Tail Base (Orange), Nose (Red)
custom_colors = [
    0.20, 0.00, 1.00;  % Front Right (Dark Blue/Indigo)
    0.00, 0.60, 1.00;  % Front Left (Sky Blue)
    0.00, 0.80, 0.75;  % Hind Right (Teal)
    0.60, 1.00, 0.40;  % Hind Left (Bright Green)
    1.00, 0.60, 0.00;  % Tail Base (Orange)
    1.00, 0.00, 0.00   % Nose (Red)
];

%% ---------- PANEL C: DLC Quality Analysis ----------
magfac = mD.magfac;
ff = makeFigureRowsCols(110, [3 5 4 1.5], 'RowsCols', [1 2], ...
    'spaceRowsCols', [0.05 0.1], 'rightUpShifts', [0.1 0.1]);

% --- C1: Bar Graph (% Frames > 0.9) ---
subplot(1,2,1);
hold on;
for i = 1:6
    % Plot each bar individually to apply specific color
    bar(i, percent_high_conf(i), 'FaceColor', custom_colors(i,:), 'EdgeColor', 'none');
end
set(gca, 'XTick', 1:6, 'XTickLabel', bodyparts, 'XTickLabelRotation', 45);
ylabel('% Frames (p > 0.9)');
ht = title('DLC Tracking Reliability'); set(ht,'FontWeight','Normal');
ylim([0 105]); 
box off; 
format_axes(gca);

% --- C2: Likelihood CDF (Multi-line) ---
subplot(1,2,2);
hold on;
for i = 1:6
    [f, x] = ecdf(likelihoods(:, i));
    plot(x, f, 'LineWidth', 2, 'Color', custom_colors(i,:));
end

% Reference line at 0.9 threshold
line([p_thresh p_thresh], [0 1], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1);

xlabel('Likelihood');
ylabel('Cumulative Probability');
% legend(bodyparts, 'Location', 'northwest', 'FontSize', 7, 'Box', 'off');
ht = title('Likelihood CDF'); set(ht,'FontWeight','Normal');
grid on; 
box off; 
format_axes(gca);

% Save to your designated PDF folder
save_pdf(ff.hf, mD.pdf_folder, 'DLC_Quality_Analysis_Colored.pdf', 600);

%%
%% ---------- Setup & Data Loading ----------
file_path = fullfile(animal(1).pdir, 'video_20251216_165824DLC_resnet50_gcamp16declimbDec18shuffle1_185000_filtered.csv');
opts = detectImportOptions(file_path);
opts.DataLines = [4, Inf]; 
opts.VariableNamingRule = 'preserve';
tbl = readtable(file_path, opts);

% Define indices for X and Y based on your CSV structure
% x: Col 2, 5, 8, 11, 14, 17 | y: Col 3, 6, 9, 12, 15, 18
x_idx = [2, 5, 8, 11, 14, 17];
y_idx = [3, 6, 9, 12, 15, 18];
bodyparts = {'Front Right', 'Front Left', 'Hind Right', 'Hind Left', 'Tail Base', 'Nose'};


% ---------- 1. Re-establish Time Base ----------
% If time_sec is missing, we reconstruct it from the known sampling rate
% Based on your previous data, your signals have 40745 samples.
N_sig = 40748; 

% If you have your sampling frequency (e.g., 30 fps or 60 fps)
% If unknown, we can infer it if 'b.fs' exists from your earlier code
if exist('b','var') && isfield(b,'fs')
    fs = b.fs;
else
    fs = 30; % Defaulting to 30 fps; change this to your actual video FPS
end

% Create time vector in seconds
time_sec = (0:N_sig-1)' / fs;
tmin = time_sec / 60; % Convert to minutes for plotting

% Extract and truncate to match your signal time base (40745 samples)
N_target = length(time_sec);
X_coords = table2array(tbl(1:N_target, x_idx));
Y_coords = table2array(tbl(1:N_target, y_idx));

% Use the custom colors from the mouse markers
custom_colors = [
    0.20, 0.00, 1.00;  % Front Right (Blue/Indigo)
    0.00, 0.60, 1.00;  % Front Left (Sky Blue)
    0.00, 0.80, 0.75;  % Hind Right (Teal)
    0.60, 1.00, 0.40;  % Hind Left (Bright Green)
    1.00, 0.60, 0.00;  % Tail Base (Orange)
    1.00, 0.00, 0.00   % Nose (Red)
];

%% ---------- Plotting: Trajectories vs Time ----------
magfac = mD.magfac;
tmin = time_sec / 60; % Time in minutes for the x-axis

ff_traj = makeFigureRowsCols(111, [3 5 4 4], 'RowsCols', [2 1], ...
    'spaceRowsCols', [0.1 0.05], 'rightUpShifts', [0.1 0.1]);

% --- Subplot 1: X-Coordinates vs Time ---
subplot(2,1,1);
hold on;
for i = 1:6
    plot(tmin, X_coords(:,i), 'Color', custom_colors(i,:), 'LineWidth', 1);
end
ylabel('X Pixel Value');
title('Horizontal Position (X) over Time');
% legend(bodyparts, 'Location', 'eastoutside', 'FontSize', 7, 'Box', 'off');
box off; format_axes(gca);
xlim([0 3]);

% --- Subplot 2: Y-Coordinates vs Time ---
subplot(2,1,2);
hold on;
for i = 1:6
    plot(tmin, Y_coords(:,i), 'Color', custom_colors(i,:), 'LineWidth', 1);
end
ylabel('Y Pixel Value');
xlabel('Time (min)');
title('Vertical Position (Y) over Time');
box off; format_axes(gca);
xlim([0 3]);
% Save the trajectory plot
save_pdf(ff_traj.hf, mD.pdf_folder, 'DLC_Trajectories.pdf', 600);

%%
T = animal(1).b.led_sig.(key);
t_paws = T.time/60;
air_paws = double(T.is_on);