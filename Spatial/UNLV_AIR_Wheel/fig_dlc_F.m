function fig_dlc_F

vp = evalin('base','vp');

vp_l = evalin('base','vp_labeled');
vf = evalin('base','vf');
v = evalin('base','v');
mD = evalin('base','mData'); colors = mD.colors; sigColor = mD.sigColor; axes_font_size = mD.axes_font_size;
mData = mD;
animal = evalin('base','animal');
n = 0;
%%

%% ---------- Load and Clean Face DLC Data ----------
% Using the path structure for the face tracking file
file_path = fullfile(animal(1).pdir, 'face_1440x1080_60_20251216_165824DLC_resnet50_face16Dec18shuffle1_125000_filtered.csv');

% Set up import options to handle the triple-header
opts = detectImportOptions(file_path);
opts.DataLines = [4, Inf]; % Numeric data begins at row 4
opts.VariableNamingRule = 'preserve';
tbl = readtable(file_path, opts);

% Defined bodyparts from the face project
bodyparts = {'Nose', 'nostril1', 'nostril2', 'Mouth', 'bodypart3', 'objectA', ...
             'wisker1', 'wisker2', 'wisker3', 'wisker4', 'wisker5', ...
             'wisker6', 'wisker7', 'wisker8', 'wisker9', 'wisker10'};

% Mapping likelihood columns (1-based indexing):
% Each bodypart has 3 columns (x, y, likelihood) starting after the frame index (Col 1).
% Likelihoods are in columns: 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49
lik_indices = 4:3:49; 

% Extract likelihood matrix [Frames x 16 bodyparts]
likelihoods = table2array(tbl(:, lik_indices));

% Quick Quality Check: Percentage of frames with confidence > 0.9
p_thresh = 0.9;
percent_high_conf = mean(likelihoods > p_thresh) * 100;

fprintf('DLC Face Tracking Loaded: %d frames for %d bodyparts.\n', size(likelihoods, 1), length(bodyparts));
%% ---------- Calculate Quality Metrics ----------
p_thresh = 0.9;
% Percentage of frames where p > 0.9 for each body part
percent_high_conf = sum(likelihoods > p_thresh, 1) / size(likelihoods, 1) * 100;

%% ---------- 1. Load Face DLC Data (16 Parts) ----------
file_path = fullfile(animal(1).pdir, 'face_1440x1080_60_20251216_165824DLC_resnet50_face16Dec18shuffle1_125000_filtered.csv');
opts = detectImportOptions(file_path);
opts.DataLines = [4, Inf]; 
opts.VariableNamingRule = 'preserve';
tbl = readtable(file_path, opts);

% Define bodyparts in the exact order found in your CSV
bodyparts = {'Nose', 'Nostril1', 'Nostril2', 'Mouth', 'M-Side1', 'M-Side1', ...
             };

% Likelihood columns: 4, 7, 10, 13, 16, 19, 22... up to 49
lik_indices = 4:3:49; 
likelihoods = table2array(tbl(:, lik_indices));
percent_high_conf = mean(likelihoods > 0.9, 1) * 100;

%% ---------- Updated Colors for Visible Face Parts (RGB) ----------
% These values match the 6 markers in image_b6c7bd.jpg exactly:
% Nose (Purple), nostril1 (Cyan), nostril2 (Magenta), Mouth (Blue), 
% bodypart3 (Green), objectA (Orange)
custom_colors_vis = [
    0.65, 0.40, 1.00;  % 1 Nose (Purple)
    0.00, 0.90, 0.95;  % 2 nostril1 (Cyan)
    1.00, 0.00, 1.00;  % 3 nostril2 (Magenta/Vibrant Pink)
    0.15, 0.40, 1.00;  % 4 Mouth (Blue)
    0.40, 0.95, 0.20;  % 5 bodypart3 (Green)
    1.00, 0.55, 0.05   % 6 objectA (Orange)
];

%% ---------- 3. Quality Analysis Plot ----------
magfac = mD.magfac;
ff = makeFigureRowsCols(110, [3 5 5 1.5], 'RowsCols', [1 2], ...
    'spaceRowsCols', [0.05 0.1], 'rightUpShifts', [0.1 0.1]);

% --- C1: Bar Graph ---
subplot(1,2,1); hold on;
for i = 1:length(bodyparts)
    bar(i, percent_high_conf(i), 'FaceColor', custom_colors_vis(i,:), 'EdgeColor', 'none');
end
set(gca, 'XTick', 1:6, 'XTickLabel', bodyparts, 'XTickLabelRotation', 45, 'FontSize', 6);
ylabel('% Frames (p > 0.9)'); title('Face Tracking Reliability');
ylim([0 105]); box off; format_axes(gca);

% --- C2: Likelihood CDF ---
subplot(1,2,2); hold on;
for i = 1:length(bodyparts)
    [f, x] = ecdf(likelihoods(:, i));
    plot(x, f, 'LineWidth', 0.25, 'Color', custom_colors(i,:));
end
xline(0.9, '--', 'Color', [0.5 0.5 0.5]);
xlabel('Likelihood'); ylabel('Cumulative Prob'); title('Likelihood CDF');
grid on; box off; format_axes(gca);

save_pdf(ff.hf, mD.pdf_folder, 'DLC_Face_Quality_Calibrated.pdf', 600);
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
    fs = 60; % Defaulting to 30 fps; change this to your actual video FPS
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
%%
T = animal(1).b.led_sig.("paws");
t_paws = T.time/60;
air_paws = double(T.is_on);
%% ---------- Plotting: Trajectories vs Time ----------
magfac = mD.magfac;
tmin = time_sec / 60; % Time in minutes for the x-axis

ff = makeFigureRowsCols(111, [3 5 6.9 1.5], 'RowsCols', [1 2], ...
    'spaceRowsCols', [0.1 5], 'rightUpShifts', [0.1 0.2],'widthHeightAdjustment',[0 -300]);
MY = 1920; ysp = 0.15285; mY = -2.5; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.35*magfac; widths = [3.12 3.12 2.85 1]*magfac; gap = 0.25*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
% axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 2];


% --- Subplot 1: X-Coordinates vs Time ---
axes(ff.h_axes(1,1))
hold on;
for i = 1:6
    miny(i) = min(X_coords(:,i)); maxy(i) = max(X_coords(:,i));
    plot(tmin, X_coords(:,i), 'Color', custom_colors(i,:), 'LineWidth', 0.25);
end
ylim([200 1400])
ylims = ylim;
plot(t_paws,air_paws*ylims(2),'k','LineWidth',0.25)
ylabel('X Pixel Value');xlabel('Time (min)');
% title('Horizontal Position (X) over Time');
% legend(bodyparts, 'Location', 'eastoutside', 'FontSize', 7, 'Box', 'off');
% box off; format_axes(gca);
xlim([0 3]);
onsets = find_rising_edge(air_paws,0.5,-1);
offsets = find_falling_edge(air_paws,-0.5,1);
ylims = ylim;
[TLx TLy] = ds2nfu(time_sec(onsets(1)),ylims(2)-0);
% axes(ff.h_axes(1,1));ylims = ylim;
[BLx BLy] = ds2nfu(time_sec(onsets(1)),ylims(1));
aH = (TLy - BLy);
len = sum(find(time_sec(onsets)<3,1,'last'));
for ii = 1:len%gth(onsets)
    [BRx BRy] = ds2nfu(time_sec(offsets(ii)),ylims(1));
    [BLx BLy] = ds2nfu(time_sec(onsets(ii)),ylims(1));
    aW = (BRx-BLx);
    annotation('rectangle',[BLx BLy aW aH],'facealpha',0.2,'linestyle','none','facecolor','k');
end
box off
format_axes(gca)

% --- Subplot 2: Y-Coordinates vs Time ---
axes(ff.h_axes(1,2))
hold on;
for i = 1:6
    plot(tmin, Y_coords(:,i), 'Color', custom_colors(i,:), 'LineWidth', 0.25);
end
ylim([10 1000])
ylims = ylim;
plot(t_paws,air_paws*ylims(2),'k','LineWidth',0.25)

ylabel('Y Pixel Value');
xlabel('Time (min)');
% title('Vertical Position (Y) over Time');
% box off; format_axes(gca);
xlim([0 3]);

onsets = find_rising_edge(air_paws,0.5,-1);
offsets = find_falling_edge(air_paws,-0.5,1);

ylims = ylim;
[TLx TLy] = ds2nfu(time_sec(onsets(1)),ylims(2)-0);
% axes(ff.h_axes(1,1));ylims = ylim;
[BLx BLy] = ds2nfu(time_sec(onsets(1)),ylims(1));
aH = (TLy - BLy);
len = sum(find(time_sec(onsets)<3,1,'last'));
for ii = 1:len%gth(onsets)
    [BRx BRy] = ds2nfu(time_sec(offsets(ii)),ylims(1));
    [BLx BLy] = ds2nfu(time_sec(onsets(ii)),ylims(1));
    aW = (BRx-BLx);
    annotation('rectangle',[BLx BLy aW aH],'facealpha',0.2,'linestyle','none','facecolor','k');
end
box off
format_axes(gca)

% Save the trajectory plot
save_pdf(ff.hf, mD.pdf_folder, 'DLC_Trajectories.pdf', 600);


%% ---------- 1. Calculate Speeds for All Body Parts ----------
% Ensure fs is defined (sampling rate)
if ~exist('fs','var'), fs = 60; end 
dt = 1/fs;

% Initialize speed matrix [Frames x 6 bodyparts]
DLC_speeds = nan(size(X_coords));

for i = 1:6
    % Calculate velocity (difference between frames)
    dx = diff(X_coords(:,i)) * fs; % Change in pixels per second
    dy = diff(Y_coords(:,i)) * fs;
    
    % Collapse X and Y into Speed (Magnitude of Velocity Vector)
    % We pad with a NaN at the start to maintain vector length
    DLC_speeds(2:end, i) = sqrt(dx.^2 + dy.^2);
end

% physical_dist_cm = 5.0;  
% pixel_dist_px = 200; % Measure this from a still frame using 'imdistline'
px_to_cm = 0.0114; 

% Update the speeds by multiplying pixels by the calibration factor
DLC_speeds_cm = DLC_speeds * px_to_cm;

%% ---------- 2. Plot: Speeds vs Time (Panel F Style) ----------
magfac = mD.magfac;
tmin = time_sec / 60;
magfac = mD.magfac;
ff = makeFigureRowsCols(107,[3 5 6.5 1.5],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.2 0.22],'widthHeightAdjustment',[10 -250]);
MY = 100; ysp = 0.15285; mY = -2.5; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.35*magfac; widths = [6.4 1 2.85 1]*magfac; gap = 0.115*magfac; adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 2];

% --- Subplot 1: Individual Speed Traces ---
axes(ff.h_axes(1,1))
hold on;
for i = 1:6
    plot(tmin, DLC_speeds_cm(:,i), 'Color', custom_colors(i,:), 'LineWidth', 0.8);
end
ylabel('Speed (pixels/s)');
title('DLC Speed: All Body Parts');
% legend(bodyparts, 'Location', 'eastoutside', 'FontSize', 7, 'Box', 'off');
box off; format_axes(gca);
ylim([0 200]);xlim([0 3]);
ylims = ylim;
plot(t_paws,air_paws*ylims(2),'k')


onsets = find_rising_edge(air_paws,0.5,-1);
offsets = find_falling_edge(air_paws,-0.5,1);
xlims = xlim;

[TLx TLy] = ds2nfu(time_sec(onsets(1)),ylims(2)-0);
% axes(ff.h_axes(1,1));ylims = ylim;
[BLx BLy] = ds2nfu(time_sec(onsets(1)),ylims(1));
aH = (TLy - BLy);
len = sum(find(time_sec(onsets)<xlims(2),1,'last'));
for ii = 1:len%gth(onsets)
    [BRx BRy] = ds2nfu(time_sec(offsets(ii)),ylims(1));
    [BLx BLy] = ds2nfu(time_sec(onsets(ii)),ylims(1));
    aW = (BRx-BLx);
    annotation('rectangle',[BLx BLy aW aH],'facealpha',0.2,'linestyle','none','facecolor','k');
end
box off
format_axes(gca)

save_pdf(ff.hf, mD.pdf_folder, 'DLC_Speeds_Combined.pdf', 600);

%% ---------- 3. Statistical Comparison (Air-ON vs Air-OFF) ----------
% Align with air epochs
is_on = air_paws(1:length(avg_speed)) > 0;
speed_on  = avg_speed(is_on & isfinite(avg_speed));
speed_off = avg_speed(~is_on & isfinite(avg_speed));

if ~isempty(speed_on) && ~isempty(speed_off)
    [~, p_speed] = ttest2(speed_on, speed_off);
    fprintf('Global DLC Speed: Air-ON vs OFF, p = %.2e\n', p_speed);
end


%% ===================== DLC SPEED + AIR ON/OFF =====================

% Inputs assumed from your script:
%   X_coords (N x 6), Y_coords (N x 6)
%   time_sec (N x 1) in seconds (DLC time base)
%   T = animal(1).b.led_sig.paws;  T.time, T.is_on  (LED time base)
%   bodyparts (1x6 cell), custom_colors (6x3)

% ---- 0) Basic checks
N = size(X_coords,1);
assert(numel(time_sec)==N, 'time_sec length must match X_coords rows.');

% ---- 1) Resample/align air signal to DLC time base
% LED time is in seconds already (but you divided by 60 earlier for t_paws)
t_led  = double(T.time(:));                % seconds
air_led = double(strcmpi(string(T.is_on), "True") | double(T.is_on)); % robust

% If T.is_on is already logical numeric, above still works.
% Make sure time vectors are monotonic
[~, idxSort] = sort(t_led);
t_led = t_led(idxSort);
air_led = air_led(idxSort);

% Resample LED air to DLC time points
air_dlc = interp1(t_led, air_led, time_sec, 'previous', 0);
air_dlc(isnan(air_dlc)) = 0;
air_dlc = air_dlc > 0.5;

% ---- 2) Compute per-bodypart speed in pixels/s
fs_dlc = 1/median(diff(time_sec));  % inferred FPS from DLC time base
dx = [zeros(1,6); diff(X_coords,1,1)];
dy = [zeros(1,6); diff(Y_coords,1,1)];
speed_pix_s = sqrt(dx.^2 + dy.^2) * fs_dlc * 0.0114 ;   % N x 6

% Optional smoothing (helps readability)
speed_pix_s_sm = movmedian(speed_pix_s, 5, 1);

% ---- 3) Find air epochs on DLC time base
air = air_dlc(:);
on_idx  = find(diff([0; air])== 1);   % rising edges
off_idx = find(diff([air; 0])==-1);   % falling edges
nTr = min(numel(on_idx), numel(off_idx));
on_idx = on_idx(1:nTr);
off_idx = off_idx(1:nTr);

% sanity: ensure each onset precedes offset
good = off_idx > on_idx;
on_idx = on_idx(good);
off_idx = off_idx(good);
nTr = numel(on_idx);

% ---- 4) Trial-wise paired comparison: Air-ON vs preceding Air-OFF window
% Define a "pre" window immediately before each onset (e.g., 1.0 s)
preWin_s = 1.0;
preSamp = max(1, round(preWin_s * fs_dlc));

mean_on  = nan(nTr,6);
mean_off = nan(nTr,6);

for k = 1:nTr
    idx_on  = on_idx(k):off_idx(k);
    idx_off = max(1, on_idx(k)-preSamp):on_idx(k)-1;

    if isempty(idx_off) || numel(idx_on) < 3
        continue
    end

    mean_on(k,:)  = mean(speed_pix_s_sm(idx_on,:), 1, 'omitnan');
    mean_off(k,:) = mean(speed_pix_s_sm(idx_off,:), 1, 'omitnan');
end

% Drop trials with NaNs
validTr = all(isfinite(mean_on),2) & all(isfinite(mean_off),2);
mean_on  = mean_on(validTr,:);
mean_off = mean_off(validTr,:);
nValid = size(mean_on,1);

% Paired t-test per bodypart
p_t = nan(1,6);
tstat = nan(1,6);
for i = 1:6
    [~, p_t(i), ~, stats] = ttest(mean_on(:,i), mean_off(:,i)); % paired
    tstat(i) = stats.tstat;
end

% ---- 5) Plot 1: Speed time-series with shaded air epochs
figure(501); clf
tmin = time_sec/60;

for i = 1:6
    subplot(3,2,i); hold on

    % Shaded air windows
    yl = [0, max(speed_pix_s_sm(:,i), [], 'omitnan')*1.05 + eps];
    ylim(yl)

    for k = 1:nTr
        x1 = time_sec(on_idx(k))/60;
        x2 = time_sec(off_idx(k))/60;
        patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], ...
              'k', 'FaceAlpha', 0.12, 'EdgeColor', 'none');
    end

    plot(tmin, speed_pix_s_sm(:,i), 'Color', custom_colors(i,:), 'LineWidth', 0.5);
    xlim([0 3]) % match your 3-min window
    xlabel('Time (min)')
    ylabel('Speed (px/s)')
    title(sprintf('%s', bodyparts{i}))
    box off
end

sgtitle('DLC point speed over time (shaded = Air ON)')

% ---- 6) Plot 2: Onset-aligned speed (mean ± SD)
win_pre  = 2;  % seconds before onset
win_post = 5;  % seconds after onset
Lpre  = round(win_pre  * fs_dlc);
Lpost = round(win_post * fs_dlc);
t_rel = (-Lpre:Lpost)'/fs_dlc;

% Build onset-aligned matrices: (time x trials x bodyparts)
M = numel(t_rel);
speed_onset = nan(M, nTr, 6);

for k = 1:nTr
    c = on_idx(k);
    idx = (c-Lpre):(c+Lpost);
    if idx(1) < 1 || idx(end) > N
        continue
    end
    speed_onset(:,k,:) = speed_pix_s_sm(idx,:);
end

figure(502); clf
for i = 1:6
    subplot(3,2,i); hold on
    X = squeeze(speed_onset(:,:,i)); % M x nTr
    mu = mean(X, 2, 'omitnan');
    sd = std(X, 0, 2, 'omitnan');

    plot(t_rel, mu, 'Color', custom_colors(i,:), 'LineWidth', 1.5);
    plot(t_rel, mu+sd, '--', 'Color', custom_colors(i,:), 'LineWidth', 0.75);
    plot(t_rel, mu-sd, '--', 'Color', custom_colors(i,:), 'LineWidth', 0.75);
    xline(0, 'k-');
    xlabel('Time from air onset (s)')
    ylabel('Speed (px/s)')
    title(bodyparts{i})
    box off
end
sgtitle('Onset-aligned DLC speed (mean ± SD)')

% ---- 7) Plot 3: Trial-wise Air-ON vs Air-OFF (paired) bar + dots
figure(503); clf
% tiledlayout(1,6,'Padding','compact','TileSpacing','tight');
%%
magfac = mD.magfac;
ff = makeFigureRowsCols(107,[3 5 6.75 1.5],'RowsCols',[1 6],'spaceRowsCols',[0.01 0.04],'rightUpShifts',[0.051 0.2],...
    'widthHeightAdjustment',[-45 -500]);


for i = 1:6
    % subplot(1,6,i);
    axes(ff.h_axes(1,i));% nexttile
    hold on
    mOn  = mean(mean_on(:,i),  'omitnan');
    mOff = mean(mean_off(:,i), 'omitnan');
    seOn  = std(mean_on(:,i),  'omitnan')/sqrt(nValid);
    seOff = std(mean_off(:,i), 'omitnan')/sqrt(nValid);

    hb = bar([1 2], [mOff mOn]); % OFF then ON
    errorbar([1 2], [mOff mOn], [seOff seOn], 'k.', 'LineWidth', 1);
    set(hb,'FaceColor',custom_colors(i,:))
    % paired dots
    for k = 1:nValid
        plot([1 2], [mean_off(k,i) mean_on(k,i)], '-', 'Color', [0 0 0 0.15]);
    end

    set(gca,'XTick',[1 2],'XTickLabel',{'Air-OFF','Air-ON'});xtickangle(30)
    if i == 1
        ylabel('Mean speed (cm/s)')
    end
    % title(sprintf('%s | p=%.3g', bodyparts{i}, p_t(i)))
    ht = title(sprintf('%s | p<0.001', bodyparts{i})); set(ht,'FontWeight','Normal')
    box off
    format_axes(gca)
end
ht = sgtitle(sprintf('Paired per-trial Air-ON vs pre-Air-OFF (N=%d trials)', nValid));set(ht,'FontSize',8,'FontWeight','Normal')

save_pdf(gcf, mD.pdf_folder, 'DLC_bar_air_on_vs_off.pdf', 600);
