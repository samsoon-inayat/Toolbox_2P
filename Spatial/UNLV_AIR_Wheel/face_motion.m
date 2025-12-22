%% === Load face motion-energy CSV from animal(1).pdir ===
csvFile = fullfile(animal(1).pdir, 'face_motion_energy.csv');
Tface = readtable(csvFile);

% Expect columns like: frame, time_s, motion_energy
% If time_s looks wrong/missing, rebuild from frame and fps:
if ~ismember('time_s', Tface.Properties.VariableNames) || all(Tface.time_s==0)
    fps_face = 60; % set your face camera FPS here
    t_face = double(Tface.frame) / fps_face;
else
    t_face = double(Tface.time_s);
end

if ismember('motion_energy', Tface.Properties.VariableNames)
    me = double(Tface.motion_energy);
else
    error('Expected column "motion_energy" in face_motion_energy.csv');
end

%% === Bring DAQ air-bin onto face time base ===
% b.t is DAQ time (s), b.air_bin is logical air signal at 5 kHz
t_daq   = double(b.t(:));
air_daq = logical(b.air_bin(:));

% Nearest-neighbor mapping: for each face time, find closest DAQ sample
idx = interp1(t_daq, (1:numel(t_daq))', t_face, 'nearest', 'extrap');
idx = max(1, min(numel(t_daq), idx));
air_face = air_daq(idx);   % air state for each face frame

%% === Overall air ON vs OFF means (all frames) ===
mean_on  = mean(me(air_face), 'omitnan');
mean_off = mean(me(~air_face), 'omitnan');

fprintf('Overall mean motion-energy: Air ON = %.4f, Air OFF = %.4f\n', mean_on, mean_off);

%% === Extract ON epochs (in face frame space) and do paired per-epoch means ===
% Find transitions on the face-aligned air signal
a = air_face(:);
d = diff([a(1); a]);           % prepend first sample
on_starts = find(d==1);        % indices where OFF->ON
on_ends   = find(d==-1);       % indices where ON->OFF (this index is first OFF)

% Handle if recording starts ON
if a(1)==true && (isempty(on_starts) || on_starts(1)~=1)
    on_starts = [1; on_starts];
end
% Handle if recording ends ON
if a(end)==true
    on_ends = [on_ends; numel(a)+1]; % +1 for consistent slicing
end

% Minimum epoch duration to ignore glitches (adjust if needed)
fps_face = 60;            % confirm
minDur_s = 0.2;
minSamp  = round(minDur_s * fps_face);

% Build cleaned epochs
epochs = [];
for k = 1:min(numel(on_starts), numel(on_ends))
    s = on_starts(k);
    e = on_ends(k);            % first OFF index
    if e <= s, continue; end
    if (e - s) >= minSamp
        epochs = [epochs; s e]; %#ok<AGROW>
    end
end

% Paired: each ON epoch vs immediately preceding OFF epoch
m_on  = [];
m_off = [];

prev_end = 1; % previous ON end (in face index space); off epoch starts here
for k = 1:size(epochs,1)
    s_on = epochs(k,1);
    e_on = epochs(k,2);

    % OFF epoch is from prev_end to s_on-1
    s_off = prev_end;
    e_off = s_on;

    % Only keep if OFF epoch long enough
    if (e_off - s_off) >= minSamp
        m_on(end+1,1)  = mean(me(s_on:e_on-1), 'omitnan'); %#ok<SAGROW>
        m_off(end+1,1) = mean(me(s_off:e_off-1), 'omitnan'); %#ok<SAGROW>
    end

    prev_end = e_on; % next OFF starts after this ON ends
end

fprintf('Paired epochs retained: %d\n', numel(m_on));

%% === Paired t-test across epochs (exploratory) ===
if numel(m_on) >= 3
    [h,p,ci,stats] = ttest(m_on, m_off); % paired by construction
    fprintf('Paired t-test (ON vs preceding OFF): t(%d)=%.3f, p=%.4g\n', ...
        stats.df, stats.tstat, p);
else
    warning('Not enough paired epochs for t-test.');
end

%% === Optional: quick visualization ===
figure; 
subplot(2,1,1);
plot(t_face, me); hold on;
yyaxis right; plot(t_face, double(air_face)); ylim([-0.1 1.1]);
xlabel('Time (s)'); ylabel('Air (0/1)'); title('Face motion-energy with DAQ air state');

subplot(2,1,2);
bar([mean(m_off,'omitnan'), mean(m_on,'omitnan')]);
set(gca,'XTickLabel',{'Air OFF','Air ON'});
ylabel('Mean motion-energy'); title(sprintf('Per-epoch paired means (n=%d)', numel(m_on)));
