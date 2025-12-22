function fig_montage

vp = evalin('base','vp');
vf = evalin('base','vf');
v = evalin('base','v');
mD = evalin('base','mData'); colors = mD.colors; sigColor = mD.sigColor; axes_font_size = mD.axes_font_size;
mmData = mD;
animal = evalin('base','animal');
n = 0;
%%
magfac = mD.magfac;
ff = makeFigureRowsCols(107,[3 5 6.75 2.5],'RowsCols',[3 8],'spaceRowsCols',[0.01 0.001],'rightUpShifts',[0.0 0.0],...
    'widthHeightAdjustment',[0 -20]);
MY = 2; ysp = 0.15285; mY = -2.5; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
% stp = 0.1*magfac; widths = [6.15 1 2.85 1]*magfac; gap = 0.115*magfac;
% adjust_axes(ff,[mY MY],stp,widths,gap,{''});
% axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 2];

%
fps = 60;
N = 4; % frames before/after (5 = ~83 ms at 60 Hz). Consider 30 for ~0.5 s.

led_sig = animal(1).b.led_sig.paws;

% --- LED signal for face camera (preferred) ---
t_led = led_sig.time(:);          % seconds
air   = logical(led_sig.is_on(:));% 0/1

% Find first rising edge
d = diff([air(1); air]);
i_on = find(d==1, 1, 'first');
if isempty(i_on), error('No rising edge found in led_sig.is_on'); end
t_on = t_led(i_on);

% Convert time to face frame index (0-based vs 1-based careful)
frame_on = round(t_on * fps) + 1;  % assuming frame 1 at t=0

totalFrames = floor(vp.Duration * vp.FrameRate);

frame_on = max(1, min(totalFrames, frame_on));

framesIdx = (frame_on-N+1):(frame_on+N);
framesIdx = framesIdx(framesIdx>=1 & framesIdx<=totalFrames);
tss = 1000*(t_led(framesIdx)-t_led(framesIdx(5)));
% Grab frames
imgs = cell(numel(framesIdx),1);
for k = 1:numel(framesIdx)
    % VideoReader is 1-based; read by frame number:
    imgs{k} = read(vp, framesIdx(k));
    axes(ff.h_axes(1,k));
    imagesc(imgs{k});axis off
    set_axes_top_text_no_line(ff.hf,ff.h_axes(1,k),sprintf('%.2f ms',tss(k)),[0.01 0.035 0 0]);
end

%

led_sig = animal(1).b.led_sig.face;

% --- LED signal for face camera (preferred) ---
t_led = led_sig.time(:);          % seconds
air   = logical(led_sig.is_on(:));% 0/1

% Find first rising edge
d = diff([air(1); air]);
i_on = find(d==1, 1, 'first');
if isempty(i_on), error('No rising edge found in led_sig.is_on'); end
t_on = t_led(i_on);

% Convert time to face frame index (0-based vs 1-based careful)
frame_on = round(t_on * fps) + 1;  % assuming frame 1 at t=0

totalFrames = floor(vp.Duration * vp.FrameRate);

frame_on = max(1, min(totalFrames, frame_on));

framesIdx = (frame_on-N+1):(frame_on+N);
framesIdx = framesIdx(framesIdx>=1 & framesIdx<=totalFrames);

% Grab frames
imgs = cell(numel(framesIdx),1);
for k = 1:numel(framesIdx)
    % VideoReader is 1-based; read by frame number:
    imgs{k} = read(vf, framesIdx(k));
    axes(ff.h_axes(2,k));
    imagesc(imgs{k});
    axis off
end

%

led_sig = animal(1).b.led_sig.pupil;

% --- LED signal for face camera (preferred) ---
t_led = led_sig.time(:);          % seconds
air   = logical(led_sig.is_on(:));% 0/1

% Find first rising edge
d = diff([air(1); air]);
i_on = find(d==1, 1, 'first');
if isempty(i_on), error('No rising edge found in led_sig.is_on'); end
t_on = t_led(i_on);

% Convert time to face frame index (0-based vs 1-based careful)
frame_on = round(t_on * fps) + 1;  % assuming frame 1 at t=0

totalFrames = floor(vp.Duration * vp.FrameRate);

frame_on = max(1, min(totalFrames, frame_on));

framesIdx = (frame_on-N+1):(frame_on+N);
framesIdx = framesIdx(framesIdx>=1 & framesIdx<=totalFrames);

% Grab frames
imgs = cell(numel(framesIdx),1);
for k = 1:numel(framesIdx)
    % VideoReader is 1-based; read by frame number:
    imgs{k} = read(v, framesIdx(k));
    axes(ff.h_axes(3,k));
    imagesc(imgs{k});axis off
end

save_pdf(ff.hf,mData.pdf_folder,'montage.pdf',600);

%%

magfac = mD.magfac;
ff = makeFigureRowsCols(107,[3 5 6.75 1.5],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.2 0.22],...
    'widthHeightAdjustment',[10 -250]);
MY = 2; ysp = 0.15285; mY = -2.5; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.35*magfac; widths = [6.3 1 2.85 1]*magfac; gap = 0.115*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 2];
hold on;
cams = {'paws','face','pupil'};
voffset = [1 1.1 1.2];
tcolors = mD.colors(1:3);
for c = 1:numel(cams)
    cam = cams{c};
    % === 2) LED – paws camera (60 Hz) ===
    T = animal(1).b.led_sig.(cam);
    t_paws = T.time/60;
    air_paws = double(T.is_on);
    plot(t_paws, air_paws * voffset(c), 'color',tcolors{c}, 'LineWidth', 0.51);
    tmax(c) = max(t_paws);
end
% === Formatting ===
xlabel('Time (min)');
ylabel('Air state');
% title('LED-based air synchronization across DAQ and video streams');
% 
% yticks([0 1.2 2.4 3.6]);
% yticklabels({'Paws','Face','Pupil'});
xlim([0 max(tmax)])
ylim([0 1.5])
legs = {'Paws','Face','Pupil',[1 0.52 1.4 1]};
% legend({'Paws','Face','Pupil'},'Location','northeast');
box off;
format_axes(gca)
putLegendH(gca,legs,tcolors)
save_pdf(ff.hf,mData.pdf_folder,'alignment.pdf',600);
%% for evaluating errors in ON duration and OFF duration of air stimulus

[dur_daq, t_on_daq, t_off_daq, dur_daq_off] = air_durations(animal(1).b.t,animal(1).b.air_bin );
[dur_paws, t_on_paws, t_off_paws, dur_paws_off] = air_durations(animal(1).b.led_sig.paws.time,animal(1).b.led_sig.paws.is_on);
[dur_face, t_on_face, t_off_face, dur_face_off] = air_durations(animal(1).b.led_sig.face.time,animal(1).b.led_sig.face.is_on);
[dur_pupil, t_on_pupil, t_off_pupil, dur_pupil_off] = air_durations(animal(1).b.led_sig.pupil.time,animal(1).b.led_sig.pupil.is_on);


S = struct();

S.daq.on   = dur_daq(:);
S.daq.off  = dur_daq_off(:);

S.paws.on  = dur_paws(:);
S.paws.off = dur_paws_off(:);

S.face.on  = dur_face(:);
S.face.off = dur_face_off(:);

S.pupil.on  = dur_pupil(:);
S.pupil.off = dur_pupil_off(:);

% common # of trials across all streams (ON durations)
N_on = min([numel(S.daq.on), numel(S.paws.on), numel(S.face.on), numel(S.pupil.on)]);
% common # of trials across all streams (OFF durations)
N_off = min([numel(S.daq.off), numel(S.paws.off), numel(S.face.off), numel(S.pupil.off)]);

fields = ["daq","paws","face","pupil"];
for f = fields
    S.(f).on  = S.(f).on(1:N_on);
    S.(f).off = S.(f).off(1:N_off);
end

E = struct();
streams = ["paws","face","pupil"];

for s = streams
    E.(s).on_ms  = 1000*(S.(s).on  - S.daq.on);
    E.(s).off_ms = 1000*(S.(s).off - S.daq.off);
end
%%
fprintf('\n=== Duration agreement vs DAQ (LED - DAQ) ===\n');

for s = streams
    % ON
    mu_on  = mean(E.(s).on_ms);
    sd_on  = std(E.(s).on_ms);
    rmse_on = sqrt(mean(E.(s).on_ms.^2));
    r_on   = corr(S.(s).on, S.daq.on, 'Rows','complete');

    % OFF
    mu_off  = mean(E.(s).off_ms);
    sd_off  = std(E.(s).off_ms);
    rmse_off = sqrt(mean(E.(s).off_ms.^2));
    r_off   = corr(S.(s).off, S.daq.off, 'Rows','complete');

    fprintf('%s:\n', upper(s));
    fprintf('  ON  error: mean = %+7.2f ms, SD = %7.2f ms, RMSE = %7.2f ms, r = %.3f (N=%d)\n', ...
        mu_on, sd_on, rmse_on, r_on, N_on);
    fprintf('  OFF error: mean = %+7.2f ms, SD = %7.2f ms, RMSE = %7.2f ms, r = %.3f (N=%d)\n', ...
        mu_off, sd_off, rmse_off, r_off, N_off);
end

%%
figure('Color','w'); 
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

% --- ON scatters ---
for i = 1:numel(streams)
    s = streams(i);
    nexttile(i);
    plot(S.daq.on, S.(s).on, 'o'); hold on;
    xlim([min(S.daq.on) max(S.daq.on)]);
    ylim([min(S.(s).on) max(S.(s).on)]);
    lo = min([xlim ylim]); hi = max([xlim ylim]);
    plot([lo hi],[lo hi],'k--'); % unity
    xlabel('DAQ ON duration (s)'); ylabel(sprintf('%s ON duration (s)', upper(s)));
    title(sprintf('%s vs DAQ (ON)', upper(s)));
    box off;
end

% --- OFF scatters ---
for i = 1:numel(streams)
    s = streams(i);
    nexttile(i+3);
    plot(S.daq.off, S.(s).off, 'o'); hold on;
    xlim([min(S.daq.off) max(S.daq.off)]);
    ylim([min(S.(s).off) max(S.(s).off)]);
    lo = min([xlim ylim]); hi = max([xlim ylim]);
    plot([lo hi],[lo hi],'k--'); % unity
    xlabel('DAQ OFF duration (s)'); ylabel(sprintf('%s OFF duration (s)', upper(s)));
    title(sprintf('%s vs DAQ (OFF)', upper(s)));
    box off;
end
%%
figure('Color','w');
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

for i = 1:numel(streams)
    s = streams(i);
    nexttile(i);
    histogram(E.(s).on_ms, 20);
    xlabel('ON duration error (ms)'); ylabel('Count');
    title(sprintf('%s (ON)', upper(s)));
    box off;
end

for i = 1:numel(streams)
    s = streams(i);
    nexttile(i+3);
    histogram(E.(s).off_ms, 20);
    xlabel('OFF duration error (ms)'); ylabel('Count');
    title(sprintf('%s (OFF)', upper(s)));
    box off;
end

%%






function [durations, on_t, off_t, durationsoff] = air_durations(t, air_bin)
% Compute air-on durations from a binary air signal

air_bin = air_bin(:) > 0;
t       = t(:);

d = diff([0; air_bin; 0]);

on_idx  = find(d == 1);
off_idx = find(d == -1) - 1;

n = min(numel(on_idx), numel(off_idx));
on_idx  = on_idx(1:n);
off_idx = off_idx(1:n);

on_t  = t(on_idx);
off_t = t(off_idx);

durations = off_t - on_t;
durationsoff = on_t(2:end) - off_t(1:(end-1));

