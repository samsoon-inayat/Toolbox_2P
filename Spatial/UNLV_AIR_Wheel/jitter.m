% --- DAQ edges ---
ev_daq = get_edges(b.t, b.air_bin);

% --- VIDEO edges ---
ev_vid = get_edges(led_paws.time, led_paws.is_on);

% Align by first rising edge
t0_daq = ev_daq.rise_t(1);
t0_vid = ev_vid.rise_t(1);

shift = t0_daq - t0_vid;           % add this to video time to align to DAQ
t_vid_aligned = led_paws.time + shift;

ev_vidA = get_edges(t_vid_aligned, led_paws.is_on);


%%
maxAbsDiff = 0.200;  % seconds tolerance for matching (adjust if needed)

M_on  = match_edges(ev_daq.rise_t, ev_vidA.rise_t, maxAbsDiff);
M_off = match_edges(ev_daq.fall_t, ev_vidA.fall_t, maxAbsDiff);

fprintf('ON edges matched:  %d/%d\n', M_on.n_matched, numel(ev_daq.rise_t));
fprintf('OFF edges matched: %d/%d\n', M_off.n_matched, numel(ev_daq.fall_t));

% Jitter stats (video - DAQ)
jit_on  = M_on.dt(~isnan(M_on.dt));
jit_off = M_off.dt(~isnan(M_off.dt));

fprintf('ON jitter:  mean = %.3f ms, SD = %.3f ms\n',  1e3*mean(jit_on),  1e3*std(jit_on));
fprintf('OFF jitter: mean = %.3f ms, SD = %.3f ms\n',  1e3*mean(jit_off), 1e3*std(jit_off));


%%
K = min(numel(ev_daq.rise_t), numel(ev_vid.rise_t));
K = min(K, 30);  % use first 30 rises (or all)

% Use median difference as robust shift
shift_med = median(ev_daq.rise_t(1:K) - ev_vid.rise_t(1:K));

t_vid_aligned = led_paws.time + shift_med;
ev_vidA = get_edges(t_vid_aligned, led_paws.is_on);

% Then compute jitter as above

%%
% Use ON jitter for drift check
t_ref = ev_daq.rise_t(~isnan(M_on.dt));
dt    = M_on.dt(~isnan(M_on.dt));

p = polyfit(t_ref - t_ref(1), dt, 1);   % dt = p1*time + p2
drift_rate = p(1);                      % sec drift per sec

fprintf('Estimated drift rate: %.3f ms/min\n', 1e3*drift_rate*60);
