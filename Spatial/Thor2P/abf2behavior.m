function [behavior, frameinterval] = abf2behavior(expt);
%

yy = expt.animal(1:4);
mm = expt.animal(5:6);
dd = expt.animal(7:8);
[~, exptname] = fileparts(expt.id);
filename = [expt.dirs.behavior '\', yy '_', mm, '_', dd, '_', exptname, '.abf'];

[d,si] = abfload(filename);

frame_onset = find(diff(d(:,1) >= 5) == 1);
frame_onset = frame_onset(2:end-1);  % neglect the first and last frame pulses
frameinterval = (frame_onset(end)-frame_onset(1))*si/(length(frame_onset)-1);
frameinterval = frameinterval * 1e-6;

behavior = struct;

ts = (0:size(d,1)-1) * si * 1e-6;

ch_a = find(diff(d(:,2) >=3) == 1);
ch_b = find(diff(d(:,3) >=3) == 1);
rwd = find(diff(d(:,5) >=5) == -1);  % falling edge triggers reward

if numel(rwd) == 0
    return;
end

ch_a = ch_a(ch_a >= rwd(1) & ch_a <= rwd(end));

ts = ts(ch_a+1);

pos_cum = (0:length(ch_a)-1) * 3.14/5 * 0.1;  % 500 pulses per revolution, convert to cm
pos_norm = zeros(1,length(ch_a));
pos_raw = zeros(1,length(ch_a));
trial = zeros(1,length(ch_a));
n_trial = length(rwd);

for i_trial = n_trial:-1:2
	if i_trial == n_trial
		[~,idx_st] = min(abs(ch_a-rwd(i_trial-1)));
		[~,idx_end] = min(abs(ch_a-rwd(i_trial)));
	else
		idx_end = idx_st-1;
		[~,idx_st] = min(abs(ch_a-rwd(i_trial-1)));		
	end
	pos_raw(idx_st:idx_end) = pos_cum(idx_st:idx_end)-pos_cum(idx_st);
	pos_norm(idx_st:idx_end) = pos_raw(idx_st:idx_end)./pos_raw(idx_end);
	trial(idx_st:idx_end) = i_trial-1;
end

speed = diff(pos_cum)./diff(ts);  % cm/s
ok = diff(ts) < 3.14/5/10;  % time intervals less than 0.3 sec (speed < 1 cm/s, ) are omitted
speed(~ok) = 0;
speed(circshift(~ok,-1,2)) = 0;
speed = [speed(1),speed];

behavior.ts = ts;
behavior.pos_norm = pos_norm;
behavior.trial = trial;
behavior.speed = speed;
behavior.pos_cum = pos_cum;
behavior.pos_raw = pos_raw;

end

