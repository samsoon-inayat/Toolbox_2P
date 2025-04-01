function [data,data1] = get_unlv_analysis_data (ei)
%%
for ani = 1:5
   data{ani} = return_data_animal(ei{ani},ani);
   data1{ani} = return_data_animal_O(ei{ani},ani);
end




function udata = return_data_animal(sei,ani)
pl = 1;
try
    firing_rate = sei.tP.deconv.spSigAll;
catch
    firing_rate = sei.plane{1}.tP.deconv.spSigAll;
end
ca_signal = get_calcium_data(sei,pl);
frames_f = sei.plane{1}.b.frames_f(1:size(firing_rate,2))';
speed = sei.b.fSpeed(frames_f);
air = sei.b.air_puff_raw(frames_f);
belt = sei.b.photo_sensor_raw(frames_f);
light = sei.b.stim_raw(frames_f);
ts = sei.b.ts(frames_f); ts = ts - ts(1); tsm = ts/60;
ds = sei.b.dist(frames_f); ds = ds - ds(1);
tso = sei.b.ts;
dso = sei.b.dist;
frf = logical(zeros(size(tso)));
frf(frames_f) = 1;
bnb = zeros(size(ts));
bnbts = [0.08510 4.09609 19.756 23.9481;
    0.1 4.425 17.02 21.14;
    0.08 4.13 18.92 22.8;
    0.1 4.5 18.3 22.4;
    0.105 4.13 17.51 21.57];
bnb(tsm >= bnbts(ani,1) & tsm <= bnbts(ani,2)) = 1;
bnb(tsm >= bnbts(ani,3) & tsm <= bnbts(ani,4)) = 1;

%%
% Example stimulus event data (your input signal)
stimulus_event = air;  % Replace with your actual signal array

% Define the half-threshold (half of the maximum value)
threshold = max(stimulus_event) / 2;

% Convert the signal to binary based on the threshold
air_events = stimulus_event > threshold;  % Values greater than threshold become 1, others become 0
air = air_events;

% Example stimulus event data (your input signal)
stimulus_event = light;  % Replace with your actual signal array

% Define the half-threshold (half of the maximum value)
threshold = max(stimulus_event) / 2;

% Convert the signal to binary based on the threshold
light_events = stimulus_event > threshold;  % Values greater than threshold become 1, others become 0
light = light_events;

% Example stimulus event data (your input signal)
stimulus_event = belt;  % Replace with your actual signal array

% Define the half-threshold (half of the maximum value)
threshold = max(stimulus_event) / 2;

% Convert the signal to binary based on the threshold
belt_events = stimulus_event < threshold;  % Values greater than threshold become 1, others become 0
belt = belt_events;


[xoff,yoff] = get_ca_motion_correction_data(sei,1);
animal_motion = sqrt(xoff.^2 + yoff.^2);

prepost = 10000;

inmat = logical(zeros(size(sei.b.ts))); 

C = inmat; 
st = sei.b.stim_r(1)-prepost; en = sei.b.stim_r(10)+prepost;
C(st:en) = 1; C = C(frames_f);
C1 = C;

C = inmat;
st = sei.b.air_puff_r(1)-prepost; en = sei.b.air_puff_f(10)+prepost;
C(st:en) = 1; C = C(frames_f);
C2 = C;

C = inmat;
st = sei.b.air_puff_r(11)-prepost; en = sei.b.air_puff_f(20)+prepost;
C(st:en) = 1; C = C(frames_f);
C3 = C;

C = inmat;
st = sei.b.air_puff_r(21)-prepost; en = sei.b.air_puff_f(30)+prepost;
C(st:en) = 1; C = C(frames_f);
C4 = C;

C = inmat;
st = sei.b.air_puff_r(31)-prepost; en = sei.b.air_puff_f(40)+prepost;
C(st:en) = 1; C = C(frames_f);
C5 = C;

C = inmat; 
st = sei.b.stim_r(21)-prepost; en = sei.b.stim_r(30)+prepost;
C(st:en) = 1; C = C(frames_f);
C6 = C;

C = inmat;
st = sei.b.air_puff_r(41)-prepost; en = sei.b.air_puff_f(50)+prepost;
C(st:en) = 1; C = C(frames_f);
C7 = C;

udata.firing_rate = firing_rate;
udata.ca_signal = ca_signal;
udata.bnb = logical(bnb);
udata.air = air_events;
udata.light = light_events;
udata.belt = belt_events;
udata.ts = ts;
udata.tsm = tsm;
udata.ds = ds;
udata.speed = speed;
udata.animal_motion = animal_motion;
udata.C1 = C1;
udata.C2 = C2;
udata.C3 = C3;
udata.C4 = C4;
udata.C5 = C5;
udata.C6 = C6;
udata.C7 = C7;
udata.tso = tso;
udata.dso = dso;
udata.frf = frf;


function udata = return_data_animal_O(sei,ani)
pl = 1;
try
    firing_rate = sei.tP.deconv.spSigAll;
catch
    firing_rate = sei.plane{1}.tP.deconv.spSigAll;
end
ca_signal = get_calcium_data(sei,pl);
frames_f = sei.plane{1}.b.frames_f(1:size(firing_rate,2))';
frames_f = frames_f(1:size(firing_rate,2));
speed = sei.b.fSpeed;
air = sei.b.air_puff_raw;
belt = sei.b.photo_sensor_raw;
light = sei.b.stim_raw;
ts = sei.b.ts; tsm = ts/60;
ds = sei.b.dist; 
accel = [0 diff(speed)./diff(ts)];
jerk = [0 diff(accel)./diff(ts)];
snap = [0 diff(jerk)./diff(ts)];

dads = [0 diff(accel)./diff(speed)];


daccel = [0 diff(speed)./diff(ds)];
djerk = [0 diff(daccel)./diff(ds)];
dsnap = [0 diff(djerk)./diff(ds)];

frf = logical(zeros(size(ts)));
frf(frames_f) = 1;
bnb = zeros(size(ts));
bnbts = [0.08510 4.09609 19.756 23.9481;
    0.1 4.425 17.02 21.14;
    0.08 4.13 18.92 22.8;
    0.1 4.5 18.3 22.4;
    0.105 4.13 17.51 21.57];
bnb(tsm >= bnbts(ani,1) & tsm <= bnbts(ani,2)) = 1;
bnb(tsm >= bnbts(ani,3) & tsm <= bnbts(ani,4)) = 1;

%%
% Example stimulus event data (your input signal)
stimulus_event = air;  % Replace with your actual signal array

% Define the half-threshold (half of the maximum value)
threshold = max(stimulus_event) / 2;

% Convert the signal to binary based on the threshold
air_events = stimulus_event > threshold;  % Values greater than threshold become 1, others become 0
air = air_events;

% Example stimulus event data (your input signal)
stimulus_event = light;  % Replace with your actual signal array

% Define the half-threshold (half of the maximum value)
threshold = max(stimulus_event) / 2;

% Convert the signal to binary based on the threshold
light_events = stimulus_event > threshold;  % Values greater than threshold become 1, others become 0
light = light_events;

% Example stimulus event data (your input signal)
stimulus_event = belt;  % Replace with your actual signal array

% Define the half-threshold (half of the maximum value)
threshold = max(stimulus_event) / 2;

% Convert the signal to binary based on the threshold
belt_events = stimulus_event < threshold;  % Values greater than threshold become 1, others become 0
belt = belt_events;

animal_motion = zeros(size(ts));
[xoff,yoff] = get_ca_motion_correction_data(sei,1);
animal_motion(frames_f) = sqrt(xoff.^2 + yoff.^2);



prepost = 10000;

inmat = logical(zeros(size(sei.b.ts))); 

C = inmat; 
s_time = sei.b.ts(sei.b.stim_r(2))-sei.b.ts(sei.b.stim_f(1));
st = sei.b.stim_r(1)-prepost; en = sei.b.stim_f(10)+s_time*(1e6/sei.b.si)+prepost;
C(st:en) = 1;
C1 = C;

C = inmat;
s_time = sei.b.ts(sei.b.air_puff_r(2))-sei.b.ts(sei.b.air_puff_f(1));
st = sei.b.air_puff_r(1)-prepost; en = sei.b.air_puff_f(10)+s_time*(1e6/sei.b.si)+prepost;
C(st:en) = 1; 
C2 = C;

C = inmat;
s_time = sei.b.ts(sei.b.air_puff_r(12))-sei.b.ts(sei.b.air_puff_f(11));
st = sei.b.air_puff_r(11)-prepost; en = sei.b.air_puff_f(20)+s_time*(1e6/sei.b.si)+prepost;
C(st:en) = 1;
C3 = C;

C = inmat;
s_time = sei.b.ts(sei.b.air_puff_r(22))-sei.b.ts(sei.b.air_puff_f(21));
st = sei.b.air_puff_r(21)-prepost; en = sei.b.air_puff_f(30)+s_time*(1e6/sei.b.si)+prepost;
C(st:en) = 1; 
C4 = C;

C = inmat;
s_time = sei.b.ts(sei.b.air_puff_r(32))-sei.b.ts(sei.b.air_puff_f(31));
st = sei.b.air_puff_r(31)-prepost; en = sei.b.air_puff_f(40)+s_time*(1e6/sei.b.si)+prepost;
C(st:en) = 1; 
C5 = C;

C = inmat; 
s_time = sei.b.ts(sei.b.stim_r(22))-sei.b.ts(sei.b.stim_f(21));
st = sei.b.stim_r(21)-prepost; en = sei.b.stim_f(30)+s_time*(1e6/sei.b.si)+prepost;
C(st:en) = 1; 
C6 = C;

C = inmat;
s_time = sei.b.ts(sei.b.air_puff_r(42))-sei.b.ts(sei.b.air_puff_f(41));
st = sei.b.air_puff_r(41)-prepost; en = sei.b.air_puff_f(50)+s_time*(1e6/sei.b.si)+prepost;
C(st:en) = 1;
C7 = C;

air_trials = zeros(size(sei.b.ts));
air_trials = build_trials(air_trials,ts,C2 & air_events,10);
air_trials = build_trials(air_trials,ts,C3 & air_events,15);
air_trials = build_trials(air_trials,ts,C4 & air_events,15);
air_trials = build_trials(air_trials,ts,C5 & air_events,15);
air_trials = build_trials(air_trials,ts,C7 & air_events,10);
air_trials_on = air_trials;

air_trials = zeros(size(sei.b.ts));
air_trials = build_trials_off(air_trials,ts,C2 & air_events,10);
air_trials = build_trials_off(air_trials,ts,C3 & air_events,15);
air_trials = build_trials_off(air_trials,ts,C4 & air_events,15);
air_trials = build_trials_off(air_trials,ts,C5 & air_events,15);
air_trials = build_trials_off(air_trials,ts,C7 & air_events,10);
air_trials_off = air_trials;


udata.firing_rate = firing_rate;
udata.ca_signal = ca_signal;
udata.bnb = logical(bnb);
udata.air = air_events;
udata.air_trials_on = air_trials_on;
udata.air_trials_off = air_trials_off;
udata.light = light_events;
udata.belt = belt_events;
udata.ts = ts;
udata.tsm = tsm;
udata.ds = ds;
udata.speed = speed;
udata.accel = accel;
udata.animal_motion = animal_motion;
udata.C1 = C1;
udata.C2 = C2;
udata.C3 = C3;
udata.C4 = C4;
udata.C5 = C5;
udata.C6 = C6;
udata.C7 = C7;
udata.frf = frf;


function Cs = build_configurations(sei)
Cs = zeros(7,size(sei.b.ts,2));
contexts = sei.plane{1}.contexts;
for ii = 1:7
    tcontext = contexts(ii);
    startM = tcontext.markers
end

function trials = build_trials(trials,ts,series,lasttrialdur)
redges = find_rising_edge(series,0.5,0.05);
fedges = find_falling_edge(series,-0.5,0.05);

if length(redges) ~= 10
    error('something wrong');
end
redges(11) = fedges(10) + (lasttrialdur/(ts(2)-ts(1)));

for ii = 1:10
    trials(redges(ii):(fedges(ii)-1)) = ii;
end


function trials = build_trials_off(trials,ts,series,lasttrialdur)
redges = find_rising_edge(series,0.5,0.05);
fedges = find_falling_edge(series,-0.5,0.05);

if length(redges) ~= 10
    error('something wrong');
end
redges(11) = fedges(10) + (lasttrialdur/(ts(2)-ts(1)));

for ii = 1:10
    trials(fedges(ii):(redges(ii+1)-1)) = ii;
end

n = 0;