function out = extract_led_from_roi(T, fps, varargin)
% extract_led_from_roi
%
% Extracts a robust binary LED ON/OFF signal from ROI intensity traces.
%
% INPUTS
%   T    : table with fields
%          - frame
%          - roi_mean_gray
%          - bg_mean_gray
%   fps  : frames per second (e.g. 60)
%
% OPTIONAL NAME-VALUE PAIRS
%   'MedianWindow' : window for median smoothing (default = 7)
%   'MinDuration'  : minimum ON/OFF duration in seconds (default = 0.2)
%   'OnFrac'       : fraction between low/high percentile for ON threshold (default = 0.6)
%   'OffFrac'      : fraction between low/high percentile for OFF threshold (default = 0.4)
%
% OUTPUT (struct)
%   out.time        : time vector (s)
%   out.signal_raw  : roi - background
%   out.signal_s    : smoothed signal
%   out.is_on       : logical LED ON/OFF vector
%   out.thresholds  : [thr_off thr_on]

% -------------------------------
% Parse inputs
% -------------------------------
p = inputParser;
addRequired(p, 'T');
addRequired(p, 'fps', @(x) isnumeric(x) && x>0);
addParameter(p, 'MedianWindow', 7);
addParameter(p, 'MinDuration', 0.2);
addParameter(p, 'OnFrac', 0.6);
addParameter(p, 'OffFrac', 0.4);
parse(p, T, fps, varargin{:});

w        = p.Results.MedianWindow;
minDur_s = p.Results.MinDuration;
onFrac   = p.Results.OnFrac;
offFrac  = p.Results.OffFrac;

% -------------------------------
% Time vector
% -------------------------------
time = double(T.frame) / fps;

% -------------------------------
% LED signal (ROI – background)
% -------------------------------
roi = double(T.roi_mean_gray);
bg  = double(T.bg_mean_gray);

sig_raw = roi - bg;
sig_s   = movmedian(sig_raw, w);

% -------------------------------
% Adaptive thresholds
% -------------------------------
lo = prctile(sig_s, 20);
hi = prctile(sig_s, 80);

thr_on  = lo + onFrac  * (hi - lo);
thr_off = lo + offFrac * (hi - lo);

% -------------------------------
% Hysteresis binarization
% -------------------------------
is_on = false(size(sig_s));
state = false;

for i = 1:numel(sig_s)
    if ~state && sig_s(i) >= thr_on
        state = true;
    elseif state && sig_s(i) <= thr_off
        state = false;
    end
    is_on(i) = state;
end

% -------------------------------
% Remove short glitches
% -------------------------------
minSamples = round(minDur_s * fps);

is_on = bwareaopen(is_on, minSamples);         % remove short ONs
is_on = ~bwareaopen(~is_on, minSamples);       % remove short OFFs

% -------------------------------
% Output
% -------------------------------
out.time        = time;
out.signal_raw  = sig_raw;
out.signal_s    = sig_s;
out.is_on       = is_on;
out.thresholds  = [thr_off thr_on];

end
