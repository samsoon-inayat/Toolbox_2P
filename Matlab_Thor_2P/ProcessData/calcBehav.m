function b = calcBehav(b,varargin)

fields_b = fields(b);

for ii = 1:length(fields_b)
    thisfield = fields_b{ii};
    cmdTxt = sprintf('b.%s = double(b.%s);',thisfield,thisfield);
    eval(cmdTxt);
%     cmdTxt = sprintf('len(ii) = numel(b.%s);',thisfield);
%     eval(cmdTxt);
end


numSamples = b.number_of_samples;
si = b.si;

ts = (0:(numSamples-1)) * si * 1e-6;
% frameinterval = (frame_onset(end)-frame_onset(1))*si/(length(frame_onset)-1);
ch_ar = b.ch_a_r;
ch_af = b.ch_a_f;
ch_br = b.ch_b_r;
ch_bf = b.ch_b_f;
photo_sensor = b.photo_sensor_f;
%%%
% ch_a = ch_ar';
% 
% % these are channel a times of the encoder ... i.e. the times at which the pulses of the encoder occured
% cha_times = ch_a*si*1e-6; 
% % From counting the encoder pulses we can find how much sDist is covered. 
% % Here we are making a sDist signal for the total sDist covered
% % There are 500 pulses per revolution i.e. 2pi radians. The radius of the
% % wheel is 5 cm ... so we find the radians per pulse and then multiply that
% % by radius to find linear sDist
% path_dist = (0:(length(ch_a)-1)) * 2*pi/500 * 5; % in cm
% 
% % Now for each pulse of channel_a we have the position starting from the first pulse and we have time for each pulse
% % then we find speed using Speed = ds/dt
% speed = diff(path_dist)./diff(cha_times); % in cm/sec
% 
% 
% speed = [0 speed];
% speed_time = cha_times;
% 
% % Find sDist as a function of time
% % this could also be used to find speed later on
% sDist = zeros(size(ts)); % initialize with all zeros
% for ii = 2:length(ch_a)
%     % find the start index where the previous pulse appears
%     st = ch_a(ii-1); 
%     % find the end index just before the occurence of current pulse
%     se = ch_a(ii)-1;
%     % the sDist covered during this time is constant as follows
%     sDist(st:se) = 5*(ii-1)*2*pi/500; % 500 pulses per revolution, 5cm radius of wheel
% end
% sDist((se+1):end) = 5*(ii)*2*pi/500;
% 
% allSpeed = zeros(size(ts));
% allSpeed(ch_a) = speed;
% 
% b.speed = allSpeed;%/max(speed);
%%%

b.ts = ts;
% b.sDist = sDist;
% There are 2000 pulses for both channel A and B combined including rising and falling edges per revolution i.e. 2pi radians. The radius of the
% wheel is 5 cm ... so we find the radians per pulse and then multiply that
% by radius to find linear sDist
if ~isfield(b,'encoderCount')
    b.encoderCount = b.dist;
end
b.dist = b.encoderCount * 2*pi/2000 * 5; % in cm
b.speed =  diff(b.dist)./diff(b.ts); % in cm/sec
b.speed = double([0 b.speed]);
b.speed = removeSpeedOutliers(b.speed);
% b.speed(b.speed < 0) = NaN;
% b.speed = fillmissing(b.speed,'linear',2,'EndValues','nearest');

samplingRate = floor(1/(b.si*1e-6));
coeffs = ones(1, samplingRate)/samplingRate;
b.fSpeed = filter(coeffs, 1, b.speed);
% figure(100);clf;plot(b.ts,b.speed);hold on;
% plot(b.ts,b.fSpeed);
% [b.eSpeed,~] = envelope(b.speed);

n = 0;


