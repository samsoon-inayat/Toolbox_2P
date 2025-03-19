function out = get_continuous_variables(data_an,air,conf,ctrl)
ts = data_an.ts;ds = data_an.ds;sp = data_an.speed;tm = data_an.animal_motion; 
ac = diff(sp)./diff(ts); ac = [0 ac];
cmdTxt = sprintf('csel = data_an.air & data_an.%s;',conf); eval(cmdTxt);
redges = find_rising_edge(csel,0.5,0.05);
fedges = find_falling_edge(csel,-0.5,0.05); 
if length(fedges) == 11 && fedges(1) == 1
    fedges(1) = [];
end
atts = []; atds = []; atsp = []; attm = []; aFR = []; atac = [];otts = [];
% figure(100);clf;plot(ts,csel);pause(0.05);
% air = 0;
if strcmp(air,'ON')
    sevent = redges; eevent = fedges; 
    for ii = 1:length(redges)
        trial = sevent(ii):eevent(ii);
        tts = ts(trial)-ts(trial(1)); tds = ds(trial)-ds(trial(1)); tsp = sp(trial); ttm = tm(trial); tac = ac(trial);
        FR = data_an.firing_rate(:,trial);
        atts = [atts tts]; atds = [atds tds]; atsp = [atsp tsp]; attm = [attm ttm]; aFR = [aFR FR]; atac = [atac tac];
        otts = [otts ts(trial)];
    end
end
if strcmp(air,'OFF')
    sevent = fedges; eevent = [redges(2:end) ceil(fedges(end)+mean(redges(2:end)-fedges(1:(length(fedges)-1))))]; % for air off
    for ii = 1:length(redges)
        trial = sevent(ii):eevent(ii);
        tts = ts(trial)-ts(trial(1)); tds = ds(trial)-ds(trial(1)); tsp = sp(trial); ttm = tm(trial); tac = ac(trial);
        FR = data_an.firing_rate(:,trial);
        atts = [atts tts]; atds = [atds tds]; atsp = [atsp tsp]; attm = [attm ttm]; aFR = [aFR FR]; atac = [atac tac];
        otts = [otts ts(trial)];
    end
end

% Creating the model for firing rate as a function of speed, distance, time, and motion
% X = [atts', atsp', atds'];%, atac', attm'];    % Predictor matrix (continuous variables)
out.time = atts';
out.dist = atds';
out.speed = atsp';
out.FR = aFR';