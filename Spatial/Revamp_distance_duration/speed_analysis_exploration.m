%% speed analysis exploration
an = 3;
cn = 1;
pn = 1;

resp_cells_speed = resp_speed{an,4};
signals = ei{an}.plane{pn}.tP.deconv.spSigAll;
% signals = get_calcium_data(ei{an},pn);
ids_resp_cells_speed = find(resp_cells_speed);
%% plot speed and firing rate scatter plot to visually see if there is any correlation
frames = ei{an}.plane{pn}.b.frames_f;
speed = ei{an}.b.fSpeed(frames);
for ii = 1:length(ids_resp_cells_speed)
    tid = ids_resp_cells_speed(ii);
    FR = signals(tid,:);
    speed = speed(1:length(FR));
    speedgz = speed(speed>0 & FR > 0);
    FRgz = FR(speed>0 & FR > 0);

    figure(100);clf;hold on;
    % plot(FR);
    % plot(speed);
    scatterhist(speedgz,FRgz);
    pause
end

%% plot speed and firing rate together to visually see if there is any correlation
for ii = 1:length(ids_resp_cells_speed)
    tid = ids_resp_cells_speed(ii);
    FR = signals(tid,:);
    speed = speed(1:length(FR));
    speedgz = speed(speed>0);
    FRgz = FR(speed>0);

    figure(100);clf;hold on;
    plot(FR);
    plot(speed);
%     scatter(speedgz,FRgz,'.');
    pause
end

%% Grouping of FR with respect to speed bins
frames = ei{an}.plane{pn}.b.frames_f;
speed = ei{an}.b.fSpeed(frames);
speed = speed(1:size(signals,2));

min_speed = 0; max_speed = 30;
bin_incr = 1;
edges = min_speed:bin_incr:max_speed;
cens = edges - (edges(2)/2); cens(1)= [];
speedgz = speed(speed>0);
G = discretize(speedgz,edges);
mFR = [];
for ii = 1:length(ids_resp_cells_speed)
    tid = ids_resp_cells_speed(ii);
    FR = signals(tid,:);
    FRgz = FR(speed>0);
    mFR(ii,:) = splitapply(@mean,FRgz',G');
    [mFRF,mdl,rs(ii,:)] = do_gauss_fit(cens,mFR(ii,:),[],[]);
    figure(1000);clf;hold on;
    plot(cens,mFR(ii,:));
    plot(cens,mFRF);
    
end