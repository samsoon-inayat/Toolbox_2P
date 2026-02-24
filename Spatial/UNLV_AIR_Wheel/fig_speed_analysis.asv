function fig_speed_analysis

% vp = evalin('base','vp');
% vf = evalin('base','vf');
% v = evalin('base','v');
mD = evalin('base','mData'); colors = mD.colors; sigColor = mD.sigColor; axes_font_size = mD.axes_font_size;
mData = mD;
animal = evalin('base','animal');
n = 0;

%%
b = animal(1).b;
b.dist = b.encoderCount * pi * 32/b.countsPerRev; % in cm
b.speed =  diff(b.dist)./diff(b.t); % in cm/sec
b.speed = double([0;b.speed]);
% b.speed = removeSpeedOutliers(b.speed);
% b.speed(b.speed < 0) = NaN;
% b.speed = fillmissing(b.speed,'linear',2,'EndValues','nearest');

samplingRate = 5000;
coeffs = ones(1, samplingRate)/samplingRate;
b.fSpeed = filter(coeffs, 1, b.speed);
%
figure(100);clf;
plot(b.t,b.fSpeed);
hold on;
plot(b.t,max(b.fSpeed)*(b.air_bin)/max(b.air_raw))
% xlim([0 200])
xlabel('Time (sec)');
ylabel('Speed cm/s');

%% with annotation


b = animal(1).b;
b.dist = b.encoderCount * pi * 32/b.countsPerRev; % in cm
b.speed =  diff(b.dist)./diff(b.t); % in cm/sec
b.speed = double([0;b.speed]);
% b.speed = removeSpeedOutliers(b.speed);
% b.speed(b.speed < 0) = NaN;
% b.speed = fillmissing(b.speed,'linear',2,'EndValues','nearest');

samplingRate = 5000;
coeffs = ones(1, samplingRate)/samplingRate;
b.fSpeed = filter(coeffs, 1, b.speed);
%
figure(100);clf;
plot(b.t,b.fSpeed);
hold on;
% plot(b.t,max(b.fSpeed)*(b.air_bin)/max(b.air_raw))
xlabel('Time (sec)');
ylabel('Speed cm/s');
ylim([0 11])
ylims = ylim;
xlims = [120 300];
xlim(xlims)
air_bin = b.air_bin;
onsets = find_rising_edge(air_bin,0.5,-1);
offsets = find_falling_edge(air_bin,-0.5,1);
time_syncM = b.t;
[TLx TLy] = ds2nfu(time_syncM(onsets(1)),ylims(2)-0);
% axes(ff.h_axes(1,1));ylims = ylim;
[BLx BLy] = ds2nfu(time_syncM(onsets(1)),ylims(1));
aH = (TLy - BLy);
len = sum(find(time_syncM(onsets)<3,1,'last'));
for ii = 1:length(onsets)
    if b.t(onsets(ii)) > xlims(2) || b.t(onsets(ii)) < xlims(1)
        continue;
    end
    [BRx BRy] = ds2nfu(time_syncM(offsets(ii)),ylims(1));
    [BLx BLy] = ds2nfu(time_syncM(onsets(ii)),ylims(1));
    aW = (BRx-BLx);
    annotation('rectangle',[BLx BLy aW aH],'facealpha',0.2,'linestyle','none','facecolor','k');
end


%% Figure raw data
magfac = mD.magfac;
ff = makeFigureRowsCols(107,[3 5 6.75 1.6],'RowsCols',[2 1],'spaceRowsCols',[0.2 -0.02],'rightUpShifts',[0.15 0.22],...
    'widthHeightAdjustment',[0 -225]);
MY = 2; ysp = 0.15285; mY = -2.5; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.3*magfac; widths = [6.4 1 2.85 1]*magfac; gap = 0.115*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 2];

axes(ff.h_axes(1,1));
plot(b.tm,b.fSpeed);
hold on;
% plot(b.tm,max(b.fSpeed)*(b.air_bin)/max(b.air_raw))
box off;
ylabel('Speed cm/s');
xlims = [0 b.tm(end)];
xlim([0 b.tm(end)]);
ylim([0 14.5])
% ylim([0 11])
ylims = ylim;

xlim(xlims)
air_bin = b.air_bin;
onsets = find_rising_edge(air_bin,0.5,-1);
offsets = find_falling_edge(air_bin,-0.5,1);
time_syncM = b.tm;
[TLx TLy] = ds2nfu(time_syncM(onsets(1)),ylims(2)-0);
% axes(ff.h_axes(1,1));ylims = ylim;
[BLx BLy] = ds2nfu(time_syncM(onsets(1)),ylims(1));
aH = (TLy - BLy);
len = sum(find(time_syncM(onsets)<3,1,'last'));
for ii = 1:length(onsets)
    if b.tm(onsets(ii)) > xlims(2) || b.tm(onsets(ii)) < xlims(1)
        continue;
    end
    [BRx BRy] = ds2nfu(time_syncM(offsets(ii)),ylims(1));
    [BLx BLy] = ds2nfu(time_syncM(onsets(ii)),ylims(1));
    aW = (BRx-BLx);
    annotation('rectangle',[BLx BLy aW aH],'facealpha',0.2,'linestyle','none','facecolor','k');
end
format_axes(gca);

axes(ff.h_axes(2,1));
plot(b.tm,b.fSpeed);
hold on;
% plot(b.tm,max(b.fSpeed)*(b.air_bin)/max(b.air_raw))
box off;
ylabel('Speed cm/s');
xlims = [2 6];
xlim([0 b.tm(end)]);
ylim([0 14.5])
% ylim([0 11])
ylims = ylim;
xlim(xlims)
air_bin = b.air_bin;
onsets = find_rising_edge(air_bin,0.5,-1);
offsets = find_falling_edge(air_bin,-0.5,1);
time_syncM = b.tm;
[TLx TLy] = ds2nfu(time_syncM(onsets(1)),ylims(2)-0);
% axes(ff.h_axes(1,1));ylims = ylim;
[BLx BLy] = ds2nfu(time_syncM(onsets(1)),ylims(1));
aH = (TLy - BLy);
len = sum(find(time_syncM(onsets)<3,1,'last'));
for ii = 1:length(onsets)
    if b.tm(onsets(ii)) > xlims(2) || b.tm(onsets(ii)) < xlims(1)
        continue;
    end
    [BRx BRy] = ds2nfu(time_syncM(offsets(ii)),ylims(1));
    [BLx BLy] = ds2nfu(time_syncM(onsets(ii)),ylims(1));
    aW = (BRx-BLx);
    annotation('rectangle',[BLx BLy aW aH],'facealpha',0.2,'linestyle','none','facecolor','k');
end
format_axes(gca);

xlabel('Time (min)');

save_pdf(ff.hf,mD.pdf_folder,'raw_speed.pdf',600);  

%%
led_paws = b.led_sig.paws;

figure(100);clf;
plot(led_paws.time, led_paws.is_on); hold on
% stairs(led_paws.time, ...
%        double(led_paws.is_on) * max(led_paws.signal_s), ...
%        'LineWidth', 1.5)
xlabel('Time (s)')
ylabel('LED signal')
% legend('Smoothed ROI signal','LEfiD ON')



%% Pre-Post Analysis to see how speed changes with the onset of air

win_pre  = 2;  % seconds
win_post = 5;  % seconds
Npre  = round(win_pre  * b.fs);
Npost = round(win_post * b.fs);


trials = [];
for i = 1:length(b.Air_r)
    idx = b.Air_r(i);
    if idx > Npre && idx + Npost <= length(b.fSpeed)
        trials(:,i) = b.fSpeed(idx-Npre : idx+Npost);
    end
end

t_evt = linspace(-win_pre, win_post, size(trials,1));


magfac = mD.magfac;
ff = makeFigureRowsCols(107,[3 5 1.75 1.5],'RowsCols',[1 1],'spaceRowsCols',[0.02 -0.02],'rightUpShifts',[0.15 0.22],...
    'widthHeightAdjustment',[10 -250]);
MY = 2; ysp = 0.15285; mY = -2.5; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.25*magfac; widths = [1.5 1 2.85 1]*magfac; gap = 0.115*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 2];
plot(t_evt, mean(trials,2), 'k', 'LineWidth', 2); hold on
plot(t_evt, mean(trials,2) + std(trials,[],2), '--k')
plot(t_evt, mean(trials,2) - std(trials,[],2), '--k')
xlabel('Time from air onset (s)')
ylabel('Speed (cm/s)')
xlim([-2 5.5]);
ylim([-1 10])
box off;
format_axes(gca);
save_pdf(ff.hf,mD.pdf_folder,'bar_graph.pdf',600);  



speed_on  = b.fSpeed(b.air_bin == 1);
speed_off = b.fSpeed(b.air_bin == 0);

[h,p] = ttest2(speed_on, speed_off);

%%
win_pre  = 2;  % seconds
win_post = 6;  % seconds
Npre  = round(win_pre  * b.fs);
Npost = round(win_post * b.fs);


trials = [];
for i = 1:length(b.Air_f)
    idx = b.Air_f(i);
    if idx > Npre && idx + Npost <= length(b.fSpeed)
        trials(:,i) = b.fSpeed(idx-Npre : idx+Npost);
    end
end

t_evt = linspace(-win_pre, win_post, size(trials,1));


magfac = mD.magfac;
ff = makeFigureRowsCols(107,[3 5 1.75 1.5],'RowsCols',[1 1],'spaceRowsCols',[0.02 -0.02],'rightUpShifts',[0.15 0.22],...
    'widthHeightAdjustment',[10 -250]);
MY = 2; ysp = 0.15285; mY = -2.5; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.25*magfac; widths = [1.5 1 2.85 1]*magfac; gap = 0.115*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 2];
plot(t_evt, mean(trials,2), 'k', 'LineWidth', 2); hold on
plot(t_evt, mean(trials,2) + std(trials,[],2), '--k')
plot(t_evt, mean(trials,2) - std(trials,[],2), '--k')
xlabel('Time from air offset (s)')
ylabel('Speed (cm/s)')
xlim([-2 win_post+0.5]);
ylim([-1 10])
box off;
format_axes(gca);
save_pdf(ff.hf,mD.pdf_folder,'bar_graph.pdf',600);  

%%
% --- Inputs ---
speed   = b.fSpeed;        % speed vector (cm/s)
air_bin = b.air_bin(:);   % logical vector (0/1)
t       = b.t(:);         % time vector (s)

fs = b.fs;                % sampling rate (Hz)

% --- Identify air ON and OFF indices ---
idx_on  = air_bin == 1;
idx_off = air_bin == 0;

% --- Mean speed ---
mean_speed_on  = mean(speed(idx_on),  'omitnan');
mean_speed_off = mean(speed(idx_off), 'omitnan');

fprintf('Mean speed (Air ON):  %.2f cm/s\n', mean_speed_on);
fprintf('Mean speed (Air OFF): %.2f cm/s\n', mean_speed_off);

% --- Bar plot (panel D) ---
figure(100);clf; hold on
bar([1 2], [mean_speed_off mean_speed_on], 0.6)
scatter(1, mean_speed_off, 60, 'k', 'filled')
scatter(2, mean_speed_on,  60, 'k', 'filled')

set(gca, 'XTick', [1 2], ...
         'XTickLabel', {'Air OFF','Air ON'}, ...
         'FontSize', 12)

ylabel('Mean speed (cm/s)')
% title('Mean running speed during Air OFF vs Air ON')
box off


[h, p] = ttest2(speed(idx_on), speed(idx_off));
fprintf('t-test p-value: %.3e\n', p);

%%
%% plot distributions Air On vs Air Off
magfac = mD.magfac;
ff = makeFigureRowsCols(107,[3 5 1.75 1.5],'RowsCols',[1 1],'spaceRowsCols',[0.02 -0.02],'rightUpShifts',[0.15 0.22],...
    'widthHeightAdjustment',[10 -250]);
MY = 2; ysp = 0.15285; mY = -2.5; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.3*magfac; widths = [1.35 1 2.85 1]*magfac; gap = 0.115*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 2];
hold on;
distD = {speed(idx_on),speed(idx_off)};

tcolors = {'b','m'};
[distDo,allVals,allValsG] = plotDistributions(distD);
minBin = min(allVals);
maxBin = max(allVals);
incr = 1;
% [ha,hb,hca] = plotDistributions(allValsG,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'do_mean','No');
[ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'do_mean','No');
set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
%     changePosition(gca,[0.129 0.15 -0.09 -0.13]);
ylim([0 100]); xlim([minBin maxBin]); %xlim([minBin 0.5]);
put_axes_labels(ha,{{'Speed (cm/s)'},[0 0 0]},{{'%'},[0 0 0]});
format_axes(ha);
[ks.h,ks.p,ks.ks2stat] = kstest2(allValsG{1},allValsG{2});
ks.DF1 = length(allValsG{1}); ks.DF2 = length(allValsG{2});
ht = set_axes_top_text_no_line(gcf,ha,'KS-Test',[0.1 -0.01 0.1 0]);set(ht,'FontSize',7);
titletxt = sprintf('%s',getNumberOfAsterisks(ks.p));
ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.33 -0.01 0.1 0]);set(ht,'FontSize',9);
legend('Air-On','Air-Off','Location','SouthEast')
% titletxt = sprintf('n = %d,',length(allValsG{1}));
% ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.015 -0.45 0 0]);set(ht,'FontSize',7,'Color','k');
% titletxt = sprintf('%d',length(allValsG{2}));
% ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.067 -0.45 0 0]);set(ht,'FontSize',7,'Color','r');
save_pdf(ff.hf,mData.pdf_folder,'firing_rate.pdf',600);

%%
%% rest vs motion FR average

air_on_idx  = b.Air_r;   % air onset indices
air_off_idx = b.Air_f;   % air offset indices

nTrials = numel(air_on_idx);

meanSpeed_ON  = nan(nTrials,1);
meanSpeed_OFF = nan(nTrials,1);

for k = 1:nTrials
    % Air ON window
    idx_on = air_on_idx(k):air_off_idx(k);
    meanSpeed_ON(k) = mean(b.fSpeed(idx_on), 'omitnan');
    medSpeed_ON(k) = median(b.fSpeed(idx_on), 'omitnan');

    % Preceding Air OFF window
    if k == 1
        idx_off = 1:(air_on_idx(k)-1);
    else
        idx_off = air_off_idx(k-1):(air_on_idx(k)-1);
    end

    meanSpeed_OFF(k) = mean(b.fSpeed(idx_off), 'omitnan');
    medSpeed_OFF(k) = median(b.fSpeed(idx_off), 'omitnan');
end

x_on  = medSpeed_ON(:);
x_off = medSpeed_OFF(:);
[p,h] = signrank(x_on, x_off);

    
tcolors = {'b','m'};
    data_C = [meanSpeed_ON meanSpeed_OFF];
    [within,dvn,xlabels] = make_within_table({'St'},[2]);
    dataT = make_between_table({data_C},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd'}});
%     ra.ranova
print_for_manuscript(ra)
   magfac = mData.magfac;
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1:3),1,2);

tcolors = {'b','m'};
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[10 4 1.25 1.5],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.2 0.2],'widthHeightAdjustment',[-550 -280]);
MY = 9; ysp = 1.5; mY = 0; ystf = 1; ysigf = 0.5;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),ra,{'St','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
% make_bars_hollow(hbs(2))
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Air-On','Air-Off'});xtickangle(30);
ylabel({'Avg. Speed'});
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0 0]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%% after review
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%

%% Parameters

nAnimals = numel(animal);

meanSpeed_ON_animal  = nan(nAnimals,1);
meanSpeed_OFF_animal = nan(nAnimals,1);

medSpeed_ON_animal  = nan(nAnimals,1);
medSpeed_OFF_animal = nan(nAnimals,1);

% Loop over animals
for a = 1:nAnimals
    
    b = animal(a).b;
    % 
    % --- recompute distance & speed (safe) ---
    % b.dist = b.encoderCount * pi * 32 / b.countsPerRev; % cm
    % b.speed = diff(b.dist)./diff(b.t); % cm/sec
    % b.speed = double([0; b.speed]);
    % removeOutliersFlag = true;   % <-- toggle here if needed
    % if removeOutliersFlag
    % 
    %     % --- robust outlier detection ---
    %     % Uses median absolute deviation (very stable)
    %     outlier_idx = isoutlier(b.fSpeed,'median');
    % 
    %     % Optional: also catch absurd physical speeds
    %     maxPhysSpeed = 200; % cm/s (adjust if needed)
    %     outlier_idx = outlier_idx | abs(b.fSpeed) > maxPhysSpeed;
    % 
    %     % --- replace with interpolation (keeps length same) ---
    %     b.fSpeed(outlier_idx) = NaN;
    %     b.fSpeed = fillmissing(b.fSpeed,'linear','EndValues','nearest');
    % 
    % end
    % % --- smooth speed (same as your code) ---
    % samplingRate = b.fs;
    % coeffs = ones(1, samplingRate)/samplingRate;
    % b.fSpeed = filter(coeffs, 1, b.speed);

    % --- air indices ---
    air_on_idx  = b.Air_r(:);
    air_off_idx = b.Air_f(:);


    % =========================================
    % FIX AIR ON/OFF ORDERING (robust)
    % =========================================
    
    % Remove any OFF that occurs before first ON
    air_off_idx(air_off_idx < air_on_idx(1)) = [];
    
    % Remove any ON that occurs after last OFF
    air_on_idx(air_on_idx > air_off_idx(end)) = [];
    
    % Make lengths equal
    nPairs = min(numel(air_on_idx), numel(air_off_idx));
    air_on_idx  = air_on_idx(1:nPairs);
    air_off_idx = air_off_idx(1:nPairs);
    
    % Safety: enforce ON < OFF for every trial
    badPairs = air_off_idx <= air_on_idx;
    
    if any(badPairs)
        warning('Fixing misordered air events...')
        air_on_idx(badPairs)  = [];
        air_off_idx(badPairs) = [];
    end
    
    nTrials = numel(air_on_idx);
    if nAnimals == 1
        start_trial = 1;
    else
        start_trial = 5;
    end
    meanSpeed_ON  = nan(length(start_trial:(nTrials)),1);
    meanSpeed_OFF = nan(length(start_trial:(nTrials)),1);

    k = 0;
   % Trial loop
    for ki = start_trial:(nTrials)
        k = k + 1;
        % ---------- AIR ON phase ----------
        idx_on = air_on_idx(k):air_off_idx(k);
        idx_on(idx_on < 1 | idx_on > numel(b.fSpeed)) = [];

        meanSpeed_ON(k) = mean(b.fSpeed(idx_on),'omitnan');
        medSpeed_ON(k) = median(b.fSpeed(idx_on),'omitnan');

        % ---------- AIR OFF phase ----------
        if k == 1
            idx_off = 1:(air_on_idx(k)-1);
        else
            idx_off = air_off_idx(k-1):(air_on_idx(k)-1);
        end

        idx_off(idx_off < 1 | idx_off > numel(b.fSpeed)) = [];

        meanSpeed_OFF(k) = mean(b.fSpeed(idx_off),'omitnan');
        medSpeed_OFF(k) = median(b.fSpeed(idx_off),'omitnan');
    end

    % Animal-level means
    meanSpeed_ON_animal(a)  = mean(meanSpeed_ON,'omitnan');
    meanSpeed_OFF_animal(a) = mean(meanSpeed_OFF,'omitnan');
    medSpeed_ON_animal(a)  = median(medSpeed_ON,'omitnan');
    medSpeed_OFF_animal(a) = median(medSpeed_OFF,'omitnan');

end
[meanSpeed_ON_animal medSpeed_ON_animal meanSpeed_OFF_animal medSpeed_OFF_animal]
[p, h, stats] = signrank(medSpeed_ON_animal, medSpeed_OFF_animal);
%%
tcolors = {'b','m'};
    data_C = [meanSpeed_ON_animal meanSpeed_OFF_animal];
    % data_C = [meanSpeed_ON meanSpeed_OFF];
    [within,dvn,xlabels] = make_within_table({'St'},[2]);
    dataT = make_between_table({data_C},dvn);
    ra = RMA(dataT,within,{0.05,{'hsd'}});
%     ra.ranova
print_for_manuscript(ra)
    % data_C = [medSpeed_ON_animal medSpeed_OFF_animal];
    % [p,h] = signrank(medSpeed_ON_animal, medSpeed_OFF_animal);
    % x_on  = medSpeed_ON_animal(:);
    % x_off = medSpeed_OFF_animal(:);
    % [p,h] = signrank(x_on, x_off);
    % meds = [median(x_on) median(x_off)];
    % 
    % q_on  = prctile(x_on,[25 75]);
    % q_off = prctile(x_off,[25 75]);
    % 
    % err_low  = [meds(1)-q_on(1), meds(2)-q_off(1)];
    % err_high = [q_on(2)-meds(1), q_off(2)-meds(2)];
    % 
    % semVar = [err_low; err_high]; % for asymmetric errorbars
    % mVar   = meds;

   magfac = mData.magfac;
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1:3),1,2);

tcolors = {'b','m'};
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[10 4 1.25 1.5],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.2 0.2],'widthHeightAdjustment',[-550 -280]);
MY = 8; ysp = 1.5; mY = 0; ystf = 1; ysigf = 0.5;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
MY = 7; ysp = 0.75; mY = 0; ystf = 1; ysigf = 0.5;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80

% [xdata,mVar,semVar,combs,p,h] = get_vals_RMA(mData,ra,{'St','hsd'},[1 2],'no');

% [hbs,maxY] = plotBarsWithSigLines(mVar,semVar,[1 2],[h p],'colors',tcolors,'sigColor','k',...
% 'ySpacing',ysp,'sigTestName','','sigLineWidth',0.25,'BaseValue',0.01,'capsize',1,...
% 'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',mData.asterisk_font_size,'barWidth',0.5,'sigLinesStartYFactor',ystf,'sigAsteriskyshift',ysigf);

 [hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),ra,{'St','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
% make_bars_hollow(hbs(2))
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Air-On','Air-Off'});xtickangle(30);
ylabel({'Avg. Speed'});
set_axes_limits(gca,[xdata(1)-0.75 xdata(end)+0.75],[mY MY]); 

% ================================
% Add paired animal dots
% ================================

x_on  = meanSpeed_ON_animal(:);
x_off = meanSpeed_OFF_animal(:);

hold on

jitter = 0.05;  % small horizontal spread
rng(1);         % reproducible

for i = 1:length(x_on)

    % jittered x positions
    x1 = xdata(1) + jitter;
    x2 = xdata(2) - jitter;

    % connecting line
    plot([x1 x2], [x_on(i) x_off(i)], '-', ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 0.75);

    % dots
    plot(x1, x_on(i), '.', ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 5);

    plot(x2, x_off(i), '.', ...
        'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 5);
end


% x_on  = medSpeed_ON_animal(:);
% x_off = medSpeed_OFF_animal(:);
% 
% hold on
% 
% jitter = -0.4;  % small horizontal spread
% rng(1);         % reproducible
% 
% for i = 1:length(x_on)
% 
%     % jittered x positions
%     x1 = xdata(1) + jitter;
%     x2 = xdata(2) + jitter;
% 
%     % connecting line
%     plot([x1 x2], [x_on(i) x_off(i)], '-', ...
%         'Color', [0.6 0.6 0.6], 'LineWidth', 0.75);
% 
%     % dots
%     plot(x1, x_on(i), 'o', ...
%         'MarkerFaceColor', 'w', ...
%         'MarkerEdgeColor', 'k', ...
%         'MarkerSize', 5);
% 
%     plot(x2, x_off(i), 'o', ...
%         'MarkerFaceColor', 'w', ...
%         'MarkerEdgeColor', 'k', ...
%         'MarkerSize', 5);
% end
% 

hold off
format_axes(gca);

% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0 0]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);

%%
%% PARAMETERS
win_pre  = 2;   % seconds
win_post = 5;   % seconds

nAnimals = numel(animal);
animal_traces = [];

for a = 1:nAnimals

    b = animal(a).b;

    removeOutliersFlag = false;   % <-- toggle here if needed
    if removeOutliersFlag

        % --- robust outlier detection ---
        % Uses median absolute deviation (very stable)
        outlier_idx = isoutlier(b.fSpeed,'median');

        % Optional: also catch absurd physical speeds
        maxPhysSpeed = 100; % cm/s (adjust if needed)
        outlier_idx = outlier_idx | abs(b.fSpeed) > maxPhysSpeed;

        % --- replace with interpolation (keeps length same) ---
        b.fSpeed(outlier_idx) = NaN;
        b.fSpeed = fillmissing(b.fSpeed,'linear','EndValues','nearest');

    end

    Npre  = round(win_pre  * b.fs);
    Npost = round(win_post * b.fs);

    trials = [];

    for i = 5:length(b.Air_r)
        idx = b.Air_r(i);

        if idx > Npre && idx + Npost <= length(b.fSpeed)
            trials(:,end+1) = b.fSpeed(idx-Npre : idx+Npost);
        end
    end

    % ✅ average within animal FIRST
    animal_mean_trace(:,a) = mean(trials,2,'omitnan');

end

t_evt = linspace(-win_pre, win_post, size(animal_mean_trace,1));

grand_mean = mean(animal_mean_trace,2,'omitnan');
grand_sem  = std(animal_mean_trace,[],2,'omitnan') ./ sqrt(nAnimals);


%%
magfac = mD.magfac;
ff = makeFigureRowsCols(107,[3 5 1.75 1.5], ...
    'RowsCols',[1 1], ...
    'spaceRowsCols',[0.02 -0.02], ...
    'rightUpShifts',[0.15 0.22], ...
    'widthHeightAdjustment',[10 -250]);

MY = 10; 
ysp = 0.15285; 
mY = -1;

stp = 0.25*magfac; 
widths = [1.5 1 2.85 1]*magfac; 
gap = 0.115*magfac;

adjust_axes(ff,[mY MY],stp,widths,gap,{''});

hold on

% % ✅ shaded SEM band (better than dashed)
% fill([t_evt fliplr(t_evt)], ...
%      [grand_mean'+grand_sem' fliplr(grand_mean'-grand_sem')], ...
%      [0.7 0.7 0.7], ...
%      'EdgeColor','none', ...
%      'FaceAlpha',0.4);
% 
% % ✅ mean line
shadedErrorBar(t_evt,grand_mean,grand_sem);
plot(t_evt, grand_mean, 'k', 'LineWidth', 2);

xlabel('Time from air onset (s)')
ylabel('Speed (cm/s)')
xlim([-2 5.5]);
ylim([0 8])

box off
format_axes(gca);

save_pdf(ff.hf,mD.pdf_folder,'air_onset_speed_ALL_ANIMALS.pdf',600);

%%

%% PARAMETERS
win_pre  = 2;   % seconds
win_post = 6;   % seconds

nAnimals = numel(animal);
animal_mean_trace = [];

for a = 1:nAnimals

    b = animal(a).b;
    
    removeOutliersFlag = false;   % <-- toggle here if needed
    if removeOutliersFlag

        % --- robust outlier detection ---
        % Uses median absolute deviation (very stable)
        outlier_idx = isoutlier(b.fSpeed,'median');

        % Optional: also catch absurd physical speeds
        maxPhysSpeed = 100; % cm/s (adjust if needed)
        outlier_idx = outlier_idx | abs(b.fSpeed) > maxPhysSpeed;

        % --- replace with interpolation (keeps length same) ---
        b.fSpeed(outlier_idx) = NaN;
        b.fSpeed = fillmissing(b.fSpeed,'linear','EndValues','nearest');

    end


    Npre  = round(win_pre  * b.fs);
    Npost = round(win_post * b.fs);

    trials = [];

    for i = 5:length(b.Air_f)

        idx = b.Air_f(i);

        if idx > Npre && idx + Npost <= length(b.fSpeed)
            trials(:,end+1) = b.fSpeed(idx-Npre : idx+Npost);
        end
    end

    % average within animal FIRST
    animal_mean_trace(:,a) = mean(trials,2,'omitnan');

end

t_evt = linspace(-win_pre, win_post, size(animal_mean_trace,1));
grand_mean = mean(animal_mean_trace,2,'omitnan');
grand_sem  = std(animal_mean_trace,[],2,'omitnan') ./ sqrt(nAnimals);
%%
magfac = mD.magfac;

ff = makeFigureRowsCols(107,[3 5 1.75 1.5], ...
    'RowsCols',[1 1], ...
    'spaceRowsCols',[0.02 -0.02], ...
    'rightUpShifts',[0.15 0.22], ...
    'widthHeightAdjustment',[10 -250]);

MY = 10;
mY = -1;

stp = 0.25*magfac;
widths = [1.5 1 2.85 1]*magfac;
gap = 0.115*magfac;

adjust_axes(ff,[mY MY],stp,widths,gap,{''});

hold on

% ⭐ population mean
shadedErrorBar(t_evt,grand_mean,grand_sem);
plot(t_evt, grand_mean, 'k', 'LineWidth', 2);

xlabel('Time from air offset (s)')
ylabel('Speed (cm/s)')
xlim([-2 win_post+0.5])
ylim([0 8])

box off
format_axes(gca);

save_pdf(ff.hf, mD.pdf_folder, 'air_offset_speed_ALL_ANIMALS.pdf', 600);

%%
%% =========================================
% MULTI-ANIMAL COLLECTION
% =========================================

nAnimals = numel(animal);

mean_speed_on_anim  = nan(nAnimals,1);
mean_speed_off_anim = nan(nAnimals,1);

all_on_speeds  = [];
all_off_speeds = [];

for a = 1:nAnimals

    b = animal(a).b;

    removeOutliersFlag = false;   % <-- toggle here if needed
    if removeOutliersFlag

        % --- robust outlier detection ---
        % Uses median absolute deviation (very stable)
        outlier_idx = isoutlier(b.fSpeed,'median');

        % Optional: also catch absurd physical speeds
        maxPhysSpeed = 100; % cm/s (adjust if needed)
        outlier_idx = outlier_idx | abs(b.fSpeed) > maxPhysSpeed;

        % --- replace with interpolation (keeps length same) ---
        b.fSpeed(outlier_idx) = NaN;
        b.fSpeed = fillmissing(b.fSpeed,'linear','EndValues','nearest');

    end

    speed   = b.fSpeed(:);
    air_bin = b.air_bin(:);

    speed(1:(b.Air_r(5)-5000)) = [];
    air_bin(1:(b.Air_r(5)-5000)) = [];

    idx_on  = air_bin == 1;
    idx_off = air_bin == 0;
    minvals_on(a) = min(speed(idx_on));
    minvals_off(a) = min(speed(idx_off));
    maxvals_on(a) = max(speed(idx_on));
    maxvals_off(a) = max(speed(idx_off));

    % ---------- per-animal means (IMPORTANT) ----------
    mean_speed_on_anim(a)  = mean(speed(idx_on),'omitnan');
    mean_speed_off_anim(a) = mean(speed(idx_off),'omitnan');

    % ---------- pooled distributions (OK for KS) ----------
    all_on_speeds  = [all_on_speeds;  speed(idx_on)];
    all_off_speeds = [all_off_speeds; speed(idx_off)];

    speed_animal_on{a} = speed(idx_on);
    speed_animal_off{a} = speed(idx_off);

end

%% plot distributions Air On vs Air Off
magfac = mD.magfac;
ff = makeFigureRowsCols(107,[3 5 1.75 1.5],'RowsCols',[1 1],'spaceRowsCols',[0.02 -0.02],'rightUpShifts',[0.15 0.22],...
    'widthHeightAdjustment',[10 -250]);
MY = 2; ysp = 0.15285; mY = -2.5; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.3*magfac; widths = [1.35 1 2.85 1]*magfac; gap = 0.115*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 2];
hold on;
distD = {speed(idx_on),speed(idx_off)};
distD = [speed_animal_on',speed_animal_off'];

tcolors = {'b','m'};
[distDo,allVals,allValsG] = plotDistributions(distD);
minBin = min(allVals);
maxBin = max(allVals);
incr = 1;
% [ha,hb,hca] = plotDistributions(allValsG,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'do_mean','No');
[ha,hb,hca] = plotDistributions(distDo,'colors',tcolors,'maxY',100,'min',minBin,'incr',incr,'max',maxBin,'do_mean','Yes');
set(gca,'FontSize',6,'FontWeight','Bold','TickDir','out','xcolor','k','ycolor','k');
%     changePosition(gca,[0.129 0.15 -0.09 -0.13]);
ylim([0 100]); xlim([minBin maxBin]); %xlim([minBin 0.5]);
put_axes_labels(ha,{{'Speed (cm/s)'},[0 0 0]},{{'%'},[0 0 0]});
format_axes(ha);
[ks.h,ks.p,ks.ks2stat] = kstest2(allValsG{1},allValsG{2});
ks.DF1 = length(allValsG{1}); ks.DF2 = length(allValsG{2});
ht = set_axes_top_text_no_line(gcf,ha,'KS-Test',[0.0351 -0.01 0.1 0]);set(ht,'FontSize',7);
titletxt = sprintf('%s',getNumberOfAsterisks(ks.p));
ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.1 -0.1 0.1 0]);set(ht,'FontSize',9);
% legend('Air-On','Air-Off','Location','SouthEast')
% titletxt = sprintf('n = %d,',length(allValsG{1}));
% ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.015 -0.45 0 0]);set(ht,'FontSize',7,'Color','k');
% titletxt = sprintf('%d',length(allValsG{2}));
% ht = set_axes_top_text_no_line(gcf,ha,titletxt,[0.067 -0.45 0 0]);set(ht,'FontSize',7,'Color','r');
save_pdf(ff.hf,mData.pdf_folder,'firing_rate.pdf',600);