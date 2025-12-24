function fig_speed_analysis

vp = evalin('base','vp');
vf = evalin('base','vf');
v = evalin('base','v');
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
xlabel('Time (sec)');
ylabel('Speed cm/s');

%% Figure raw data
magfac = mD.magfac;
ff = makeFigureRowsCols(107,[3 5 6.5 1.5],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.15 0.22],...
    'widthHeightAdjustment',[10 -200]);
MY = 2; ysp = 0.15285; mY = -2.5; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.3*magfac; widths = [6.15 1 2.85 1]*magfac; gap = 0.115*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 2];
plot(b.tm,b.fSpeed);
hold on;
plot(b.tm,max(b.fSpeed)*(b.air_bin)/max(b.air_raw))
box off;
xlabel('Time (min)');
ylabel('Speed cm/s');
xlim([0 b.tm(end)]);
ylim([0 14.5])
format_axes(gca);

save_pdf(ff.hf,mD.pdf_folder,'bar_graph.pdf',600);  

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
    meanSpeed_ON(k) = mean(b.speed(idx_on), 'omitnan');

    % Preceding Air OFF window
    if k == 1
        idx_off = 1:(air_on_idx(k)-1);
    else
        idx_off = air_off_idx(k-1):(air_on_idx(k)-1);
    end

    meanSpeed_OFF(k) = mean(b.speed(idx_off), 'omitnan');
end



    
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


