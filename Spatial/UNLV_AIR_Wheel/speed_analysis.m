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
ff = makeFigureRowsCols(107,[3 5 2.5 1.5],'RowsCols',[1 1],'spaceRowsCols',[0.02 -0.02],'rightUpShifts',[0.15 0.22],...
    'widthHeightAdjustment',[10 -250]);
MY = 2; ysp = 0.15285; mY = -2.5; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.25*magfac; widths = [2.2 1 2.85 1]*magfac; gap = 0.115*magfac;
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
ff = makeFigureRowsCols(107,[3 5 2.5 1.5],'RowsCols',[1 1],'spaceRowsCols',[0.02 -0.02],'rightUpShifts',[0.15 0.22],...
    'widthHeightAdjustment',[10 -250]);
MY = 2; ysp = 0.15285; mY = -2.5; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.25*magfac; widths = [2.2 1 2.85 1]*magfac; gap = 0.115*magfac;
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