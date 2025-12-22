%% Load face optical flow
csvFile = fullfile(animal(1).pdir, 'face_opticalflow_roi.csv');
Tflow = readtable(csvFile);

t_face = double(Tflow.time_s);
flow   = double(Tflow.flow_mag);   % change if needed

flow_s = movmedian(flow, 7);

%%
led_air = extract_led_from_roi(b.led.face, 60);
%% LED air signal (from video LED extraction)
t_led  = double(led_air.time(:));
air_led = logical(led_air.is_on(:));

% Map LED air state onto face frames (nearest neighbor)
idx = interp1(t_led, 1:numel(t_led), t_face, 'nearest', 'extrap');
idx = max(1, min(numel(air_led), idx));
air_face = air_led(idx);

%%
figure;
yyaxis left
plot(t_face, flow_s, 'k');
ylabel('Face motion (optical flow)')

yyaxis right
plot(t_face, double(air_face), 'r');
ylim([-0.1 1.1])
ylabel('Air (ON/OFF)')

xlabel('Time (s)')
title('Face motion and air stimulus')

%%
flow_on  = flow(air_face);
flow_off = flow(~air_face);

figure;
histogram(flow_off, 30, 'Normalization','probability'); hold on
histogram(flow_on,  30, 'Normalization','probability');
legend({'Air OFF','Air ON'})
xlabel('Face motion (optical flow)')
ylabel('Probability')
title('Distribution of face motion')

%%
figure;
boxplot([flow_off; flow_on], ...
        [zeros(size(flow_off)); ones(size(flow_on))], ...
        'Labels',{'Air OFF','Air ON'});
ylabel('Face motion (optical flow)')

%%
%% Identify air ON epochs in face-time space
a = air_face(:);
d = diff([a(1); a]);

on_starts = find(d==1);
on_ends   = find(d==-1);

% Edge cases
if a(1)==1
    on_starts = [1; on_starts];
end
if a(end)==1
    on_ends = [on_ends; numel(a)+1];
end

fps_face = 60;
minDur_s = 0.2;
minSamp  = round(minDur_s * fps_face);

m_on  = [];
m_off = [];
prev_end = 1;

for k = 1:min(numel(on_starts), numel(on_ends))
    s_on = on_starts(k);
    e_on = on_ends(k)-1;

    if (e_on - s_on + 1) < minSamp
        continue
    end

    s_off = prev_end;
    e_off = s_on-1;

    if (e_off - s_off + 1) < minSamp
        prev_end = on_ends(k);
        continue
    end

    m_on(end+1,1)  = mean(flow(s_on:e_on),'omitnan'); %#ok<SAGROW>
    m_off(end+1,1) = mean(flow(s_off:e_off),'omitnan'); %#ok<SAGROW>

    prev_end = on_ends(k);
end

%%
%% Paired t-test
[h,p,~,stats] = ttest(m_on, m_off);

fprintf('Face motion Air ON vs OFF: t(%d)=%.3f, p=%.4g\n', ...
    stats.df, stats.tstat, p);

figure;
bar([mean(m_off) mean(m_on)]); hold on
errorbar([1 2], [mean(m_off) mean(m_on)], ...
         [std(m_off)/sqrt(numel(m_off)) std(m_on)/sqrt(numel(m_on))], ...
         'k','LineStyle','none')
set(gca,'XTickLabel',{'Air OFF','Air ON'})
ylabel('Mean face motion')
title(sprintf('Paired comparison (n=%d epochs)', numel(m_on)))
