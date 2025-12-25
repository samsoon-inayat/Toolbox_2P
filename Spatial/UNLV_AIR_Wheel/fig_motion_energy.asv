function fig_motion_energy()


vp = evalin('base','vp');
vf = evalin('base','vf');
v = evalin('base','v');
mD = evalin('base','mData'); colors = mD.colors; sigColor = mD.sigColor; axes_font_size = mD.axes_font_size;
mData = mD;
animal = evalin('base','animal');
n = 0;
%% OF speed and motion energy with respect to air on and air off

entry = animal(1);
animal_id = entry.ID;
session_date = entry.date;
pdir = entry.pdir;

% Targets to plot
targets = {'pupil'};
num_targets = length(targets);

magfac = mD.magfac;
ff = makeFigureRowsCols(107,[3 5 4.5 1.5],'RowsCols',[1 1],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.2 0.22],'widthHeightAdjustment',[10 -250]);
MY = 2; ysp = 0.15285; mY = -2.5; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.35*magfac; widths = [3.85 1 2.85 1]*magfac; gap = 0.115*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 2];

for j = 1:num_targets
    key = targets{j};
    % subplot(num_targets, 1, j);
    
    % Get the video filename from the structure
    if isfield(entry.video.mp4, key)
        full_mp4_path = entry.video.mp4.(key);
        [~, base_name, ~] = fileparts(full_mp4_path);
        
         T = animal(1).b.led_sig.(key);
        t_paws = T.time/60;
        air_paws = double(T.is_on);
        

        % Define possible CSV paths (Reduced vs Original)
        csv_path = fullfile(pdir, [base_name, '_reduced_OF.csv']);
        if ~exist(csv_path, 'file')
            csv_path = fullfile(pdir, [base_name, '_OF.csv']);
        end
        
        if exist(csv_path, 'file')
            % 1. Read Data
            data = readtable(csv_path);
            
            % 2. Calculations
            time_sec = (data.time_ms / 1000)/60;
            speed = sqrt(data.avg_u.^2 + data.avg_v.^2);
            motion_energy = data.motion_energy;
            
            % 3. Plot Speed (Left Axis)
            yyaxis left
            plot(time_sec, speed * 0.0114 * 60, 'Color', [0, 0.447, 0.741], 'LineWidth', 0.251);
            ylabel(sprintf('%s Avg. Speed (cm/s)','OF')) % upper(key)));
            ax = gca;
            axy = ax;
            ax.YColor = [0, 0.447, 0.741];
            
            % 4. Plot Motion Energy (Right Axis)
            yyaxis right
            % plot(t_paws, air_paws*20, 'color',[0, 0.447, 0.741], 'LineWidth', 0.51);hold on
            plot(time_sec, motion_energy, 'Color', [0.85, 0.325, 0.098], 'LineStyle', '--', 'LineWidth', 0.251);
            ylabel('Energy');
            ax.YColor = [0.85, 0.325, 0.098];
            
            grid on;
            if j == num_targets
                xlabel('Time (min)');
            end
        else
            text(0.5, 0.5, ['CSV Not Found: ', key], 'HorizontalAlignment', 'center');
        end
    end
end
xlim([0 3])
onsets = find_rising_edge(air_paws,0.5,-1);
offsets = find_falling_edge(air_paws,-0.5,1);

axes(axy);ylims = ylim;
[TLx TLy] = ds2nfu(time_sec(onsets(1)),ylims(2)-0);
axes(ff.h_axes(1,1));ylims = ylim;
[BLx BLy] = ds2nfu(time_sec(onsets(1)),ylims(1));
aH = (TLy - BLy);
len = sum(find(time_sec(onsets)<3,1,'last'));
for ii = 1:len%gth(onsets)
    [BRx BRy] = ds2nfu(time_sec(offsets(ii)),ylims(1));
    [BLx BLy] = ds2nfu(time_sec(onsets(ii)),ylims(1));
    aW = (BRx-BLx);
    annotation('rectangle',[BLx BLy aW aH],'facealpha',0.2,'linestyle','none','facecolor','k');
end
box off
format_axes(gca)
axy.YColor = [0.85, 0.325, 0.098];

save_pdf(ff.hf,mData.pdf_folder,'alignment.pdf',600);

%%

time_sec = time_sec(:);
speed = speed(:);
motion_energy = motion_energy(:);
N_sig = numel(time_sec); 

% Align air_paws
if numel(air_paws) >= N_sig
    air_paws_aligned = air_paws(1:N_sig);
else
    air_paws_aligned = false(N_sig, 1);
    air_paws_aligned(1:numel(air_paws)) = air_paws;
end
air_paws_aligned = air_paws_aligned > 0;

dt = median(diff(time_sec), 'omitnan');
fps = 1/dt;
%% Pre-Post Analysis to see how speed changes with the onset of air

win_pre  = 2;  % seconds
win_post = 5;  % seconds
Npre  = round(win_pre  * 60);
Npost = round(win_post * 60);

signal = speed;

trials = [];
for i = 1:length(onsets)
    idx = onsets(i);
    if idx > Npre && idx + Npost <= length(signal)
        trials(:,i) = signal(idx-Npre : idx+Npost);
    end
end

t_evt = linspace(-win_pre, win_post, size(trials,1));


magfac = mD.magfac;
ff = makeFigureRowsCols(107,[3 5 1.75 1.5],'RowsCols',[1 1],'spaceRowsCols',[0.02 -0.02],'rightUpShifts',[0.15 0.22],...
    'widthHeightAdjustment',[10 -250]);
MY = 2; ysp = 0.15285; mY = -2.5; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.3*magfac; widths = [1.5 1 2.85 1]*magfac; gap = 0.115*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 2];
plot(t_evt, mean(trials,2), 'color',[0, 0.447, 0.741] , 'LineWidth', 2); hold on
plot(t_evt, mean(trials,2) + std(trials,[],2), 'color',[0, 0.447, 0.741],'LineStyle','--');
plot(t_evt, mean(trials,2) - std(trials,[],2), 'color',[0, 0.447, 0.741],'LineStyle','--');
xlabel('Time from air onset (s)')
ylabel('OF Avg. Speed (cm/s)')
xlim([-2 5.5]);
% ylim([-1 10])
box off;
format_axes(gca);
save_pdf(ff.hf,mD.pdf_folder,'bar_graph.pdf',600);  

%% Pre-Post Analysis to see how speed changes with the onset of air

win_pre  = 2;  % seconds
win_post = 5;  % seconds
Npre  = round(win_pre  * 60);
Npost = round(win_post * 60);

signal = motion_energy;

trials = [];
for i = 1:length(onsets)
    idx = onsets(i);
    if idx > Npre && idx + Npost <= length(signal)
        trials(:,i) = signal(idx-Npre : idx+Npost);
    end
end

t_evt = linspace(-win_pre, win_post, size(trials,1));


magfac = mD.magfac;
ff = makeFigureRowsCols(107,[3 5 1.75 1.5],'RowsCols',[1 1],'spaceRowsCols',[0.02 -0.02],'rightUpShifts',[0.15 0.22],...
    'widthHeightAdjustment',[10 -250]);
MY = 2; ysp = 0.15285; mY = -2.5; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.3*magfac; widths = [1.5 1 2.85 1]*magfac; gap = 0.115*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 2];
plot(t_evt, mean(trials,2), 'color',[0.85, 0.325, 0.098] , 'LineWidth', 2); hold on
plot(t_evt, mean(trials,2) + std(trials,[],2), 'color',[0.85, 0.325, 0.098],'LineStyle','--');
plot(t_evt, mean(trials,2) - std(trials,[],2), 'color',[0.85, 0.325, 0.098],'LineStyle','--');
xlabel('Time from air onset (s)')
ylabel('Energy')
xlim([-2 5.5]);
% ylim([-1 10])
box off;
format_axes(gca);
save_pdf(ff.hf,mD.pdf_folder,'bar_graph.pdf',600);  


%%
win_pre  = 2;  % seconds
win_post = 6;  % seconds
Npre  = round(win_pre  * 60);
Npost = round(win_post * 60);


trials = [];
for i = 1:length(offsets)
    idx = offsets(i);
    if idx > Npre && idx + Npost <= length(speed)
        trials(:,i) = speed(idx-Npre : idx+Npost);
    end
end

t_evt = linspace(-win_pre, win_post, size(trials,1));


magfac = mD.magfac;
ff = makeFigureRowsCols(107,[3 5 1.75 1.5],'RowsCols',[1 1],'spaceRowsCols',[0.02 -0.02],'rightUpShifts',[0.15 0.22],...
    'widthHeightAdjustment',[10 -250]);
MY = 2; ysp = 0.15285; mY = -2.5; titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
stp = 0.3*magfac; widths = [1.5 1 2.85 1]*magfac; gap = 0.115*magfac;
adjust_axes(ff,[mY MY],stp,widths,gap,{''});
axes_title_shifts_line = [0 0.55 0 0]; axes_title_shifts_text = [0.02 0.1 0 0]; xs_gaps = [1 2];
plot(t_evt, mean(trials,2), 'color',[0, 0.447, 0.741] , 'LineWidth', 2); hold on
plot(t_evt, mean(trials,2) + std(trials,[],2), 'color',[0, 0.447, 0.741],'LineStyle','--');
plot(t_evt, mean(trials,2) - std(trials,[],2), 'color',[0, 0.447, 0.741],'LineStyle','--');
xlabel('Time from air offset (s)')
ylabel('OF Avg. Speed (cm/s)')
xlim([-2 win_post+0.5]);
% ylim([-1 10])
box off;
format_axes(gca);
save_pdf(ff.hf,mD.pdf_folder,'bar_graph.pdf',600);  

%%
%% rest vs motion FR average

air_on_idx  = onsets;   % air onset indices
air_off_idx = offsets;   % air offset indices

nTrials = numel(air_on_idx);

meanSpeed_ON  = nan(nTrials,1);
meanSpeed_OFF = nan(nTrials,1);

for k = 1:nTrials
    % Air ON window
    idx_on = air_on_idx(k):air_off_idx(k);
    meanSpeed_ON(k) = mean(speed(idx_on), 'omitnan');

    % Preceding Air OFF window
    if k == 1
        idx_off = 1:(air_on_idx(k)-1);
    else
        idx_off = air_off_idx(k-1):(air_on_idx(k)-1);
    end

    meanSpeed_OFF(k) = mean(speed(idx_off), 'omitnan');
end



    
tcolors = {'b','c'};
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

tcolors = {'b','c'};
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[10 4 1.25 1.5],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.25 0.2],'widthHeightAdjustment',[-550 -280]);
MY = 0.35; ysp = 0.0925; mY = 0; ystf = 0.09251; ysigf = 0.015;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),ra,{'St','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
% make_bars_hollow(hbs(2))
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Air-On','Air-Off'});xtickangle(30);
ylabel({'Avg. OF Avg. Speed'});
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0 0]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);
%% rest vs motion FR average

air_on_idx  = onsets;   % air onset indices
air_off_idx = offsets;   % air offset indices

nTrials = numel(air_on_idx);

meanSpeed_ON  = nan(nTrials,1);
meanSpeed_OFF = nan(nTrials,1);

for k = 1:nTrials
    % Air ON window
    idx_on = air_on_idx(k):air_off_idx(k);
    meanSpeed_ON(k) = mean(motion_energy(idx_on), 'omitnan');

    % Preceding Air OFF window
    if k == 1
        idx_off = 1:(air_on_idx(k)-1);
    else
        idx_off = air_off_idx(k-1):(air_on_idx(k)-1);
    end

    meanSpeed_OFF(k) = mean(motion_energy(idx_off), 'omitnan');
end



    
tcolors = {'b','c'};
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

tcolors = {'r','m'};
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[10 4 1.25 1.5],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.25 0.2],'widthHeightAdjustment',[-550 -280]);
MY = 4.5; ysp = 0.925; mY = 0; ystf = 0.9251; ysigf = 0.15;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),ra,{'St','hsd',0.05},[1 2],tcolors,[mY MY ysp ystf ysigf],mData);
% make_bars_hollow(hbs(2))
format_axes(gca);
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'Air-On','Air-Off'});xtickangle(30);
ylabel({'Avg. Energy'});
% set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,2,{'Pooled'},{[0 0]});
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('C1 - AOn'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);


%% ---------- 1. Parameter Setup ----------
win_pre  = 5;  % seconds
win_post = 5;  % seconds
fs       = 60; % Sampling rate (frames per second)
Npre     = round(win_pre  * fs);
Npost    = round(win_post * fs);

signal   = pupil_area_sync; % Your data vector

% Initialize arrays for mean values
pre_means  = []; 
post_means = [];

% ---------- 2. Extraction Loop ----------
for i = 1:length(onsets)
    idx = onsets(i);
    
    % Ensure the window is within the bounds of the signal
    if idx > Npre && idx + Npost <= length(signal)
        
        % Extract the pre-event segment (5s before onset)
        pre_segment = signal(idx - Npre : idx - 1);
        
        % Extract the post-event segment (5s starting at onset)
        post_segment = signal(idx : idx + Npost);
        
        % Calculate and store the mean for this specific trial
        pre_means(end+1)  = mean(pre_segment, 'omitnan');
        post_means(end+1) = mean(post_segment, 'omitnan');
        
    end
end

% Display results for verification
fprintf('Extracted means for %d valid trials.\n', length(pre_means));

