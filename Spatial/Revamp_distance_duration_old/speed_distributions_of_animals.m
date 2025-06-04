% Define binning parameters
min_speed = 1; 
max_speed = 40; 
bin_incr = 1;
speed_bins = min_speed:bin_incr:max_speed;

% Load speed data from all animals
num_animals = length(udataT); % Assuming udataT contains data for all animals
all_histograms = zeros(num_animals, length(speed_bins)-1);

ff = makeFigureRowsCols(2020,[0.5 0.5 1 1],'RowsCols',[1 1],...
        'spaceRowsCols',[0.15 0.065],'rightUpShifts',[0.29 0.3],'widthHeightAdjustment',...
        [-340 -400]);    set(gcf,'color','w'); set(gcf,'Position',[10 4 1 1]);
hold on;
% Plot histograms for each animal

for ii = 1:length(ei)
    an_speeds{ii} = ei{ii}.b.fSpeed;
end

for i = 1:num_animals
    speed_data = an_speeds{i}; % Assuming 'speed' is the speed variable
    speed_hist = histcounts(speed_data, speed_bins);
    all_histograms(i, :) = speed_hist / sum(speed_hist); % Normalize
    plot(speed_bins(1:end-1), all_histograms(i, :), 'DisplayName', ['Animal ' num2str(i)],'LineWidth',0.25);
end

% Compute and plot the average histogram
avg_hist = mean(all_histograms, 1);
sem_avg = std(all_histograms)/sqrt(5);
% plot(speed_bins(1:end-1), avg_hist, 'k', 'LineWidth', 2, 'DisplayName', 'Average');
shadedErrorBar(speed_bins(1:end-1),avg_hist,sem_avg,{'color','k','linewidth',1.25},0.7);


% Add labels and legend
xlabel('Speed (cm/s)');
ylabel({'Norm. Freq.'});
% title('Speed Distribution Across Animals');
% legend show;
hold off;
format_axes(gca)
save_pdf(ff.hf,mData.pdf_folder,'speed_distributions.pdf',600);


