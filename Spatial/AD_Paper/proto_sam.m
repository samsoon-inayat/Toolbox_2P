%% load them data
clear all
load NB_decoding_data_all_cells.mat
% load NB_decoding_data_place_cells.mat

%% define parameters
% I find these parameters show the best results
circ = false;
bins = 140;
dt = 25;
cv = true;
norm = false;
mle = false;

%% pick an example
n = XsC{1, 1}; % neuronal time series
x = YsC{1, 1}; % distance

%% linearize
trials = repmat(1:size(x, 1), [size(x, 2) 1]);
trials = trials(:);
x = x';
x = x(:);
n = permute(n, [2 1 3]);
n = reshape(n, [numel(x) 1 size(n, 3)]);
n = squeeze(n);

%% run decoding
md = sam_mape(x, n, trials, 'circ', circ, 'bins', bins, 'dt', dt, 'cv', cv, 'norm', norm, 'mle', mle); % I don't know what the sampling rate is but this delta-t seems to give highest accuracy

figure
plot(md.err)
ylabel('mean decoding error (cm)')
xlabel('time (frame)')

figure
hold on
plot(md.x)
plot(md.decoded)
legend('real', 'decoded')
ylabel('distance (cm)')
xlabel('time (frame)')

%% run over all sessions
for ii = 1:size(XsC,1)
    for jj = 1:size(XsC,2)
        n = XsC{ii, jj}; % neuronal time series
        x = YsC{ii, jj}; % distance

        trials = repmat(1:size(x, 1), [size(x, 2) 1]);
        trials = trials(:);
        x = x';
        x = x(:);
        n = permute(n, [2 1 3]);
        n = reshape(n, [numel(x) 1 size(n, 3)]);
        n = squeeze(n);

        mdC(ii, jj) = sam_mape(x, n, trials, 'circ', circ, 'bins', bins, 'dt', dt, 'cv', cv, 'norm', norm, 'mle', mle);
    end
end

for ii = 1:size(XsA,1)
    for jj = 1:size(XsA,2)
        n = XsA{ii, jj}; % neuronal time series
        x = YsA{ii, jj}; % distance

        trials = repmat(1:size(x, 1), [size(x, 2) 1]);
        trials = trials(:);
        x = x';
        x = x(:);
        n = permute(n, [2 1 3]);
        n = reshape(n, [numel(x) 1 size(n, 3)]);
        n = squeeze(n);

        mdA(ii, jj) = sam_mape(x, n, trials, 'circ', circ, 'bins', bins, 'dt', dt, 'cv', cv, 'norm', norm, 'mle', mle);
    end
end


%% plot pretty figs
sem = @(x) std(x, 0, 2) / sqrt(size(x, 2));

errC = arrayfun(@(x) x.err, mdC, 'UniformOutput', false);
errC = arrayfun(@(x) cell2mat(errC(:, x)'), 1:size(errC, 2), 'UniformOutput', false);
errA = arrayfun(@(x) x.err, mdA, 'UniformOutput', false);
errA = arrayfun(@(x) cell2mat(errA(:, x)'), 1:size(errA, 2), 'UniformOutput', false);

muC = cellfun(@(x) mean(x, 2), errC, 'UniformOutput', false);
muA = cellfun(@(x) mean(x, 2), errA, 'UniformOutput', false);
devC = cellfun(@(x) sem(x), errC, 'UniformOutput', false);
devA = cellfun(@(x) sem(x), errA, 'UniformOutput', false);

figure
for ii = 1:length(muA)
    ax(ii) = subplot(1, length(muA), ii);
    hold on
    errorbar(muC{ii}, devC{ii}, 'k')
    errorbar(muA{ii}, devC{ii}, 'r')
end
legend('control', 'AD')
linkaxes(ax, 'y');
ylabel(ax(1), 'average decoding error (cm)')

