function sam_mape
%% load them data
% clear all
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
fileName = fullfile(mData.pd_folder,sprintf('NB_decoding_data_place_cells.mat'));
load(fileName);
% load NB_decoding_data_place_cells.mat

%% define parameters
% I find these parameters show the best results
circ = false;
bins = 150;
dt = 60;
sig = 10;
penalty = 1e-2;
cv = true;
mle = false;

%% pick an example
% animal 5 has most accurate results
n = XsC{5, 1}; % neuronal time series
x = YsC{5, 1}; % distance

%% linearize
trials = repmat(1:size(x, 1), [size(x, 2) 1]);
trials = trials(:);
x = x';
x = x(:);
n = permute(n, [2 1 3]);
n = reshape(n, [numel(x) 1 size(n, 3)]);
n = squeeze(n);

%% run decoding
md = sam_mape(x, n, trials, 'circ', circ, 'bins', bins, 'dt', dt, 'sig', sig, 'cv', cv, 'mle', mle, 'penalty', penalty); % I don't know what the sampling rate is but this delta-t seems to give highest accuracy
md.mse

figure
plot(md.xbins(1:end-1), md.err)
ylabel('mean decoding error (cm)')
xlabel('distance (cm)')

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

        mdC(ii, jj) = sam_mape(x, n, trials, 'circ', circ, 'bins', bins, 'dt', dt, 'sig', sig, 'cv', cv, 'mle', mle, 'penalty', penalty);
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

        mdA(ii, jj) = sam_mape(x, n, trials, 'circ', circ, 'bins', bins, 'dt', dt, 'sig', sig, 'cv', cv, 'mle', mle, 'penalty', penalty);
    end
end

%% plot for individual animals
figure
for ii = 1:size(mdC, 1)
    for jj = 1:size(mdC, 2)
        ax(ii, jj) = subplot(size(mdC, 1), size(mdC, 2), (ii - 1) * size(mdC, 2) + jj);
        plot(mdC(ii, jj).x)
        hold on
        plot(mdC(ii, jj).decoded)
    end
end
linkaxes(ax, 'xy')

figure
for ii = 1:size(mdA, 1)
    for jj = 1:size(mdA, 2)
        ax(ii, jj) = subplot(size(mdA, 1), size(mdA, 2), (ii - 1) * size(mdA, 2) + jj);
        plot(mdA(ii, jj).x)
        hold on
        plot(mdA(ii, jj).decoded)
    end
end
linkaxes(ax, 'xy')


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
    errorbar(md.xbins(1:end-1), muC{ii}, devC{ii}, 'k')
    errorbar(md.xbins(1:end-1), muA{ii}, devC{ii}, 'r')
end
legend('control', 'AD')
linkaxes(ax, 'y');
ylabel(ax(1), 'average decoding error (cm)')
xlabel(ax(1), 'distance (cm)')



%% load them data
clear all
load NB_Decoding_Dist_Time/NB_decoding_data_all_cells.mat
% load NB_decoding_data_place_cells.mat

% missing data point
Dist{1, 1}(:, end) = [];
Space{1, 1}(:, end) = [];
T{1, 1}(:, end) = [];

%% define parameters
% I find these parameters show the best results
circ = false;
bins = 150;
dt = 60;
sig = 10;
penalty = 1e-2;
cv = true;
mle = false;

%% run over all sessions
for ii = 1:size(Rasters,1)
    for jj = 1:size(Rasters,2)
        n = Rasters{ii, jj}; % neuronal time series
        d = Dist{ii, jj}; % distance
        x = Space{ii, jj}; % distance
        t = T{ii, jj}; % distance

        trials = repmat(1:size(x, 1), [size(x, 2) 1]);
        trials = trials(:);
        x = x';
        x = x(:);
        d = d';
        d = d(:);
        t = t';
        t = t(:);
        n = permute(n, [2 1 3]);
        n = reshape(n, [numel(x) 1 size(n, 3)]);
        n = squeeze(n);

        mdX(ii, jj) = sam_mape(x, n, trials, 'circ', circ, 'bins', bins, 'dt', dt, 'sig', sig, 'cv', cv, 'mle', mle, 'penalty', penalty);
        mdD(ii, jj) = sam_mape(d, n, trials, 'circ', circ, 'bins', bins, 'dt', dt, 'sig', sig, 'cv', cv, 'mle', mle, 'penalty', penalty);
        mdT(ii, jj) = sam_mape(t, n, trials, 'circ', circ, 'bins', round(range(t) / .1), 'dt', dt, 'sig', sig, 'cv', cv, 'mle', mle, 'penalty', penalty);
    end
end

%% individual time decoding
figure
for ii = 1:size(mdT, 1)
    for jj = 1:size(mdT, 2)
        ax(ii, jj) = subplot(size(mdT, 1), size(mdT, 2), (ii - 1) * size(mdT, 2) + jj);
        plot(mdT(ii, jj).x)
        hold on
        plot(mdT(ii, jj).decoded)
    end
end
linkaxes(ax, 'xy')


%% individual time decoding
figure
for ii = 1:size(mdD, 1)
    for jj = 1:size(mdD, 2)
        ax(ii, jj) = subplot(size(mdD, 1), size(mdD, 2), (ii - 1) * size(mdD, 2) + jj);
        plot(mdD(ii, jj).x)
        hold on
        plot(mdD(ii, jj).decoded)
    end
end
linkaxes(ax, 'xy')


%% plot pretty figs
sem = @(x) std(x, 0, 2) / sqrt(size(x, 2));

errX = arrayfun(@(x) x.err, mdX, 'UniformOutput', false);
errX = arrayfun(@(x) cell2mat(errX(:, x)'), 1:size(errX, 2), 'UniformOutput', false);
errD = arrayfun(@(x) x.err, mdD, 'UniformOutput', false);
errD = arrayfun(@(x) cell2mat(errD(:, x)'), 1:size(errD, 2), 'UniformOutput', false);

l = min(arrayfun(@(x) length(x.err), mdT), [], 'all');
errT = arrayfun(@(x) x.err(1:l), mdT, 'UniformOutput', false);
errT = arrayfun(@(x) cell2mat(errT(:, x)'), 1:size(errT, 2), 'UniformOutput', false);

muX = cellfun(@(x) mean(x, 2), errX, 'UniformOutput', false);
devX = cellfun(@(x) sem(x), errX, 'UniformOutput', false);
muD = cellfun(@(x) mean(x, 2), errD, 'UniformOutput', false);
devD = cellfun(@(x) sem(x), errD, 'UniformOutput', false);
muT = cellfun(@(x) mean(x, 2), errT, 'UniformOutput', false);
devT = cellfun(@(x) sem(x), errT, 'UniformOutput', false);

figure
for ii = 1:length(muX)
    ax(1, ii) = subplot(3, length(muX), ii);
    errorbar(mdX(1).xbins(1:end-1), muX{ii}, devX{ii}, 'k')
    ax(2, ii) = subplot(3, length(muD), 3 + ii);
    errorbar(mdX(1).xbins(1:end-1), muD{ii}, devD{ii}, 'k')
    ax(3, ii) = subplot(3, length(muT), 6 + ii);
    errorbar(linspace(0, l*.1, l), muT{ii}, devT{ii}, 'k')
end
linkaxes(ax(1:2, :), 'xy');
linkaxes(ax(3, :), 'xy');
ylabel(ax(1, 1), 'average decoding error (cm)')
ylabel(ax(2, 1), 'average decoding error (cm)')
ylabel(ax(3, 1), 'average decoding error (sec)')
xlabel(ax(2, 1), 'distance (cm)')
xlabel(ax(3, 1), 'time (sec)')

title(ax(1, 2), 'space')
title(ax(2, 2), 'distance')
title(ax(3, 2), 'time')
