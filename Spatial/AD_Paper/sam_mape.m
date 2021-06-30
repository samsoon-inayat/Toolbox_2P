function md = sam_mape(x, n, trials, varargin)
% Maximum a posteriori estimation (aka Bayesian Decoding by
% neuroscientists) for Samsoon Inayat et al.
%
% Decoding error computed by Leave-one-out cross-validation (LOO-CV)
%
% Inputs:
%   x           distance/position vector
%   n           deconvolved firing rates matrix (r:time x c:neurons)
%   trials      trials vector of same length as x
% 'Name', 'Value' pairs:
%   'dt', 4         delta-t parameter for decoder (see Zhang et al. 1998)
%   'bins', 50      number of distance bins
%   'circ', false   compute circular error
%   'cv', true      leave-one-out cross-validation
%   'mle', false    use Maximum likelihood estimation instead (no prior)
%   'norm', true    normalize neuronal firing rates
%
% Outputs:
%   md.decoded      decoded distance/position vector (cross-validated)
%   md.x            binned distance
%   md.err          mean absolute error as a function of position
%   md.ops          options used in getting this result
%
% Copyright (C) 2021 - HaoRan Chang
% SPDX-License-Identifier: GPL-2.0-or-later

ops = parse_inputs(varargin);

% remove NaNs
nan_idx = isnan(x);
x(nan_idx) = [];
trials(nan_idx) = [];
n(nan_idx, :) = [];

% "straight up" the vectors
x = x(:);
trials = trials(:);

% there shall be no negative firing rate
n(n < 0) = 0;
n(isnan(n)) = 0;

% n = (n - min(n)) ./ range(n) .* ( realmax .^ (1/exp(1)) ); % normalize and fill precision
if ops.norm
    n = (n - min(n)) ./ range(n) .* exp(ops.dt * 2);
end

% get rid of silent neurons
n(:, sum(n, 1) == 0) = [];

% bin distance
bin2cm = range(x) / ops.bins;
x = discretize(x, linspace(min(x), max(x), ops.bins + 1));

decoded = zeros(size(x));
if ops.cv
    for k = unique(trials')
        loo = trials == k;
        decoded(loo) = decode(x(~loo), n(~loo, :), n(loo, :), ops);
    end
else
    decoded = decode(x, n, n, ops);
end

decoded = decoded(:);
if ops.circ
    err = min([mod(x - decoded, range(x)), mod(decoded - x, range(x))], [], 2);
else
    err = abs(x - decoded);
end
mse = mean((err .* bin2cm) .^ 2);
err = accumarray(x, err, [ops.bins, 1], @mean);

md.decoded = decoded .* bin2cm;
md.x = x .* bin2cm;
md.err = err .* bin2cm;
md.mse = mse;
md.ops = ops;


function decoded = decode(x, train, test, ops)
% Gaussian smooth firing rate estimates
% Note: for 2p data, Gaussian smoothing greatly outperforms average
% windowing for whatever reason...
train = fast_smooth(train, ops.dt);
test = fast_smooth(test, ops.dt);

% make "stack"/lambda firing rates
Px = accumarray(x, 1, [ops.bins, 1]);

lambda = arrayfun(@(ii) accumarray(x, train(:, ii), [ops.bins, 1]), 1:size(train, 2), 'UniformOutput', false);
lambda = cell2mat(lambda);
lambda = lambda ./ Px;

Px = Px ./ length(x);

% conduct decoding
n = permute(test, [3 1 2]); % dimensions: distance x time x neuron
lambda = permute(lambda, [1 3 2]);

lambda(lambda == 0) = nan;

ll = sum(n .* log(lambda) - lambda, 3, 'omitnan');
if ~ops.mle
    ll = ll + log(Px);
end
[~, decoded] = max(ll, [], 1);


function ops = parse_inputs(inputs)
ops.bins = 50;
ops.dt = 4;
ops.circ = false;
ops.cv = true;
ops.mle = false;
ops.norm = true;

count = 1;
while count < length(inputs)
    switch lower(inputs{count})
        case {'bins', 'nbins', 'num_bins'}
            ops.bins = inputs{count+1};
        case {'dt', 'deltat', 'delta-t'}
            ops.dt = inputs{count+1};
        case {'circ', 'circular', 'mod', 'modulo'}
            ops.circ = logical(inputs{count+1});
        case {'cv', 'validate', 'cross-validate'}
            ops.cv = logical(inputs{count+1});
        case 'mle'
            ops.mle = logical(inputs{count+1});
        case {'norm', 'normalize'}
            ops.norm = logical(inputs{count+1});
        otherwise
            error(['''' inputs{count} ''' is undefined.']);
    end
    count = count + 2;
end


function smoothed=fast_smooth(data,sig,dim)
% Faster 1d gaussian kernel smoothing
% Accounts for edge underestimation and NaNs
% Usage: 
%   data: matrix of size n x m where rows contain observations and columns
%      contain variables (can also be vector)
%   sig: kernel standard dev
%   dim: dimension of observations (default 1)
%
%   smoothed: smoothed data
%
% by HaoRan Chang

if sig==0
    smoothed=data;
    return
end

if nargin==2
    dim=1;
end

switch dim
    case 1
    case 2
        data=data';
    otherwise
        error('WTF dude?')
end

dataSize=size(data,1);
kernelSize=ceil(10*sig);
kernelSize=kernelSize-~mod(kernelSize,2); %avoid even kernel sizes
alpha=(kernelSize-1)/sig/2;
kernel=gausswin(kernelSize,alpha);
kernel=kernel./sum(kernel);

taper=zeros(kernelSize,1);
isnan_idx = isnan(data); % now generalized
nan_idx = [repmat(taper,1,size(data,2)); ~isnan_idx; repmat(taper,1,size(data,2))];
nan_idx = nan_idx(:);
data=[repmat(taper,1,size(data,2));data;repmat(taper,1,size(data,2))];
data=data(:);
data(isnan(data))=0;

smoothed=conv(data,kernel,'same');
smoothed=reshape(smoothed,dataSize+kernelSize*2,[]);
smoothed = smoothed./ reshape( conv(nan_idx - 1, kernel,'same') + sum(kernel), dataSize+kernelSize*2, [] ); % new performance oriented (^old) with generalization
smoothed([1:kernelSize end-kernelSize+1:end],:)=[];
smoothed(isnan_idx) = nan;

switch dim
    case 1
    case 2
        smoothed=smoothed';
end