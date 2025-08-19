function S = pairwise_stats_from_X(all_cells, cols, varargin)
% PAIRWISE_STATS_FROM_X  Pairwise set stats + meta-analysis across animals.
%
% S = pairwise_stats_from_X(all_cells, cols, 'Tail','right','Weights',[])
%
% INPUTS
%   all_cells : A x M cell array. Row = animal. Column = population.
%               Each entry is a logical/0-1 vector (cells in that animal).
%   cols      : vector of column indices to analyze (e.g., 18 MI or 18 PC).
%
% Name-Value:
%   'Tail'    : 'right' (enrichment, default), 'left' (depletion), or 'both'.
%   'Weights' : A-length vector for Stouffer combining (default = ones).
%
% OUTPUT struct S fields (C = numel(cols), A = #animals):
%   .per_animal.J, .DICE, .CI_total               (C x C x A)
%   .per_animal.k, .a_only, .b_only, .n           (C x C x A)  % 2x2 counts
%   .per_animal.p_fisher                          (C x C x A)  % one-tailed or two-sided per 'Tail'
%   .per_animal.LOR                               (C x C x A)  % log-odds (0.5-corrected)
%
%   .mean.J, .sem.J     ;  .mean.DICE, .sem.DICE  ;  .mean.CI_total, .sem.CI_total
%
%   .stouffer.p    (C x C)  % combined p across animals (per 'Tail')
%   .stouffer.Z    (C x C)  % combined Z
%
%   .meta.OR       (C x C)  % pooled odds ratio (fixed-effect)
%   .meta.CI       (C x C x 2) % [low, high]
%   .meta.p_two_sided (C x C)  % meta p (two-sided)
%   .meta.Q, .meta.I2          % heterogeneity
%
% NOTES
%   - CI_total = 100 * |A∩B| / N_cells (symmetric "conjunction %").
%   - Jaccard = |A∩B| / |A∪B| ; DICE = 2|A∩B| / (2|A∩B| + |A\B| + |B\A|).
%   - For meta log-odds: fixed-effect with Haldane–Anscombe 0.5 correction.
%   - Run FDR (BH) outside on S.stouffer.p or S.meta.p_two_sided if needed.

% ---------- parse args ----------
p = inputParser;
p.addParameter('Tail','right',@(s)ischar(s) || isstring(s));
p.addParameter('Weights',[],@(v)isnumeric(v) || isempty(v));
p.parse(varargin{:});
tail = lower(string(p.Results.Tail));
if ~ismember(tail,["right","left","both"])
    error('Tail must be ''right'', ''left'', or ''both''.');
end

A = size(all_cells,1);              % animals
C = numel(cols);                    % populations (e.g., 18)
if isempty(p.Results.Weights)
    wS = ones(A,1);                 % Stouffer weights
else
    wS = p.Results.Weights(:);
    if numel(wS)~=A, error('Weights must match #animals'); end
end

% ---------- prealloc per-animal ----------
J        = nan(C,C,A);
DICE     = nan(C,C,A);
CI_total = nan(C,C,A);
k_all    = nan(C,C,A);
a_all    = nan(C,C,A);
b_all    = nan(C,C,A);
n_all    = nan(C,C,A);
pF       = nan(C,C,A);
LOR      = nan(C,C,A);

% ---------- per-animal loops ----------
for ai = 1:A
    % build X for this animal
    nCells = numel(all_cells{ai, cols(1)});
    X = false(nCells, C);
    for j = 1:C
        v = all_cells{ai, cols(j)};
        X(:,j) = logical(v(:));
        if numel(v) ~= nCells
            error('Inconsistent cell count within animal %d.', ai);
        end
    end
    N = nCells;

    for aIdx = 1:C
        Aset = X(:,aIdx);
        sumA = sum(Aset);
        for bIdx = 1:C
            Bset = X(:,bIdx);
            sumB = sum(Bset);

            k  = sum(Aset & Bset);
            ao = sum(Aset & ~Bset);
            bo = sum(~Aset & Bset);
            nn = sum(~Aset & ~Bset);

            % store counts
            k_all(aIdx,bIdx,ai) = k;
            a_all(aIdx,bIdx,ai) = ao;
            b_all(aIdx,bIdx,ai) = bo;
            n_all(aIdx,bIdx,ai) = nn;

            % symmetric overlaps
            denomJ = k + ao + bo;
            if denomJ==0
                J(aIdx,bIdx,ai) = NaN;
            else
                J(aIdx,bIdx,ai) = k / denomJ;
            end
            denomD = 2*k + ao + bo;
            if denomD==0
                DICE(aIdx,bIdx,ai) = NaN;
            else
                DICE(aIdx,bIdx,ai) = (2*k) / denomD;
            end

            % conjunction % of total cells
            CI_total(aIdx,bIdx,ai) = 100 * (k / N);

            % Fisher p (per requested tail)
            T2 = [k, ao; bo, nn];
            try
                switch tail
                    case "right"
                        [~, pF(aIdx,bIdx,ai)] = fishertest(T2,'Tail','right');
                    case "left"
                        [~, pF(aIdx,bIdx,ai)] = fishertest(T2,'Tail','left');
                    otherwise % both
                        [~, pF(aIdx,bIdx,ai)] = fishertest(T2,'Tail','both');
                end
            catch
                % fallback (rough): use hypergeometric for right/left; two-sided ~ 2*min
                K = sumA; M = sumB; % NOTE: margins vs notation in hygecdf
                % right: P(K>=k)
                p_right = 1 - hygecdf(k-1, N, K, M);
                p_left  = hygecdf(k,     N, K, M);
                switch tail
                    case "right", pF(aIdx,bIdx,ai) = p_right;
                    case "left",  pF(aIdx,bIdx,ai) = p_left;
                    otherwise,    pF(aIdx,bIdx,ai) = min(1, 2*min(p_left,p_right));
                end
            end

            % log-odds ratio (0.5 correction)
            LOR(aIdx,bIdx,ai) = log( ((k+0.5)*(nn+0.5)) / ((ao+0.5)*(bo+0.5)) );
        end
    end
end

% ---------- aggregate (equal-animal mean/SEM) ----------
S.per_animal.J        = J;
S.per_animal.DICE     = DICE;
S.per_animal.CI_total = CI_total;
S.per_animal.k        = k_all;
S.per_animal.a_only   = a_all;
S.per_animal.b_only   = b_all;
S.per_animal.n        = n_all;
S.per_animal.p_fisher = pF;
S.per_animal.LOR      = LOR;

S.mean.J        = nanmean(J,   3);
S.sem.J         = nanstd(J,[], 3) ./ sqrt(A);
S.mean.DICE     = nanmean(DICE,3);
S.sem.DICE      = nanstd(DICE,[],3) ./ sqrt(A);
S.mean.CI_total = nanmean(CI_total,3);
S.sem.CI_total  = nanstd(CI_total,[],3) ./ sqrt(A);

% ---------- Stouffer combined p ----------
S.stouffer.Z = nan(C,C);
S.stouffer.p = nan(C,C);
for aIdx = 1:C
    for bIdx = 1:C
        p_vec = squeeze(pF(aIdx,bIdx,:));
        w     = wS;
        ok = ~isnan(p_vec) & p_vec>0 & p_vec<1 & ~isnan(w);
        if ~any(ok), continue; end
        p_vec = p_vec(ok);  w = w(ok);

        switch tail
            case "right"
                z = norminv(1 - p_vec);                 % right-tailed
                Z = sum(w.*z) / sqrt(sum(w.^2));
                p_comb = 1 - normcdf(Z);
            case "left"
                z = norminv(p_vec);                     % left-tailed
                Z = sum(w.*z) / sqrt(sum(w.^2));
                p_comb = normcdf(Z);
            otherwise % both: sign by LOR, two-sided combine
                lor_vec = squeeze(LOR(aIdx,bIdx,ok));
                z_one   = norminv(1 - p_vec/2);         % two-sided -> one-sided z
                z_signed= sign(lor_vec) .* z_one;
                Z = sum(w.*z_signed) / sqrt(sum(w.^2));
                p_comb = 2*(1 - normcdf(abs(Z)));
        end
        S.stouffer.Z(aIdx,bIdx) = Z;
        S.stouffer.p(aIdx,bIdx) = p_comb;
    end
end

% ---------- Fixed-effect meta-analysis of log-OR ----------
S.meta.OR          = nan(C,C);
S.meta.CI          = nan(C,C,2);
S.meta.p_two_sided = nan(C,C);
S.meta.Q           = nan(C,C);
S.meta.I2          = nan(C,C);

for aIdx = 1:C
    for bIdx = 1:C
        a = squeeze(k_all(aIdx,bIdx,:));   % overlap
        b = squeeze(a_all(aIdx,bIdx,:));   % A only
        c = squeeze(b_all(aIdx,bIdx,:));   % B only
        d = squeeze(n_all(aIdx,bIdx,:));   % neither

        ok = ~isnan(a) & ~isnan(b) & ~isnan(c) & ~isnan(d);
        if ~any(ok), continue; end

        a=a(ok)+0.5; b=b(ok)+0.5; c=c(ok)+0.5; d=d(ok)+0.5;
        logOR_i = log((a.*d)./(b.*c));
        var_i   = 1./a + 1./b + 1./c + 1./d;
        w       = 1./var_i;

        if ~all(isfinite(logOR_i)) || ~all(isfinite(w)), continue; end

        logOR_hat = sum(w.*logOR_i) / sum(w);
        se_hat    = sqrt(1 / sum(w));
        z         = logOR_hat / se_hat;
        p2        = 2*(1 - normcdf(abs(z)));

        % heterogeneity
        Q  = sum(w .* (logOR_i - logOR_hat).^2);
        df = numel(logOR_i) - 1;
        I2 = max(0, (Q - df) / Q) * 100;

        S.meta.OR(aIdx,bIdx)          = exp(logOR_hat);
        S.meta.CI(aIdx,bIdx,1)        = exp(logOR_hat - 1.96*se_hat);
        S.meta.CI(aIdx,bIdx,2)        = exp(logOR_hat + 1.96*se_hat);
        S.meta.p_two_sided(aIdx,bIdx) = p2;
        S.meta.Q(aIdx,bIdx)           = Q;
        S.meta.I2(aIdx,bIdx)          = I2;
    end
end
end
