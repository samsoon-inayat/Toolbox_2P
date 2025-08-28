function P = stability_percentages(idv)
% Compute Stable/Dynamic/Middle percentages per label for each animal.
% Returns struct P with fields:
%   .stable   : A x 3 percentages (T,D,S)
%   .dynamic  : A x 3
%   .middle   : A x 3
%   .denom    : A x 3 (#cells that ever had that label)
%   .mean_sem : struct with group means and SEMs for each category

A = size(idv,1);
stable = nan(A,3); dynamic = nan(A,3); middle = nan(A,3); denom = nan(A,3);

for a = 1:A
    S_bit = stability_per_label(idv,a);
    [is_stable, is_dynamic] = stability_classify(S_bit);
    for b = 1:3
        mask = ~isnan(S_bit(:,b));            % cells that ever had this label
        d = sum(mask);
        denom(a,b) = d;
        if d==0, continue; end
        nS = sum(is_stable(mask,b));
        nD = sum(is_dynamic(mask,b));
        nM = d - nS - nD;
        stable(a,b)  = 100*nS/d;
        dynamic(a,b) = 100*nD/d;
        middle(a,b)  = 100*nM/d;
    end
end

% Group means/SEMs
ms = mean(stable, 1, 'omitnan');   ss = std(stable, 0, 1, 'omitnan')./sqrt(sum(~isnan(stable),1));
md = mean(dynamic,1, 'omitnan');   sd = std(dynamic,0, 1, 'omitnan')./sqrt(sum(~isnan(dynamic),1));
mm = mean(middle, 1, 'omitnan');   sm = std(middle, 0, 1, 'omitnan')./sqrt(sum(~isnan(middle),1));

P.stable   = stable;
P.dynamic  = dynamic;
P.middle   = middle;
P.denom    = denom;
P.mean_sem.stable_mean  = ms;  P.mean_sem.stable_sem  = ss;
P.mean_sem.dynamic_mean = md;  P.mean_sem.dynamic_sem = sd;
P.mean_sem.middle_mean  = mm;  P.mean_sem.middle_sem  = sm;
end
