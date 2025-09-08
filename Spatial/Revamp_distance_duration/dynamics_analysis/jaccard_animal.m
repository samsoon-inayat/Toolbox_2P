function [D_cell, D_mean, D_sem, N_pair, P_same, P_same_sem] = jaccard_animal(idv, animal_idx)
% Jaccard-based identity distance for ONE animal.
% idv: A x C cell, each {nCells x 3} logical [T D S]
% Outputs:
%   D_cell : [nCells x C x C] per-cell distances (NaN for None<->None and self)
%   D_mean : [C x C] mean over cells (omit NaN), diag=NaN
%   D_sem  : [C x C] SEM over cells (omit NaN), diag=NaN
%   N_pair : [C x C] contributing cell counts per entry
%   P_same : [C x C] agreement map (= 1 - D_mean), diag=NaN
%   P_same_sem : [C x C] SEM of agreement (same as D_sem), diag=NaN

C = size(idv,2);
nCells = size(idv{animal_idx,1},1);

% Stack to cells x 3 x C
X3 = false(nCells,3,C);
for c = 1:C
    x = logical(idv{animal_idx,c});
    if size(x,1)~=nCells || size(x,2)~=3
        error('Animal %d: size mismatch at case %d', animal_idx, c);
    end
    X3(:,:,c) = x;
end

% Per-cell Jaccard distance
D_cell = nan(nCells, C, C);
for r = 1:nCells
    for i = 1:C
        vi = X3(r,:,i)>0;
        for j = 1:C
            if i==j
                D_cell(r,i,j) = NaN;                 % force self-comparison to NaN
            else
                vj = X3(r,:,j)>0;
                U = sum(vi | vj);
                if U==0
                    D_cell(r,i,j) = NaN;             % None<->None
                else
                    K = sum(vi & vj);
                    D_cell(r,i,j) = 1 - (K / U);     % 0 same … 1 full change
                end
            end
        end
    end
end

% Mean / SEM over cells, per (i,j)
D_mean = squeeze(mean(D_cell,1,'omitnan'));
N_pair = squeeze(sum(~isnan(D_cell),1));
D_sem  = squeeze(std(D_cell,0,1,'omitnan')) ./ sqrt(max(N_pair,1));

% Ensure diagonals are NaN (redundant but explicit)
D_mean(1:C+1:end) = NaN;
D_sem (1:C+1:end) = NaN;

% Agreement and its SEM
P_same     = 1 - D_mean;
P_same_sem = D_sem;                     % shift by 1 → same SEM
end
