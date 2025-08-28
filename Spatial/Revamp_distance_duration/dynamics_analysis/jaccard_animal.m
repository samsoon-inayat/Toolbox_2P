function [D_cell, D_mean, D_sem, N_pair, P_same] = jaccard_animal(idv, animal_idx)
% Jaccard-based identity distance for ONE animal.
% idv: A x 6 cell, each {nCells x 3} logical [T D S]
% Outputs:
%   D_cell : [nCells x 6 x 6] per-cell distances (NaN for None<->None)
%   D_mean : [6 x 6] mean over cells (omit NaN)
%   D_sem  : [6 x 6] SEM over cells (omit NaN)
%   N_pair : [6 x 6] contributing cell counts per entry
%   P_same : [6 x 6] agreement map (= 1 - D_mean)

C = size(idv,2);  assert(C==6,'Expected 6 cases');
nCells = size(idv{animal_idx,1},1);

% Stack to cells x 3 x 6
X3 = false(nCells,3,C);
for c = 1:C
    x = logical(idv{animal_idx,c});
    if size(x,1)~=nCells || size(x,2)~=3
        error('Animal %d: size mismatch at case %d', animal_idx, c);
    end
    X3(:,:,c) = x;
end

% Per-cell 6x6 Jaccard distance
D_cell = nan(nCells, C, C);
for r = 1:nCells
    for i = 1:C
        vi = X3(r,:,i)>0;
        for j = 1:C
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

% Mean / SEM over cells, per (i,j)
D_mean = squeeze(mean(D_cell,1,'omitnan'));
N_pair = squeeze(sum(~isnan(D_cell),1));
D_sem  = squeeze(std(D_cell,0,1,'omitnan')) ./ sqrt(max(N_pair,1));
P_same = 1 - D_mean;
end
