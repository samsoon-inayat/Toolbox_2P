%% first find the variables time_cells, distance_cells, and speed_cells in the final_time_dist_bin...m file
variable_combs = {'FR_time','FR_dist','FR_speed'};
metric = 'MI';
metric = 'PC';
clear time_cells distance_cells speed_cells
clear t_cells_T d_cells_T s_cells_T td_cells_T ds_cells_T ts_cells_T tds_cells_T ntds_cells_T
clear t_cells_I d_cells_I s_cells_I td_cells_I ds_cells_I ts_cells_I tds_cells_I ntds_cells_I
all_cells = {};
si = [Ar_t_T Ar_i_T Ar_t_D Ar_i_D ArL_t_T ArL_i_T ArL_t_D ArL_i_D Ars_t_T Ars_i_T Ars_t_D Ars_i_D]; propsPL = get_props_Rs(o.Rs(:,si),30);
si = [Ar_t_T Ar_i_T ArL_t_T ArL_i_T Ars_t_T Ars_i_T]; si_cn_ap = [[1 1 2 2 3 3];[1 2 1 2 1 2]];
% si = [Ar_t_D Ar_i_D ArL_t_D ArL_i_D Ars_t_D Ars_i_D]; 
props = get_props_Rs(o.Rs(:,si),30);
propsT = get_props_new(outT,met_valsT,props,si_cn_ap);
propsD = get_props_new(outD,met_valsD,props,si_cn_ap);
% [SinT,MixT,AllT] = get_the_pops(propsT,propsD);

%%
MVT = met_valsT; MVD = met_valsD;
for an = 1:5
    for cni = 1:size(si_cn_ap,2)
        cn = si_cn_ap(1,cni); ap = si_cn_ap(2,cni);
        idx = 3; pvalsPC = [MVT{1}{an,cn,ap}.PC(:,idx) MVD{2}{an,cn,ap}.PC(:,idx) MVT{3}{an,cn,ap}.PC(:,idx)];
        idx = 3; pvalsMI = [MVT{1}{an,cn,ap}.MI(:,idx) MVD{2}{an,cn,ap}.MI(:,idx) MVT{3}{an,cn,ap}.MI(:,idx)];
        idvPC{an,cni} = pvalsPC < 0.05;
        idvMI{an,cni} = pvalsMI < 0.05;
    end
end

idv = idvMI;

%% === INPUT: idv is 5x6 cell, each cell: [nCells x 3] logical [T D S] ===
% rows = animals (A=5), cols = 6 cases (e.g., C1A1, C1A2, C2A1, C2A2, C3A1, C3A2)

A = size(idv,1);
C = size(idv,2);         % should be 6
assert(C==6, 'Expected 6 cases (columns) in idv');

% Storage
D_mean_PC   = cell(A,1);   % 6x6 identity-distance (mean over cells)
P_same_PC   = cell(A,1);   % 6x6 agreement = 1 - distance
N_pair_PC   = cell(A,1);   % 6x6 contributing cell counts
Ptrans_PC   = cell(A,1);   % 5x5 global adjacent transitions (row-stochastic)
PtransS_PC  = cell(A,1);   % 5x5 transitions conditioned on start in T/D/S
Ctrans_PC   = cell(A,1);   % 5x5 integer counts

mean_adj_dist_PC = nan(A,1);   % mean distance for adjacent pairs per animal

% Labels for plotting
labs5 = {'T','D','S','Mixed','None'};
caseLabs = arrayfun(@(k)sprintf('Case %d',k), 1:C, 'UniformOutput', false);

for a = 1:A
    % --- Build X3: [nCells x 3 x 6] from your idv{a,c} ---
    nCells = size(idv{a,1},1);
    X3 = false(nCells,3,C);
    for c = 1:C
        X3(:,:,c) = logical(idv{a,c});
        if size(idv{a,c},1) ~= nCells || size(idv{a,c},2) ~= 3
            error('Animal %d: size mismatch at case %d', a, c);
        end
    end

    % --- Per-cell 6x6 Jaccard distance on the 3-bit vectors ---
    D_cell = nan(nCells, C, C);
    for r = 1:nCells
        for i = 1:C
            vi = X3(r,:,i)>0;   % 1x3
            for j = 1:C
                vj = X3(r,:,j)>0;
                U = sum(vi | vj);
                if U==0
                    D_cell(r,i,j) = NaN;   % None↔None: uninformative
                else
                    K = sum(vi & vj);
                    D_cell(r,i,j) = 1 - (K / U);
                end
            end
        end
    end
    Dm = squeeze(mean(D_cell,1,'omitnan'));      % 6x6
    D_mean_PC{a} = Dm;
    P_same_PC{a} = 1 - Dm;

    % Contributing counts per (i,j)
    Nij = zeros(C,C);
    for i=1:C
        for j=1:C
            Nij(i,j) = sum(~isnan(D_cell(:,i,j)));
        end
    end
    N_pair_PC{a} = Nij;

    % Adjacent-pair mean distance (i -> i+1)
    adj_idx = sub2ind([C C], 1:C-1, 2:C);
    mean_adj_dist_PC(a) = mean(Dm(adj_idx),'omitnan');

    % --- Global adjacent transitions among {T,D,S,Mixed,None} ---
    % Convert bits to labels 1..5 per cell×case
    L = 5*ones(nCells,C);  % default None=5
    for c = 1:C
        M = X3(:,:,c);                % nCells x 3
        s = sum(M,2);                 % #bits on
        idx1 = find(s==1);
        for k = idx1.'
            % find which of T(1),D(2),S(3) is on
            L(k,c) = find(M(k,:),1,'first');
        end
        L(s>=2, c) = 4;               % Mixed=4
        % None stays 5
    end

    % Aggregate transitions over the 5 adjacent steps
    C5  = zeros(5,5);
    C5s = zeros(5,5);    % conditioned on start in T/D/S (1..3)
    for t = 1:C-1
        from = L(:,t); to = L(:,t+1);
        ok = ~isnan(from) & ~isnan(to);
        for rlab = 1:5
            idx = ok & (from==rlab);
            if any(idx)
                for clab = 1:5
                    C5(rlab,clab) = C5(rlab,clab) + sum(to(idx)==clab);
                end
            end
        end
        % starts in T/D/S only
        okS = ok & ismember(from,1:3);
        for rlab = 1:5
            idx = okS & (from==rlab);
            if any(idx)
                for clab = 1:5
                    C5s(rlab,clab) = C5s(rlab,clab) + sum(to(idx)==clab);
                end
            end
        end
    end

    % Row-normalize to probabilities
    P5  = nan(5,5);
    P5s = nan(5,5);
    for rlab = 1:5
        den  = sum(C5(rlab,:));
        denS = sum(C5s(rlab,:));
        if den>0,  P5(rlab,:)  = C5(rlab,:) / den;   end
        if denS>0, P5s(rlab,:) = C5s(rlab,:) / denS; end
    end
    Ptrans_PC{a}  = P5;
    PtransS_PC{a} = P5s;
    Ctrans_PC{a}  = C5;
end

%% === Quick plots for Animal 1 ===
a = 3;
figure('Color','w');
subplot(1,3,1);
imagesc(D_mean_PC{a}, [0 1]); axis square; colormap(parula); colorbar
set(gca,'XTick',1:C,'XTickLabel',caseLabs,'YTick',1:C,'YTickLabel',caseLabs);
xtickangle(45); title(sprintf('Animal %d: Identity distance (PC)',a));

subplot(1,3,2);
imagesc(P_same_PC{a}, [0 1]); axis square; colormap(parula); colorbar
set(gca,'XTick',1:C,'XTickLabel',caseLabs,'YTick',1:C,'YTickLabel',caseLabs);
xtickangle(45); title('Agreement = 1 - distance');

subplot(1,3,3);
imagesc(PtransS_PC{a}, [0 1]); axis square; colormap(parula); colorbar
set(gca,'XTick',1:5,'XTickLabel',labs5,'YTick',1:5,'YTickLabel',labs5);
title('Adjacent transitions (start in T/D/S)');

%% === Across-animal summaries (PC) ===
fprintf('\nMean adjacent identity distance per animal (PC):\n');
disp(mean_adj_dist_PC');

% Group mean ± SEM
mAdj = mean(mean_adj_dist_PC,'omitnan');
sAdj = std(mean_adj_dist_PC,'omitnan') / sqrt(sum(~isnan(mean_adj_dist_PC)));
fprintf('Group mean adjacent distance (PC): %.3f ± %.3f (SEM)\n', mAdj, sAdj);


%%
%% --- Quick numbers & a permutation test (Animal 1, PC) ---
a = 1; C = 6; labs5 = {'T','D','S','Mixed','None'};

% 1) Mean adjacent identity distance (i -> i+1)
adj_idx = sub2ind([C C], 1:C-1, 2:C);
madj = mean(D_mean_PC{a}(adj_idx), 'omitnan');
fprintf('Animal %d (PC): mean adjacent distance = %.3f\n', a, madj);

% 2) Adjacent transition rows (conditioned on starting in T/D/S)
for r = 1:3
    v = PtransS_PC{a}(r,:);
    fprintf('%-5s ->  T:%4.2f  D:%4.2f  S:%4.2f  Mix:%4.2f  None:%4.2f\n', ...
        labs5{r}, v(1), v(2), v(3), v(4), v(5));
end

% 3) Permutation test: is stability above chance?
%    (Null preserves per-case label frequencies by shuffling cells within each case.)
%    For distance, "more stable than chance" means OBSERVED distance is LOWER than null.
X3 = false(size(idv{a,1},1),3,C);
for c = 1:C, X3(:,:,c) = logical(idv{a,c}); end

nperm = 1000;                         % reduce to 200 for a quick run
null_adj = nan(nperm,1);
for r = 1:nperm
    Xp = X3;
    for c = 1:C
        idx = randperm(size(X3,1));
        Xp(:,:,c) = X3(idx,:,c);      % shuffle across cells, preserve marginals
    end
    % compute mean adjacent distance for Xp
    Dp = nan(C,C);
    for i = 1:C
        for j = 1:C
            % Jaccard distance averaged over cells, with None<->None skipped
            U = sum( (Xp(:,:,i)|Xp(:,:,j)), 2 );
            K = sum( (Xp(:,:,i)&Xp(:,:,j)), 2 );
            dj = nan(size(U));
            nz = U>0;
            dj(nz) = 1 - (K(nz)./U(nz));
            Dp(i,j) = mean(dj, 'omitnan');
        end
    end
    null_adj(r) = mean(Dp(adj_idx), 'omitnan');
end
p_left  = mean(null_adj <= madj);                      % stability > chance (distance lower)
z_eff   = (madj - mean(null_adj)) / std(null_adj);     % effect size (negative = more stable)

fprintf('Permutation: p_left=%.4f (stability above chance), z=%.2f\n', p_left, z_eff);
%%
% Inputs: idv (A x 6 cell), choose animal a
a = 1; C = 6;
X3 = cell(C,1);
for c=1:C, X3{c} = logical(idv{a,c}); end  % [nCells x 3]

dists = nan(C-1,1);
for t = 1:C-1
    M1 = X3{t};      M2 = X3{t+1};
    tuned_both = (sum(M1,2)>0) & (sum(M2,2)>0);   % keep only both tuned
    if any(tuned_both)
        U = sum( M1(tuned_both,:) | M2(tuned_both,:), 2 );
        K = sum( M1(tuned_both,:) & M2(tuned_both,:), 2 );
        dj = 1 - (K./U);                          % Jaccard distance per cell
        dists(t) = mean(dj,'omitnan');
    end
end
mean_adj_tuned_only = mean(dists,'omitnan');
fprintf('Animal %d (PC): mean adj distance (tuned-only) = %.3f\n', a, mean_adj_tuned_only);

ESI = (mean(null_adj) - madj) / mean(null_adj);   % % reduction vs null
fprintf('Excess stability index: %.1f%% below null\n', 100*ESI);
