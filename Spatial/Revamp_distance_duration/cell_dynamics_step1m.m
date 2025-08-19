%% STEP 1 — setup + sanity plot (Animal 1, MI)
% Requires: all_cellsnew_hm1 (5x36 cell), txl_all (1x36 cell array of names)
assert(exist('all_cellsnew_hm1','var')==1, 'all_cellsnew_hm1 not in workspace');
assert(exist('txl_all','var')==1, 'txl_all not in workspace');

tx = string(txl_all);

% ---- Build column indices for MI and PC (robust to order) ----
configs = ["C3-A1","C3-A2","C4-A1","C4-A2","C5-A1","C5-A2"];
types   = ["-T","-D","-S"];

mi18 = zeros(1,18); pc18 = zeros(1,18); k = 0;
for c = configs
    for t = types
        k = k + 1;
        mi18(k) = find(contains(tx,c) & contains(tx,"-MI-") & contains(tx,t), 1, 'first');
        pc18(k) = find(contains(tx,c) & contains(tx,"-PC-") & contains(tx,t), 1, 'first');
    end
end

% ---- Labels per cell for Animal 1, MI lens ----
a = 1;  % choose animal index (1..size(all_cellsnew_hm1,1))
X_MI = build_X(all_cellsnew_hm1, a, mi18);  % cells x 18 (6 configs x [T D S])
L_MI = labels_per_config(X_MI);             % cells x 6 labels: 1=T,2=D,3=S,4=Multi,5=None

% ---- Sanity plot: label counts per config (MI, Animal 1) ----
counts = zeros(6,5);   % rows=configs, cols=labels
for cfg = 1:6
    for lab = 1:5
        counts(cfg,lab) = sum(L_MI(:,cfg)==lab);
    end
end

figure('Color','w');
bar(counts,'stacked');
xticks(1:6); xticklabels(configs);
ylabel('# cells'); xlabel('Config (C-A)');
legend({'Time','Distance','Speed','Multi','None'}, 'Location','eastoutside');
title('Animal 1 — MI lens: cell labels per config');

disp('STEP 1 complete: MI/PC indices built, labels computed, sanity plot shown.');

%% STEP 2 — per-cell stability (Animal 1, MI) + plots
% Requires: build_X.m, labels_per_config.m in path
% Requires: all_cellsnew_hm1, txl_all in workspace (from STEP 1)

assert(exist('all_cellsnew_hm1','var')==1, 'all_cellsnew_hm1 not in workspace');
assert(exist('txl_all','var')==1, 'txl_all not in workspace');

% Rebuild MI indices if not present
if ~exist('mi18','var') || numel(mi18)~=18
    tx = string(txl_all);
    configs = ["C3-A1","C3-A2","C4-A1","C4-A2","C5-A1","C5-A2"];
    types   = ["-T","-D","-S"];
    mi18 = zeros(1,18); k = 0;
    for c = configs
        for t = types
            k = k + 1;
            mi18(k) = find(contains(tx,c) & contains(tx,"-MI-") & contains(tx,t), 1, 'first');
        end
    end
end

% Choose animal (1..size(all_cellsnew_hm1,1)); default to 1
if ~exist('a','var'), a = 1; end

% Build labels (cells x 6) for MI
X_MI = build_X(all_cellsnew_hm1, a, mi18);
L_MI = labels_per_config(X_MI);   % labels: 1=T,2=D,3=S,4=Multi,5=None

% ---------- Per-cell metrics ----------
C = size(L_MI,2);                             % 6 configs
adj_stability = mean(L_MI(:,1:C-1) == L_MI(:,2:C), 2, 'omitnan');  % 0..1 per cell
num_changes   = sum(L_MI(:,1:C-1) ~= L_MI(:,2:C), 2, 'omitnan');   % 0..5 per cell

% Optional extra summaries (print to console)
mean_adj = mean(adj_stability,'omitnan');
sem_adj  = std(adj_stability,'omitnan')/sqrt(size(L_MI,1));
frac_no_change = mean(num_changes==0,'omitnan');

fprintf('Animal %d (MI): mean adjacent stability = %.3f ± %.3f (SEM)\n', a, mean_adj, sem_adj);
fprintf('Animal %d (MI): fraction of cells with NO change = %.3f\n', a, frac_no_change);

% ---------- Plots ----------
figure('Color','w');

subplot(1,2,1);
histogram(adj_stability, 0:0.1:1);
xlabel('Adjacent stability (fraction staying same)'); ylabel('# cells');
title(sprintf('Animal %d — MI: per-cell stability', a));
xlim([0 1]); grid on; hold on
yl = ylim; plot([mean_adj mean_adj], yl, 'r--', 'LineWidth', 1);
text(mean_adj, 0.95*yl(2), sprintf(' mean=%.2f', mean_adj), 'Color','r', 'HorizontalAlignment','left');
hold off

subplot(1,2,2);
edges = -0.5:1:5.5;                             % bins for integer changes 0..5
histogram(num_changes, edges);
set(gca,'XTick',0:5);
xlabel('# label changes across 6 configs'); ylabel('# cells');
title(sprintf('Animal %d — MI: #changes', a));
grid on

disp('STEP 2 complete: per-cell metrics computed and plots shown.');

%%
%% STEP 3 — Across-animal summary of per-cell stability (MI vs PC)
% Requires: build_X.m, labels_per_config.m on path
% Requires: all_cellsnew_hm1, txl_all in workspace (from Step 1)

assert(exist('all_cellsnew_hm1','var')==1, 'all_cellsnew_hm1 not in workspace');
assert(exist('txl_all','var')==1, 'txl_all not in workspace');

A = size(all_cellsnew_hm1,1);         % # animals
tx = string(txl_all);

% Build MI/PC indices (robust to order)
configs = ["C3-A1","C3-A2","C4-A1","C4-A2","C5-A1","C5-A2"];
types   = ["-T","-D","-S"];

mi18 = zeros(1,18); pc18 = zeros(1,18); k = 0;
for c = configs
    for t = types
        k = k+1;
        mi18(k) = find(contains(tx,c) & contains(tx,"-MI-") & contains(tx,t), 1, 'first');
        pc18(k) = find(contains(tx,c) & contains(tx,"-PC-") & contains(tx,t), 1, 'first');
    end
end

% Compute mean adjacent stability per animal for MI and PC
mean_adj_MI = nan(A,1);
mean_adj_PC = nan(A,1);

for a = 1:A
    % --- MI ---
    X_MI = build_X(all_cellsnew_hm1, a, mi18);      % cells x 18
    L_MI = labels_per_config(X_MI);                 % cells x 6 (1..5)
    C = size(L_MI,2);                               % should be 6
    adj_MI = mean(L_MI(:,1:C-1) == L_MI(:,2:C), 2, 'omitnan');  % per-cell
    mean_adj_MI(a) = mean(adj_MI, 'omitnan');       % single number per animal

    % --- PC ---
    X_PC = build_X(all_cellsnew_hm1, a, pc18);
    L_PC = labels_per_config(X_PC);
    adj_PC = mean(L_PC(:,1:C-1) == L_PC(:,2:C), 2, 'omitnan');
    mean_adj_PC(a) = mean(adj_PC, 'omitnan');
end

% Basic paired stats across animals
[~, p_t]      = ttest(mean_adj_MI, mean_adj_PC);   % parametric
p_signrank    = signrank(mean_adj_MI, mean_adj_PC);% nonparametric (paired)

% Group means ± SEM
mMI  = mean(mean_adj_MI,'omitnan');
sMI  = std(mean_adj_MI,'omitnan')/sqrt(A);
mPC  = mean(mean_adj_PC,'omitnan');
sPC  = std(mean_adj_PC,'omitnan')/sqrt(A);

fprintf('Across animals:\n');
fprintf('  MI mean(adj stability)  = %.3f ± %.3f (SEM)\n', mMI, sMI);
fprintf('  PC mean(adj stability)  = %.3f ± %.3f (SEM)\n', mPC, sPC);
fprintf('  Paired t-test p = %.3g ; Wilcoxon signrank p = %.3g\n', p_t, p_signrank);

% ---------- Plot: paired dot/line per animal + group means ----------
figure('Color','w');

% Paired dots with connecting lines (one line per animal)
subplot(1,2,1);
hold on
xMI = ones(A,1)*1; xPC = ones(A,1)*2;
for a = 1:A
    plot([1 2], [mean_adj_MI(a) mean_adj_PC(a)], '-o', 'LineWidth', 1);
end
xlim([0.5 2.5]); ylim([0 1]);
set(gca,'XTick',[1 2],'XTickLabel',{'MI','PC'});
ylabel('Mean per-cell adjacent stability');
title('Per animal (paired)');

% Group means ± SEM on the right
subplot(1,2,2);
hold on
errorbar(1, mMI, sMI, 'o', 'LineWidth', 1.5); 
errorbar(2, mPC, sPC, 'o', 'LineWidth', 1.5);
xlim([0.5 2.5]); ylim([0 1]);
set(gca,'XTick',[1 2],'XTickLabel',{'MI','PC'});
ylabel('Mean per-cell adjacent stability');
title(sprintf('Group means (t p=%.3g, signrank p=%.3g)', p_t, p_signrank));
grid on

disp('STEP 3 complete: across-animal summary computed and plotted.');


%%
%% STEP 4 — Per-label (T/D/S) stability per animal: MI vs PC
% Requires: build_X.m, labels_per_config.m on path
% Requires: all_cellsnew_hm1, txl_all in workspace

assert(exist('all_cellsnew_hm1','var')==1, 'all_cellsnew_hm1 not in workspace');
assert(exist('txl_all','var')==1, 'txl_all not in workspace');

A = size(all_cellsnew_hm1,1);   % # animals
tx = string(txl_all);

% Build MI/PC indices (robust to order)
configs = ["C3-A1","C3-A2","C4-A1","C4-A2","C5-A1","C5-A2"];
types   = ["-T","-D","-S"];
mi18 = zeros(1,18); pc18 = zeros(1,18); k = 0;
for c = configs
    for t = types
        k = k+1;
        mi18(k) = find(contains(tx,c) & contains(tx,"-MI-") & contains(tx,t), 1, 'first');
        pc18(k) = find(contains(tx,c) & contains(tx,"-PC-") & contains(tx,t), 1, 'first');
    end
end

% Retention per label (1=T,2=D,3=S), per animal
ret_MI = nan(A,3);
ret_PC = nan(A,3);

for a = 1:A
    % --- MI ---
    L_MI = labels_per_config( build_X(all_cellsnew_hm1, a, mi18) );  % cells x 6
    % --- PC ---
    L_PC = labels_per_config( build_X(all_cellsnew_hm1, a, pc18) );

    % Compute retention for labels 1..3 (T/D/S): P(stay same label next step | current label)
    for lab = 1:3
        % MI
        starts_MI = sum(sum( L_MI(:,1:end-1)==lab ));
        stays_MI  = sum(sum( (L_MI(:,1:end-1)==lab) & (L_MI(:,2:end)==lab) ));
        ret_MI(a,lab) = (starts_MI>0) * (stays_MI / max(1,starts_MI));  % NaN if no starts
        if starts_MI==0, ret_MI(a,lab) = NaN; end

        % PC
        starts_PC = sum(sum( L_PC(:,1:end-1)==lab ));
        stays_PC  = sum(sum( (L_PC(:,1:end-1)==lab) & (L_PC(:,2:end)==lab) ));
        ret_PC(a,lab) = (starts_PC>0) * (stays_PC / max(1,starts_PC));
        if starts_PC==0, ret_PC(a,lab) = NaN; end
    end
end

% Group means ± SEM (ignore NaNs where a label never occurs)
mMI  = nanmean(ret_MI,1);  sMI  = nanstd(ret_MI,[],1)./sqrt(sum(~isnan(ret_MI),1));
mPC  = nanmean(ret_PC,1);  sPC  = nanstd(ret_PC,[],1)./sqrt(sum(~isnan(ret_PC),1));
labels = {'Time','Distance','Speed'};

% Paired stats per label (MI vs PC)
p_t   = nan(1,3);
p_wcx = nan(1,3);
for lab=1:3
    x = ret_MI(:,lab); y = ret_PC(:,lab);
    ok = ~isnan(x) & ~isnan(y);
    if sum(ok)>=2
        [~,p_t(lab)]   = ttest(x(ok), y(ok));
        p_wcx(lab)     = signrank(x(ok), y(ok));
    end
end

fprintf('Per-label retention (MI vs PC):\n');
for lab=1:3
    fprintf('  %-8s  MI=%.3f±%.3f  PC=%.3f±%.3f   t p=%.3g  signrank p=%.3g\n', ...
        labels{lab}, mMI(lab), sMI(lab), mPC(lab), sPC(lab), p_t(lab), p_wcx(lab));
end

% ---------- Plot 1: paired lines per animal for each label ----------
figure('Color','w'); 
for lab=1:3
    subplot(1,3,lab); hold on
    for a = 1:A
        if ~isnan(ret_MI(a,lab)) && ~isnan(ret_PC(a,lab))
            plot([1 2], [ret_MI(a,lab) ret_PC(a,lab)], '-o', 'LineWidth', 1);
        end
    end
    xlim([0.5 2.5]); ylim([0 1]);
    set(gca,'XTick',[1 2],'XTickLabel',{'MI','PC'});
    ylabel('Retention P(stay | current)');
    title(sprintf('%s  (t p=%.3g; sr p=%.3g)', labels{lab}, p_t(lab), p_wcx(lab)));
    grid on
end

% ---------- Plot 2: group means ± SEM for MI vs PC per label ----------
figure('Color','w'); hold on
x = 1:3;
dx = 0.18;
errorbar(x-dx, mMI, sMI, 'o-', 'LineWidth',1.5);
errorbar(x+dx, mPC, sPC, 'o-', 'LineWidth',1.5);
xlim([0.5 3.5]); ylim([0 1]);
set(gca,'XTick',1:3,'XTickLabel',labels);
ylabel('Retention P(stay | current)');
legend({'MI','PC'}, 'Location','eastoutside');
title('Per-label stability (mean ± SEM) across animals');
grid on

disp('STEP 4 complete: per-label retention computed and plotted.');

%%
%% STEP 5 — Animal-level stats: stability, per-label retention, transitions (MI & PC)
assert(exist('all_cellsnew_hm1','var')==1, 'all_cellsnew_hm1 not in workspace');
assert(exist('txl_all','var')==1, 'txl_all not in workspace');

% --- params ---
nperm = 2000;                   % permutations for null
labels3 = {'Time','Distance','Speed'};
lens_list = {'MI','PC'};        % analyze both lenses
rng default

% --- build column indices (robust to order) ---
tx = string(txl_all);
configs = ["C3-A1","C3-A2","C4-A1","C4-A2","C5-A1","C5-A2"];
types   = ["-T","-D","-S"];
mi18 = zeros(1,18); pc18 = zeros(1,18); k = 0;
for c = configs
    for t = types
        k = k+1;
        mi18(k) = find(contains(tx,c) & contains(tx,"-MI-") & contains(tx,t), 1, 'first');
        pc18(k) = find(contains(tx,c) & contains(tx,"-PC-") & contains(tx,t), 1, 'first');
    end
end
A = size(all_cellsnew_hm1,1);   % # animals

for a = 1:A
  fprintf('\n=== Animal %d ===\n', a);
  for LENS = 1:numel(lens_list)
    lens = lens_list{LENS};
    cols18 = strcmp(lens,'MI')*mi18 + strcmp(lens,'PC')*pc18;  % pick indices

    % ----- build labels: cells x 6 (values 1..5) -----
    X18 = build_X(all_cellsnew_hm1, a, cols18);
    L   = labels_per_config(X18);         % 1=T,2=D,3=S,4=Multi,5=None
    [nCells, C] = size(L);                % C should be 6

    % ===== 1) Overall adjacent stability vs chance =====
    obs_stab = mean(L(:,1:C-1)==L(:,2:C), 'all','omitnan');

    null_stab = zeros(nperm,1);
    for r=1:nperm
      Lp = L;
      for cfg=1:C
        Lp(:,cfg) = L(randperm(nCells), cfg);  % shuffle within each config (preserve marginals)
      end
      null_stab(r) = mean(Lp(:,1:C-1)==Lp(:,2:C), 'all','omitnan');
    end
    p_stab = mean(null_stab >= obs_stab);  % right-tailed

    % ===== 2) Per-label retention P(stay | current = T/D/S) =====
    obs_ret = nan(1,3);  % T,D,S
    for lab=1:3
      starts = sum(sum(L(:,1:C-1)==lab));
      stays  = sum(sum( (L(:,1:C-1)==lab) & (L(:,2:C)==lab) ));
      if starts>0
            obs_ret(lab) = stays/starts;
        else
            obs_ret(lab) = NaN;
        end
    end
    null_ret = nan(nperm,3);
    for r=1:nperm
      Lp = L;
      for cfg=1:C, Lp(:,cfg) = L(randperm(nCells), cfg); end
      for lab=1:3
        s = sum(sum(Lp(:,1:C-1)==lab));
        t = sum(sum((Lp(:,1:C-1)==lab) & (Lp(:,2:C)==lab)));
        if s>0
            null_ret(r,lab) = t/s;
        else
            null_ret(r,lab) = NaN;
        end
      end
    end
    p_ret = nan(1,3);
    for lab=1:3
      if ~isnan(obs_ret(lab))
        p_ret(lab) = mean(null_ret(:,lab) >= obs_ret(lab));  % right-tailed
      end
    end
    % FDR within these 3 tests (if mafdr available)
    if exist('mafdr','file')
      q_ret = nan(1,3);
      ok = ~isnan(p_ret);
      q_ret(ok) = mafdr(p_ret(ok),'BHFDR',true);
    else
      q_ret = p_ret;   % fallback: no adjustment
    end

    % ===== 3) T/D/S→T/D/S transition matrix & enrichment =====
    % Combine all adjacent steps across cells/configs, restrict to 1..3
    Pobs = nan(3,3);  % row-normalized transitions among T/D/S only
    counts = zeros(3,3);
    starts_vec = zeros(3,1);
    for lab=1:3
      starts_vec(lab) = sum(sum(L(:,1:C-1)==lab));
      for lab2=1:3
        counts(lab,lab2) = sum(sum((L(:,1:C-1)==lab) & (L(:,2:C)==lab2)));
      end
      if starts_vec(lab)>0
            Pobs(lab,:) = counts(lab,:)/starts_vec(lab);
        else
            Pobs(lab,:) = nan(1,3);
        end
    end
    % Null for each pair lab->lab2
    obs_pair = counts ./ max(1,starts_vec);
    null_pair = nan(nperm,3,3);
    for r=1:nperm
      Lp = L;
      for cfg=1:C, Lp(:,cfg) = L(randperm(nCells), cfg); end
      for lab=1:3
        s = sum(sum(Lp(:,1:C-1)==lab));
        for lab2=1:3
          t = sum(sum((Lp(:,1:C-1)==lab) & (Lp(:,2:C)==lab2)));
          if s>0
            null_pair(r,lab,lab2) = t/s;
        else
            null_pair(r,lab,lab2) = NaN;
        end
        end
      end
    end
    p_pair = nan(3,3);
    for lab=1:3
      for lab2=1:3
        v = squeeze(null_pair(:,lab,lab2));
        if ~isnan(obs_pair(lab,lab2))
          p_pair(lab,lab2) = mean(v >= obs_pair(lab,lab2));  % enrichment of that specific switch
        end
      end
    end
    % FDR across the 9 transitions (optional)
    if exist('mafdr','file')
      v = p_pair(:); ok = ~isnan(v);
      qtmp = nan(size(v)); qtmp(ok) = mafdr(v(ok),'BHFDR',true);
      q_pair = reshape(qtmp, size(p_pair));
    else
      q_pair = p_pair;
    end

    % ===== report =====
    fprintf('  [%s] overall stability = %.3f  (p=%.3g)\n', lens, obs_stab, p_stab);
    for lab=1:3
      fprintf('  [%s] retention %-8s = %.3f  (p=%.3g  q=%.3g)\n', ...
        lens, labels3{lab}, obs_ret(lab), p_ret(lab), q_ret(lab));
    end

    % ===== plots (per animal & lens) =====
    figure('Color','w','Name',sprintf('Animal %d — %s',a,lens));
    % (i) retention bars with stars
    subplot(1,2,1); hold on
    bar(1:3, obs_ret, 0.6);
    ylim([0 1]); xlim([0.5 3.5]);
    set(gca,'XTick',1:3,'XTickLabel',labels3); ylabel('P(stay | current)');
    title(sprintf('Animal %d — %s: retention (p_stab=%.3g)', a, lens, p_stab));
    % add stars for q<=0.05
    for lab=1:3
      if ~isnan(q_ret(lab)) && q_ret(lab)<=0.05
        text(lab, min(0.98, obs_ret(lab)+0.05), '*', 'HorizontalAlignment','center', 'FontSize',14);
      end
    end
    grid on

    % (ii) 3x3 transition heatmap (T/D/S only)
    subplot(1,2,2);
    imagesc(Pobs, [0 1]); axis square; colormap(parula); colorbar
    set(gca,'XTick',1:3,'XTickLabel',labels3,'YTick',1:3,'YTickLabel',labels3);
    title('Adj transitions: P(next | current)');
    % overlay dots for enriched transitions (q<=0.05)
    hold on
    [ri,ci] = find(~isnan(q_pair) & q_pair<=0.05);
    plot(ci,ri,'k.','MarkerSize',12);
    hold off
  end
end

disp('STEP 5 complete: animal-level stats computed and plotted for MI & PC.');
