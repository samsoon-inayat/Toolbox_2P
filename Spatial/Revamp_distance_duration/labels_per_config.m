function L = labels_per_config(X18)
% LABELS_PER_CONFIG  Map 0/1 membership to a single label per config.
% Input: X18 (cells×18) in order 6 configs × [T D S]
% Output: L (cells×6), labels: 1=T, 2=D, 3=S, 4=Multi, 5=None
nC = size(X18,1);
L  = 5*ones(nC,6);  % default None

for cfg = 1:6
    cols = (cfg-1)*3 + (1:3);   % [T D S] for this config
    M = X18(:,cols);            % cells × 3 logical
    n_on = sum(M,2);

    % exactly one type on -> that type (1..3)
    idx1 = find(n_on==1);
    for k = idx1.'
        L(k,cfg) = find(M(k,:), 1, 'first');
    end

    % two or three types on -> Multi (4)
    L(n_on>1, cfg) = 4;
end
end
