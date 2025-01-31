function I = mutualinfo(X, Y)
    % Compute joint histogram
    joint_hist = histcounts2(X, Y, 'Normalization', 'probability');
    
    % Marginal distributions
    px = sum(joint_hist, 2);  % Sum across rows (marginal for X)
    py = sum(joint_hist, 1);  % Sum across columns (marginal for Y)
    
    % Joint entropy
    joint_entropy = -nansum(joint_hist(:) .* log2(joint_hist(:) + eps));
    
    % Marginal entropies
    hx = -nansum(px .* log2(px + eps));
    hy = -nansum(py .* log2(py + eps));
    
    % Mutual information
    I = hx + hy - joint_entropy;
end
