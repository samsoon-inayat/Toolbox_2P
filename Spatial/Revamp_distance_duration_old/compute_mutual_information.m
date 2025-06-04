function mi_value = compute_mutual_information(x, y, num_bins)
    % Computes mutual information (MI) between two variables x and y.
    % 
    % Inputs:
    %   x - First variable (e.g., firing rate), a 1D vector
    %   y - Second variable (e.g., speed or distance), a 1D vector
    %   num_bins - Number of bins for discretization (optional)
    % 
    % Output:
    %   mi_value - Mutual information (MI) in bits
    
    if nargin < 3
        num_bins = round(sqrt(length(x))); % Adaptive bin size
    end
    
    if isempty(y)
        y = 1:length(x);
    end

    % Remove NaNs
    valid_idx = ~isnan(x) & ~isnan(y);
    x = x(valid_idx);
    y = y(valid_idx);
    
    % Check for sufficient data
    if length(x) < 10 || length(y) < 10
        mi_value = NaN;
        return;
    end

    % Discretize variables into bins
    x_binned = discretize(x, num_bins);
    y_binned = discretize(y, num_bins);

    % Compute probability distributions
    epsilon = eps;%1e-10; % Small constant to prevent log(0)
    
    P_x = histcounts(x_binned, num_bins, 'Normalization', 'probability') + epsilon;
    P_y = histcounts(y_binned, num_bins, 'Normalization', 'probability') + epsilon;
    P_xy = histcounts2(x_binned, y_binned, num_bins, 'Normalization', 'probability') + epsilon;
    
    % Convert probabilities to column vectors for element-wise multiplication
    P_x = P_x(:);
    P_y = P_y(:);
    P_xy = P_xy(:);
    P_xy1 = P_x * P_y';

    % Compute mutual information
    mi_value = nansum(P_xy .* log2(P_xy ./ P_xy1(:)));

end
