% Function to process each cell (neuron x trial matrix)
function mean_MI = process_cell(MI_matrix)
    MI_matrix(isnan(MI_matrix)) = 0;    % Replace NaNs with 0 (or use nanmean if ignoring)
    mean_MI = mean(MI_matrix, 2);       % Compute mean across trials for each neuron
end