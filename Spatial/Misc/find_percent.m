function p = find_percent(popDT,numCells)

if ~exist('numCells','var')
    numCells = (exec_fun_on_cell_mat(popDT,'length'));
end
% p = 100*(exec_fun_on_cell_mat(popDT,'sum'))./(exec_fun_on_cell_mat(popDT,'length'));
p = 100*(exec_fun_on_cell_mat(popDT,'sum'))./numCells;

