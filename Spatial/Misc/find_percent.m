function p = find_percent(popDT)

p = 100*exec_fun_on_cell_mat(popDT,'sum')./exec_fun_on_cell_mat(popDT,'length');

