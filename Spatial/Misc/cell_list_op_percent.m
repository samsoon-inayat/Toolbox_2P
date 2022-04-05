function presp = cell_list_op_percent(resp,ccs)

if exist('ccs','var')
    presp = 100*exec_fun_on_cell_mat(resp,'sum',ccs)./exec_fun_on_cell_mat(resp,'length');
else
    presp = 100*exec_fun_on_cell_mat(resp,'sum')./exec_fun_on_cell_mat(resp,'length');
end