function [dzMIt,dzMIi] = get_zMI_comp_dist_time(RsDt,RsTt,RsDi,RsTi)

ntrials = 50;
propsD = get_props_Rs(RsDt,ntrials);
propsT = get_props_Rs(RsTi,ntrials);
dzMIt = do_the_math_get(RsDt,RsTt,propsD,ntrials);
dzMIi = do_the_math_get(RsDi,RsTi,propsT,ntrials);



function dzMI = do_the_math_get(RsD,RsT,props,ntrials)
propsD = get_props_Rs(RsD,ntrials); propsT = get_props_Rs(RsT,ntrials);
dzMI = prop_op(propsD,propsT,0);
mean_diff = exec_fun_on_cell_mat(dzMI.diff_D_T,'nanmean');
mDiff = mean(mean_diff(:)); sDiff = std(mean_diff(:));
dzMI = prop_op(propsD,propsT,mDiff+sDiff);
dzMI.resp_D_g_T_and_good_FR = cell_list_op(dzMI.resp_D_g_T,props.good_FR,'and');
dzMI.resp_T_g_D_and_good_FR  = cell_list_op(dzMI.resp_T_g_D,props.good_FR,'and');
dzMI.resp_complex_and_good_FR = cell_list_op(dzMI.resp_complex,props.good_FR,'and');

