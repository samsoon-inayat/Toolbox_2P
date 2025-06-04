function [dzMIt] = get_dzMI_based_dist_time(RsDt,RsTt)

ntrials = 50;
dzMIt = do_the_math_get(RsDt,RsTt,[],ntrials);




function dzMI = do_the_math_get(RsD,RsT,props,ntrials)
propsD = get_props_Rs(RsD,ntrials); propsT = get_props_Rs(RsT,ntrials);
dzMI = prop_op(propsD,propsT,0);
% mean_diff = (exec_fun_on_cell_mat(dzMI.diff_D_T,'nanmean'));
% mDiff = mean(mean_diff(:)); sDiff = std(mean_diff(:));
% dzMI = prop_op(propsD,propsT,mDiff+sDiff);
% dzMI.resp_D_g_T_and_good_FR = cell_list_op(dzMI.resp_D_g_T,props.good_FR,'and');
% dzMI.resp_T_g_D_and_good_FR  = cell_list_op(dzMI.resp_T_g_D,props.good_FR,'and');
% dzMI.resp_complex_and_good_FR = cell_list_op(dzMI.resp_complex,props.good_FR,'and');
% dzMI.rs = prop_op_G(propsD,propsT,'rs');
% dzMI.HaFD = prop_op_G(propsD,propsT,'HaFD');
% dzMI.HiFD = prop_op_G(propsD,propsT,'HiFD');
% dzMI.mDiff = mDiff;
% dzMI.sDiff = sDiff;

