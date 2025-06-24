function props = get_the_pops(propsT,propsD)

props.cells_time = propsT.newPC.cells_time;
props.cells_dist = propsD.newPC.cells_dist;
props.cells_speed = propsT.newPC.cells_speed;




ctnums_names = {'NR','speed','dist','dist_speed','time','time_speed','time_dist','time_dist_speed'};

tProps = props;
singularly_tuned = cell_list_op(tProps.newPC.cells_time,tProps.newPC.cells_dist,'or');
singularly_tuned = cell_list_op(singularly_tuned,tProps.newPC.cells_speed,'or');
mixed_tuned = cell_list_op(tProps.newPC.cells_time_dist,tProps.newPC.cells_dist_speed,'or');
mixed_tuned = cell_list_op(mixed_tuned,tProps.newPC.cells_time_speed,'or');
all_tuned = tProps.newPC.cells_time_dist_speed;

props.newPC.cells_singularly_tuned = singularly_tuned;
props.newPC.cells_mixed_tuned = mixed_tuned;
props.newPC.cells_all_tuned = all_tuned;


singularly_tuned = cell_list_op(tProps.newMI.cells_time,tProps.newMI.cells_dist,'or');
singularly_tuned = cell_list_op(singularly_tuned,tProps.newMI.cells_speed,'or');
mixed_tuned = cell_list_op(tProps.newMI.cells_time_dist,tProps.newMI.cells_dist_speed,'or');
mixed_tuned = cell_list_op(mixed_tuned,tProps.newMI.cells_time_speed,'or');
all_tuned = tProps.newMI.cells_time_dist_speed;

props.newMI.cells_singularly_tuned = singularly_tuned;
props.newMI.cells_mixed_tuned = mixed_tuned;
props.newMI.cells_all_tuned = all_tuned;

props.newPC.names = {'NR','speed','dist','dist_speed','time','time_speed','time_dist','time_dist_speed','singularly_tuned','mixed_tuned','all_tuned'};
props.newMI.names = {'NR','speed','dist','dist_speed','time','time_speed','time_dist','time_dist_speed','singularly_tuned','mixed_tuned','all_tuned'};