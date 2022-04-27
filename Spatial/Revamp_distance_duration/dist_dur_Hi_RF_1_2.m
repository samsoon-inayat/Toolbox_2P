function dist_dur_Hi_RF_1


%% check the difference in zMI for Dis and Dur for the different cell types DurC, DisC, and DDM
while 1
    %%
    mean_dzMI = [];
    FD_Prop = dzMI_FD.diff_T_D; FT_Prop = dzMI_FT.diff_T_D; 
%     FD_Prop = dzMI_FD.rs.diff_T_D; FT_Prop = dzMI_FT.rs.diff_T_D; 
%     FD_Prop = dzMI_FD.HaFD.diff_T_D; FT_Prop = dzMI_FT.HaFD.diff_T_D; 
%     FD_Prop = dzMI_FD.HiFD.diff_T_D; FT_Prop = dzMI_FT.HiFD.diff_T_D; 
%     FD_Prop = props{rfi,1}.N_Resp_Trials; FT_Prop = props{rfi,3}.N_Resp_Trials;
%     FD_Prop = props{1}.mean_FR; FT_Prop = props{3}.mean_FR;
    for rfi = 4
        TD = FD_Prop;
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dur_cells_T)];
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dis_cells_T)];

        TD = FT_Prop;
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dur_cells_I)];
        mean_dzMI = [ mean_dzMI exec_fun_on_cell_mat(TD,'nanmean',dis_cells_I)];
    end
    [within,dvn,xlabels] = make_within_table({'TI','CT','Cond'},[2,2,3]);
    dataT = make_between_table({mean_dzMI},dvn);
    ra = RMA(dataT,within,{'hsd'});
    ra.ranova
    print_for_manuscript(ra)
    %%
    break;
end
